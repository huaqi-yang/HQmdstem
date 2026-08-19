#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract the length distribution of white atomic chains from a TEM image.

The script automatically:
  1. locates the dark TEM panel inside the white figure border;
  2. detects the black scale bar in the bottom-left label box;
  3. thresholds the panel to segment white chains;
  4. keeps every bright component by default (set --min-ratio higher to drop
     round single-atom blobs);
  5. measures each chain along its medial path (geodesic length in px),
     converts it to nm and writes:
       - <prefix>_stats.csv
       - <prefix>_distribution.png    (bar chart)
       - <prefix>_overlay.png  (chains marked in the original image)

Only numpy and Pillow are required. If matplotlib is installed it is used for
the bar chart, otherwise a Pillow-drawn chart is produced.

Example:
    python extract_chain_lengths.py "D:/gpumd/2picture/CT2/bondangelGr/wenzhang/figures/图片1.png"
"""

import argparse
import csv
import heapq
from collections import deque
from pathlib import Path

import numpy as np

try:
    from PIL import Image, ImageDraw, ImageFont

    PIL_AVAILABLE = True
except Exception:
    PIL_AVAILABLE = False


def load_gray(path):
    """Load an image as grayscale uint8 using Pillow or matplotlib."""
    if PIL_AVAILABLE:
        return np.asarray(Image.open(path).convert("L"), dtype=np.uint8)
    import matplotlib.image as mpimg

    rgb = np.asarray(mpimg.imread(str(path)))
    if rgb.ndim == 3:
        rgb = rgb[:, :, :3]
    gray = rgb @ np.array([0.2989, 0.5870, 0.1140])
    return np.clip(gray * 255.0, 0, 255).astype(np.uint8)


# ---------------------------------------------------------------------------
# image helpers
# ---------------------------------------------------------------------------

def longest_run(values, threshold=245, less=True):
    """Return (length, start, end) of the longest run satisfying the predicate."""
    mask = values < threshold if less else values >= threshold
    best_len, best_start, best_end = 0, None, None
    i, n = 0, len(mask)
    while i < n:
        if mask[i]:
            j = i
            while j < n and mask[j]:
                j += 1
            if j - i > best_len:
                best_len, best_start, best_end = j - i, i, j - 1
            i = j
        else:
            i += 1
    return best_len, best_start, best_end


def detect_panel(gray, margin_threshold=245):
    """Return (y0, y1, x0, x1) of the dark TEM panel inside white margins."""
    row_mean = gray.mean(axis=1)
    col_mean = gray.mean(axis=0)
    _, y0, y1 = longest_run(row_mean, margin_threshold)
    _, x0, x1 = longest_run(col_mean, margin_threshold)
    if y0 is None or x0 is None:
        return 0, gray.shape[0], 0, gray.shape[1]
    return y0, y1 + 1, x0, x1 + 1


def detect_scale_bar(panel, min_len=100, max_len=160, bright_frac=0.5):
    """Detect the black scale bar inside the bottom-left white label box.

    Returns (row, x0, x1, length_px) in panel coordinates, or None.
    """
    h, w = panel.shape
    top = max(0, h - 170)
    region = panel[top : max(0, h - 20), : min(300, w)]
    dark = region < 150
    candidates = []
    for yy in range(dark.shape[0]):
        row = dark[yy]
        i, n = 0, len(row)
        while i < n:
            if row[i]:
                j = i
                while j < n and row[j]:
                    j += 1
                length = j - i
                if min_len <= length <= max_len:
                    candidates.append((length, yy + top, i, j - 1))
                i = j
            else:
                i += 1

    candidates.sort(key=lambda c: -c[0])
    for length, sy, xa, xb in candidates:
        x0b, x1b = max(0, xa - 20), min(w, xb + 20)
        y0a, y1a = max(0, sy - 58), max(0, sy - 28)
        y0b, y1b = min(h, sy + 2), min(h, sy + 14)
        if y1a <= y0a or y1b <= y0b or x1b <= x0b:
            continue
        frac_above = float((panel[y0a:y1a, x0b:x1b] > 210).mean())
        frac_below = float((panel[y0b:y1b, x0b:x1b] > 210).mean())
        if frac_above > bright_frac and frac_below > bright_frac:
            return sy, xa, xb, length
    return None


def label_components(mask):
    """8-connected component labelling; returns a list of (y, x) pixel arrays."""
    h, w = mask.shape
    label = np.zeros((h, w), dtype=np.int32)
    components = []
    current = 0
    for i in range(h):
        for j in range(w):
            if mask[i, j] and label[i, j] == 0:
                current += 1
                queue = deque([(i, j)])
                label[i, j] = current
                pixels = [(i, j)]
                while queue:
                    r, c = queue.popleft()
                    for dr in (-1, 0, 1):
                        for dc in (-1, 0, 1):
                            if dr == 0 and dc == 0:
                                continue
                            nr, nc = r + dr, c + dc
                            if (
                                0 <= nr < h
                                and 0 <= nc < w
                                and mask[nr, nc]
                                and label[nr, nc] == 0
                            ):
                                label[nr, nc] = current
                                queue.append((nr, nc))
                                pixels.append((nr, nc))
                components.append(np.array(pixels, dtype=np.int64))
    return components


def geodesic_length(points):
    """Longest shortest path inside a component (8-connected, diagonal=sqrt2)."""
    n = len(points)
    if n <= 1:
        return 0.0
    ys, xs = points[:, 0], points[:, 1]
    y0, y1 = int(ys.min()), int(ys.max())
    x0, x1 = int(xs.min()), int(xs.max())
    lh, lw = y1 - y0 + 1, x1 - x0 + 1
    local = np.full((lh, lw), -1, dtype=np.int32)
    for k in range(n):
        local[int(ys[k]) - y0, int(xs[k]) - x0] = k

    adj = [[] for _ in range(n)]
    for k in range(n):
        y, x = int(ys[k]), int(xs[k])
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                if dr == 0 and dc == 0:
                    continue
                yy, xx = y + dr - y0, x + dc - x0
                if 0 <= yy < lh and 0 <= xx < lw:
                    j = int(local[yy, xx])
                    if j >= 0 and j != k:
                        weight = 1.0 if dr == 0 or dc == 0 else 2.0 ** 0.5
                        adj[k].append((j, weight))

    def dijkstra(start):
        dist = [1e18] * n
        dist[start] = 0.0
        heap = [(0.0, start)]
        while heap:
            d, u = heapq.heappop(heap)
            if d > dist[u]:
                continue
            for v, weight in adj[u]:
                nd = d + weight
                if nd < dist[v]:
                    dist[v] = nd
                    heapq.heappush(heap, (nd, v))
        farthest = max(range(n), key=lambda i: dist[i])
        return farthest, dist

    far1, _ = dijkstra(0)
    far2, dist2 = dijkstra(far1)
    return float(dist2[far2])


def local_background(gray, radius):
    """Local background estimate via a square box average (numpy only)."""
    if radius <= 0:
        return np.zeros_like(gray, dtype=np.float64)
    img = gray.astype(np.float64)
    k = 2 * radius + 1
    padded = np.pad(img, ((radius, radius), (radius, radius)), mode="edge")
    cum = np.cumsum(np.cumsum(padded, axis=0), axis=1)
    cum = np.pad(cum, ((1, 0), (1, 0)))
    h, w = img.shape
    rows = np.arange(h, dtype=np.int64)
    cols = np.arange(w, dtype=np.int64)
    total = (
        cum[rows[:, None] + k, cols[None, :] + k]
        - cum[rows[:, None], cols[None, :] + k]
        - cum[rows[:, None] + k, cols[None, :]]
        + cum[rows[:, None], cols[None, :]]
    )
    return total / float(k * k)


def load_font(size):
    """Load a readable TTF font on Windows/Linux, fall back to the default."""
    if not PIL_AVAILABLE:
        return None
    candidates = [
        r"C:\Windows\Fonts\arial.ttf",
        r"C:\Windows\Fonts\segoeui.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]
    for path in candidates:
        if Path(path).exists():
            try:
                return ImageFont.truetype(path, size)
            except Exception:
                pass
    return ImageFont.load_default()


# ---------------------------------------------------------------------------
# bar chart
# ---------------------------------------------------------------------------

def _histogram(lengths, bin_width):
    lengths = np.asarray(lengths, dtype=float)
    if lengths.size == 0:
        return np.zeros(1, dtype=int), np.array([0.0, bin_width])
    edges = np.arange(0.0, float(lengths.max()) + bin_width, bin_width)
    if len(edges) < 2:
        edges = np.array([0.0, bin_width])
    counts, edges = np.histogram(lengths, bins=edges)
    return counts, edges


CHART_STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 11,
    "axes.linewidth": 1.3,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "xtick.major.size": 3.5,
    "xtick.major.width": 1.3,
    "xtick.direction": "in",
    "ytick.major.size": 3.5,
    "ytick.major.width": 1.3,
    "ytick.direction": "in",
    "lines.linewidth": 1.9,
    "legend.fontsize": 9,
    "legend.frameon": False,
    "figure.dpi": 300,
    "savefig.dpi": 2300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
}


def _chart_with_matplotlib(lengths, bin_width, out_path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    counts, edges = _histogram(lengths, bin_width)
    centers = (edges[:-1] + edges[1:]) / 2.0

    with plt.rc_context(CHART_STYLE):
        fig, ax = plt.subplots(figsize=(8 / 2.54, 5 / 2.54))
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")
        fig.subplots_adjust(left=0.10, bottom=0.12, right=0.99, top=0.99)

        from matplotlib.colors import LinearSegmentedColormap

        cmap = LinearSegmentedColormap.from_list(
            "chain_palette", ["#81C784", "#FFEB3B", "#FF9800"])
        bar_colors = [cmap(i / max(len(centers) - 1, 1)) for i in range(len(centers))]
        ax.bar(centers, counts, width=bin_width * 0.92, color=bar_colors,
               edgecolor="#33691e", linewidth=0.5, zorder=2)
        if counts.size:
            xs = np.linspace(edges[0], edges[-1], 240)
            ys = np.interp(xs, centers, counts.astype(float))
            ax.plot(xs, ys, color="#E65100", lw=1.6,
                    solid_capstyle="butt", clip_on=False, zorder=3)
            ax.plot(centers, counts, color="#E65100", lw=0,
                    marker="o", markersize=3.2, clip_on=False, zorder=4)

        ax.set_xlabel("Chain length (nm)", labelpad=7)
        ax.set_ylabel("Number of chains", labelpad=7)
        ax.tick_params(which="both", length=4, width=1.2,
                       direction="in", right=True, top=True)
        ax.tick_params(which="minor", length=4)
        ax.grid(False)
        if counts.size:
            from matplotlib.ticker import MaxNLocator, FormatStrFormatter
            ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
            ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
            ax.set_xlim(edges[0] - 0.02, edges[-1] + 0.02)
            ax.set_ylim(0, max(float(counts.max()) * 1.15, 1.0))
            for x, y in zip(centers, counts):
                if y:
                    ax.text(x, y, str(int(y)), ha="center", va="bottom",
                            fontsize=9, color="#212121")

        base = str(out_path)
        if base.lower().endswith(".png"):
            base = base[:-4]
        for fmt in ("png", "pdf", "svg"):
            try:
                fig.savefig(f"{base}.{fmt}", format=fmt, facecolor="white")
            except Exception as exc:
                print(f"[WARN] could not save {base}.{fmt}: {exc}")
        plt.close(fig)


def _chart_with_pil(lengths, bin_width, out_path, title):
    if not PIL_AVAILABLE:
        raise ImportError("Pillow is required for the fallback bar chart")
    counts, edges = _histogram(lengths, bin_width)
    centers = (edges[:-1] + edges[1:]) / 2.0
    n = len(counts)
    max_count = max(int(counts.max()), 1)

    width, height = 900, 560
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    font_title = load_font(22)
    font_label = load_font(17)
    font_tick = load_font(14)

    left, right = 60, 920
    top, bottom = 60, 505
    plot_w = right - left
    plot_h = bottom - top
    bar_w = plot_w / n

    if title:
        draw.text(((left + right) / 2, 38), title, fill="black",
                  font=font_title, anchor="mm")

    tick_step = max(1, int(np.ceil(max_count / 6)))
    yticks = list(range(0, max_count + 1, tick_step))
    for yv in yticks:
        y = bottom - (yv / max_count) * plot_h
        draw.line([left, y, right, y], fill=(225, 225, 225), width=1)
        draw.text((left - 8, y), str(yv), fill="black", font=font_tick, anchor="rm")

    def _grad(c1, c2, t):
        return tuple(int(a + (b - a) * t) for a, b in zip(c1, c2))

    _c_green, _c_yellow, _c_orange = (129, 199, 132), (255, 235, 59), (255, 152, 0)
    def _bar_color(i):
        t = i / max(n - 1, 1)
        if t < 0.5:
            return _grad(_c_green, _c_yellow, t * 2)
        return _grad(_c_yellow, _c_orange, (t - 0.5) * 2)

    for i in range(n):
        if counts[i] == 0:
            continue
        x0 = left + i * bar_w + bar_w * 0.12
        x1 = left + (i + 1) * bar_w - bar_w * 0.12
        y0 = bottom - (counts[i] / max_count) * plot_h
        draw.rectangle([x0, y0, x1, bottom], fill=_bar_color(i),
                       outline="black", width=1)
        draw.text((x0 + (x1 - x0) / 2, y0 - 16), str(int(counts[i])),
                  fill="black", font=font_tick, anchor="mm")

    if n > 1:
        cx = np.array([left + i * bar_w + bar_w / 2 for i in range(n)], dtype=float)
        cy = bottom - (counts / max_count) * plot_h
        sx = np.linspace(left + bar_w / 2, right - bar_w / 2, 240)
        sy = np.interp(sx, cx, cy)
        draw.line([(int(x), int(y)) for x, y in zip(sx, sy)],
                  fill=(230, 81, 0), width=2, joint="curve")

    tick_step = max(1, int(np.ceil(n / 6)))
    for i in range(0, n, tick_step):
        x = left + i * bar_w + bar_w / 2
        draw.text((x, bottom + 8), f"{centers[i]:.2f}", fill="black",
                  font=font_tick, anchor="ma")

    draw.line([left, top, left, bottom], fill="black", width=2)
    draw.line([left, bottom, right, bottom], fill="black", width=2)
    draw.text(((left + right) / 2, bottom + 38), "Chain length (nm)",
              fill="black", font=font_label, anchor="mm")

    # rotated y-axis label
    label_img = Image.new("RGBA", (260, 34), (255, 255, 255, 0))
    label_draw = ImageDraw.Draw(label_img)
    label_draw.text((130, 17), "Number of chains", fill="black",
                    font=font_label, anchor="mm")
    label_img = label_img.rotate(90, expand=True)
    img.paste(label_img, (20, (top + bottom) // 2 - label_img.height // 2),
              label_img)

    if lengths is None or len(lengths) == 0:
        draw.text(((left + right) / 2, (top + bottom) / 2),
                  "No chains detected", fill="gray", font=font_title, anchor="mm")

    img.save(out_path)


def draw_chart(lengths, bin_width, out_path, title=None):
    try:
        _chart_with_matplotlib(lengths, bin_width, str(out_path))
        return "matplotlib"
    except Exception:
        _chart_with_pil(lengths, bin_width, str(out_path), title)
        return "Pillow"


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Measure white atomic chain lengths in a TEM image."
    )
    parser.add_argument("input", help="input image path (PNG/JPG/TIFF...)")
    parser.add_argument("--output-dir", default=None,
                        help="output directory (default: beside the input image)")
    parser.add_argument("--prefix", default="chain_length",
                        help="output file name prefix (default: chain_length)")
    parser.add_argument("--threshold", type=int, default=150,
                        help="absolute minimum intensity for chain pixels (default: 150)")
    parser.add_argument("--adaptive-radius", type=int, default=35,
                        help="local background blur radius in px; 0 disables adaptive threshold (default: 35)")
    parser.add_argument("--adaptive-offset", type=int, default=10,
                        help="required brightness above the local background (default: 10)")
    parser.add_argument("--supplement-floor", type=int, default=105,
                        help="low threshold for catching small dim atoms (default: 105)")
    parser.add_argument("--supplement-offset", type=int, default=5,
                        help="low adaptive offset for small dim atoms (default: 5)")
    parser.add_argument("--supplement-max-area", type=int, default=150,
                        help="max area of supplemented small components (default: 150)")
    parser.add_argument("--min-area", type=int, default=8,
                        help="minimum chain area in pixels (default: 8)")
    parser.add_argument("--min-ratio", type=float, default=1.0,
                        help="minimum length/width ratio to keep a chain "
                             "(default 1.0: keep every shape; increase to 2 to drop round atoms)")
    parser.add_argument("--bin-width", type=float, default=0.1,
                        help="bar chart bin width in nm (default: 0.1)")
    parser.add_argument("--scale-nm", type=float, default=0.5,
                        help="length of the detected scale bar in nm (default: 0.5)")
    parser.add_argument("--scale-px", type=float, default=None,
                        help="override the auto-detected scale bar length in pixels")
    parser.add_argument("--nm-per-px", type=float, default=None,
                        help="override the pixel calibration in nm/pixel")
    parser.add_argument("--margin", type=int, default=5,
                        help="pixels eroded from the panel edges (default: 5)")
    parser.add_argument("--exclude-rect", default=None,
                        help="extra region to exclude: y0,y1,x0,x1 in panel coords")
    args = parser.parse_args()

    input_path = Path(args.input)
    prefix = args.prefix.replace(' ', '_')
    if not input_path.exists():
        raise SystemExit(f"input image not found: {input_path}")

    output_dir = Path(args.output_dir) if args.output_dir else input_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- load image and locate the dark TEM panel ----
    gray = load_gray(input_path)
    py0, py1, px0, px1 = detect_panel(gray)
    panel = gray[py0:py1, px0:px1].copy()
    ph, pw = panel.shape
    print(f"panel: rows {py0}-{py1 - 1}, cols {px0}-{px1 - 1} "
          f"({pw} x {ph} px)")

    # ---- scale bar calibration ----
    scale_pos = None
    scale_px = None
    scale = detect_scale_bar(panel)
    if args.scale_px is not None:
        scale_px = float(args.scale_px)
    elif scale is not None:
        scale_pos = scale[:3]
        scale_px = float(scale[3])
    else:
        print("[WARN] scale bar not detected; using --scale-px fallback of 126 px")
        scale_px = float(args.scale_px) if args.scale_px is not None else 126.0

    if args.nm_per_px is not None:
        nm_per_px = float(args.nm_per_px)
    else:
        nm_per_px = args.scale_nm / scale_px
    print(f"scale bar: {scale_px:g} px = {args.scale_nm:g} nm "
          f"(nm/px = {nm_per_px:.6g})")

    # ---- build the chain mask ----
    if args.adaptive_radius > 0:
        local_bg = local_background(panel, args.adaptive_radius)
        mask = ((panel.astype(np.int16) - local_bg) > args.adaptive_offset) & (panel > args.threshold)
    else:
        mask = panel > args.threshold
    margin = args.margin
    mask[:margin, :] = False
    mask[-margin:, :] = False
    mask[:, :margin] = False
    mask[:, -margin:] = False

    # exclude the white scale-bar label box
    if scale_pos is not None:
        sy, xa, xb = scale_pos
        ex_y0 = max(0, sy - 52)
        ex_y1 = min(ph, sy + 20)
        ex_x1 = min(pw, xb + 45)
        mask[ex_y0:ex_y1, :ex_x1] = False
        print(f"scale label box excluded: y {ex_y0}-{ex_y1}, x 0-{ex_x1}")

    if args.exclude_rect:
        ey0, ey1, ex0, ex1 = (int(v) for v in args.exclude_rect.split(","))
        mask[ey0:ey1, ex0:ex1] = False

    # ---- segment and measure chains ----
    components = label_components(mask)
    chains = []
    chain_pixels = np.zeros_like(gray, dtype=bool)
    excluded_round = 0

    for pts in components:
        area = len(pts)
        if area < args.min_area:
            continue
        length_px = geodesic_length(pts)
        width = area / length_px
        ratio = length_px / width
        if ratio < args.min_ratio:
            excluded_round += 1
            continue

        ys, xs = pts[:, 0], pts[:, 1]
        chains.append({
            "chain_id": len(chains) + 1,
            "length_nm": length_px * nm_per_px,
            "length_px": length_px,
            "area_px": area,
            "mean_width_px": width,
            "linearity": ratio,
            "y_min": int(ys.min()) + py0,
            "y_max": int(ys.max()) + py0,
            "x_min": int(xs.min()) + px0,
            "x_max": int(xs.max()) + px0,
        })
        chain_pixels[py0 + ys, px0 + xs] = True

    # second pass: catch small dim atoms missed by the main mask
    supplemented = 0
    if args.supplement_floor > 0 and args.adaptive_radius > 0:
        supp = ((panel.astype(np.int16) - local_bg) > args.supplement_offset) & (panel > args.supplement_floor)
        supp[:margin, :] = False
        supp[-margin:, :] = False
        supp[:, :margin] = False
        supp[:, -margin:] = False
        if scale_pos is not None:
            supp[ex_y0:ex_y1, :ex_x1] = False
        for pts in label_components(supp):
            area = len(pts)
            if area < 2 or area > args.supplement_max_area:
                continue
            if chain_pixels[py0 + pts[:, 0], px0 + pts[:, 1]].any():
                continue
            length_px = geodesic_length(pts)
            width = area / length_px
            ratio = length_px / width
            ys, xs = pts[:, 0], pts[:, 1]
            chains.append({
                "chain_id": len(chains) + 1,
                "length_nm": length_px * nm_per_px,
                "length_px": length_px,
                "area_px": area,
                "mean_width_px": width,
                "linearity": ratio,
                "y_min": int(ys.min()) + py0,
                "y_max": int(ys.max()) + py0,
                "x_min": int(xs.min()) + px0,
                "x_max": int(xs.max()) + px0,
            })
            chain_pixels[py0 + ys, px0 + xs] = True
            supplemented += 1

    lengths_nm = np.array([c["length_nm"] for c in chains], dtype=float)
    mode_desc = (f"adaptive (radius {args.adaptive_radius}, offset {args.adaptive_offset})"
                 if args.adaptive_radius > 0 else "plain threshold")
    print(f"threshold: {args.threshold}, mode: {mode_desc}, min area: {args.min_area} px, "
          f"min length/width ratio: {args.min_ratio:g}")
    print(f"chains detected: {len(chains)} "
          f"(round blobs excluded: {excluded_round}, supplemented small atoms: {supplemented})")
    if len(chains):
        print(f"chain length (nm): min={lengths_nm.min():.3f}, "
              f"median={np.median(lengths_nm):.3f}, "
              f"mean={lengths_nm.mean():.3f}, max={lengths_nm.max():.3f}")

    # ---- save CSV ----
    csv_path = output_dir / f"{prefix}_stats.csv"
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "chain_id", "length_nm", "length_px", "area_px",
            "mean_width_px", "linearity", "y_min", "y_max", "x_min", "x_max",
        ])
        writer.writeheader()
        for c in chains:
            writer.writerow(c)
    print(f"saved: {csv_path}")

    # ---- save bar chart ----
    chart_path = output_dir / f"{prefix}_distribution.png"
    backend = draw_chart(lengths_nm, args.bin_width, chart_path)
    print(f"saved: {chart_path} (drawn with {backend})")
    if backend == "matplotlib":
        print(f"saved: {chart_path.with_suffix('.pdf')}")

        # ---- save overlay (needs Pillow) ----
    if not PIL_AVAILABLE:
        print("[WARN] Pillow not installed; chains overlay skipped")
    else:
        rgb = np.asarray(Image.open(input_path).convert("RGB")).copy()
        overlay = np.zeros_like(rgb)
        overlay[chain_pixels] = (255, 40, 40)
        merged = Image.fromarray(np.clip(rgb * 0.55 + overlay * 0.45, 0, 255).astype(np.uint8))
        draw = ImageDraw.Draw(merged)
        if scale_pos is not None:
            sy, xa, xb = scale_pos
            draw.rectangle(
                [px0 + xa - 4, py0 + sy - 4, px0 + xb + 4, py0 + sy + 6],
                outline=(0, 128, 255),
                width=2,
            )
            draw.text(
                (px0 + xa, py0 + sy - 20),
                f"{args.scale_nm:g} nm",
                fill=(0, 128, 255),
                font=load_font(14),
            )
        overlay_path = output_dir / f"{prefix}_overlay.png"
        merged.save(overlay_path)
        print(f"saved: {overlay_path}")


if __name__ == "__main__":
    main()
