"""
RDF merged — 4 rows × 1 col, hspace=0, shared X, no labels.
SCRAPS-4 → 8 → 16 → 32, green→red.
"""

import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 19.0, 'axes.linewidth': 2.6,
    'axes.labelsize': 19.0, 'xtick.labelsize': 15.0, 'ytick.labelsize': 15.0,
    'xtick.major.size': 7.0, 'xtick.major.width': 2.6, 'xtick.direction': 'in',
    'ytick.major.size': 7.0, 'ytick.major.width': 2.6, 'ytick.direction': 'in',
    'lines.linewidth': 3.6, 'legend.fontsize': 15.0, 'legend.frameon': False,
    'figure.dpi': 1800, 'savefig.dpi': 1800, 'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.02,
    'path.snap': True,
    'path.simplify': False,
    'agg.path.chunksize': 0,
})

COLORS = {'50K': '#4CAF50', '100K': '#8BC34A', '200K': '#FF9800', '250K': '#E53935'}
TEMPS = ['50K', '100K', '200K', '250K']

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rdf_data")
STRUCTURES = ['4', '8', '16', '32']
YLIMS = {'4': 9, '8': 13, '16': 13, '32': 13}

def read_rdf(filepath):
    r_vals, g_vals = [], []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            parts = line.split()
            if len(parts) < 3: continue
            try:
                r_vals.append(float(parts[1]))
                g_vals.append(float(parts[2]))
            except ValueError: continue
    return np.array(r_vals), np.array(g_vals)

# ── Load all data ───────────────────────────────
all_data = {}
for s in STRUCTURES:
    all_data[s] = {}
    for temp in TEMPS:
        path = os.path.join(BASE, s, temp, '1.txt')
        if os.path.exists(path):
            all_data[s][temp] = read_rdf(path)

# ── 4×1 figure (double-size canvas for high pixel density) ──
FW, FH = 16/2.54, 32/2.54  # 4x the area → ~7560px tall at 600DPI
fig, axes = plt.subplots(4, 1, figsize=(FW, FH),
                          sharex=True, sharey=False)
fig.patch.set_facecolor('white')
fig.subplots_adjust(left=0.16, bottom=0.06, right=0.96, top=0.995, hspace=0.0)

for i, (s, ax) in enumerate(zip(STRUCTURES, axes)):
    data = all_data[s]
    ax.set_facecolor('white')
    ax.grid(False)

    for temp in TEMPS:
        if temp in data:
            r, g = data[temp]
            mask = r < 2.30
            g_plot = g.copy()
            g_plot[mask] = 0.0
            ax.plot(r, g_plot, color=COLORS[temp], lw=2.3)

    ax.set_xlim(1.7, 3.2)
    ax.set_ylim(0, YLIMS[s])
    ax.set_ylabel(r'$g(r)$', fontsize=17.0)

    # Hide X labels on top 3, keep on bottom
    if i < 3:
        ax.tick_params(axis='x', labelbottom=False, which='both')
    else:
        ax.set_xlabel(r'$r$ ($\mathrm{\AA}$)')
    ax.tick_params(which='both', direction='in', top=True, right=True)

    # Legend (top-right of each subplot)
    for j, temp in enumerate(TEMPS):
        ypos = 0.93 - j * 0.09
        ax.plot([0.78, 0.86], [ypos, ypos], transform=ax.transAxes,
                color=COLORS[temp], lw=2.8, solid_capstyle='butt', clip_on=False)
        ax.text(0.88, ypos, temp, transform=ax.transAxes, color=COLORS[temp],
                fontsize=14.0, ha='left', va='center')

    # ── Inset ──
    if s == '4':
        # Main peak + shoulder, zoom 2.5-3.1A
        inset = ax.inset_axes([0.08, 0.32, 0.40, 0.60])
        for temp in TEMPS:
            if temp in data:
                r, g = data[temp]
                mask = (r >= 2.30) & (r <= 3.15)
                inset.plot(r[mask], g[mask], color=COLORS[temp], lw=3.2)
        inset.set_xlim(2.50, 3.10)
        inset.set_ylim(1.2, YLIMS[s] - 0.2)
        inset.tick_params(labelsize=12.5, direction='in', top=True, right=True, width=2.0)
        inset.grid(False)
        inset.text(0.98, 0.92, 'T' + chr(8593) + ': peak' + chr(8595) + '\nshoulder' + chr(8595),
                   transform=inset.transAxes, fontsize=12.5, ha='right', va='top')
        inset.annotate('', xy=(2.95, YLIMS[s]*0.3), xytext=(2.717, YLIMS[s]*0.60),
                       arrowprops=dict(arrowstyle='->', color='black', lw=3.4,
                                       connectionstyle='arc3,rad=-0.45'))
    else:
        # Single peak, zoom 2.45-2.7A — text LEFT, arrow right→left
        inset = ax.inset_axes([0.08, 0.28, 0.40, 0.65])
        peak_vals = []
        for temp in TEMPS:
            if temp in data:
                r, g = data[temp]
                mask = (r >= 2.30) & (r <= 2.75)
                inset.plot(r[mask], g[mask], color=COLORS[temp], lw=3.2)
                pk = g[(r >= 2.45) & (r <= 2.65)].max()
                peak_vals.append(pk)
        inset.set_xlim(2.45, 2.70)
        ytop = max(peak_vals) * 1.15 if peak_vals else YLIMS[s]
        inset.set_ylim(0, ytop)
        inset.tick_params(labelsize=12.5, direction='in', top=True, right=True, width=2.0)
        inset.grid(False)
        # Annotation top-LEFT of inset
        inset.text(0.02, 0.92, 'T' + chr(8593) + ': peak' + chr(8595) + '\nwidth' + chr(8593),
                   transform=inset.transAxes, fontsize=12.5, ha='left', va='top')
        # Arrow: top-right → bottom-left, arc toward left text
        if peak_vals:
            inset.annotate('', xy=(2.52, ytop*0.22), xytext=(2.62, ytop*0.593),
                           arrowprops=dict(arrowstyle='->', color='black', lw=3.4,
                                           connectionstyle='arc3,rad=0.45'))

# ── Save ──────────────────────────────────────────
out_dir = os.path.join(BASE, 'RDF_figures')
os.makedirs(out_dir, exist_ok=True)
for fmt in ['svg', 'pdf', 'png']:
    out = os.path.join(out_dir, f'RDF_4x1_v4.{fmt}')
    fig.savefig(out, format=fmt)
    print(f'[OK] {out}')
plt.close()
print('Done: RDF 4×1 merged')
