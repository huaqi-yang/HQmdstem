"""
Born stability − 3×2 panels, green→red for SCRAPS-4/8/16/32.
Style: matches RDF_4x1_merged (Arial, white, no grid, fig.text labels, SVG+PDF).
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
import subprocess, re, os

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 14, 'axes.linewidth': 1.2,
    'axes.labelsize': 14, 'xtick.labelsize': 12, 'ytick.labelsize': 12,
    'xtick.major.size': 5, 'xtick.major.width': 1.2, 'xtick.direction': 'in',
    'ytick.major.size': 5, 'ytick.major.width': 1.2, 'ytick.direction': 'in',
    'lines.linewidth': 2.0, 'lines.markersize': 8,
    'mathtext.fontset': 'stix',
    'legend.fontsize': 12, 'legend.frameon': False,
})

FS_LABEL = 14; FS_TICK = 12; FS_SUB = 14; FS_LEGEND = 12; FS_INNER = 13

TEMPS = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650]
SCALES = ['4', '8', '16', '32']

# Green→red for SCRAPS sizes
COLORS = ['#4CAF50', '#8BC34A', '#FF9800', '#E53935']
MARKERS = ['o', 's', '^', 'D']

BASE_DIR = os.environ.get("CT2_BASE", "")

def run_read_thermo3(scale, temp):
    work_dir = os.path.join(BASE_DIR, str(scale), f"{temp}K", "elast")
    script_path = os.path.join(work_dir, "read_thermo3.py")
    if not os.path.isdir(work_dir): raise FileNotFoundError(f"dir: {work_dir}")
    if not os.path.isfile(script_path): raise FileNotFoundError(f"script: {script_path}")
    result = subprocess.run(["python3", "-u", "read_thermo3.py"], cwd=work_dir,
                            capture_output=True, text=True, check=True)
    output = result.stdout
    keys = ['C11','C22','C33','C44','C55','C66','C12','C13','C23']
    data = {}
    for key in keys:
        m = re.search(rf"{key}:\s*([\-]?\d+\.?\d*)", output)
        data[key] = float(m.group(1)) if m else np.nan
    return data

def calc_stability(c):
    c11,c22,c33 = float(c['C11']),float(c['C22']),float(c['C33'])
    c12,c13,c23 = float(c['C12']),float(c['C13']),float(c['C23'])
    det2 = c11*c22 - c12**2
    det3 = c11*c22*c33 + 2*c12*c13*c23 - c11*c23**2 - c22*c13**2 - c33*c12**2
    return det2, det3

print("reading data...")
all_data = {}
for scale in SCALES:
    all_data[scale] = {'C11':[],'C44':[],'C55':[],'C66':[],'DET2':[],'DET3':[]}
    for temp in TEMPS:
        try:
            c = run_read_thermo3(scale, temp)
            d2,d3 = calc_stability(c)
            all_data[scale]['C11'].append(c['C11'])
            all_data[scale]['C44'].append(c['C44'])
            all_data[scale]['C55'].append(c['C55'])
            all_data[scale]['C66'].append(c['C66'])
            all_data[scale]['DET2'].append(d2)
            all_data[scale]['DET3'].append(d3)
        except Exception as e:
            print(f"  skip {scale}/{temp}K: {e}")
            for k in ['C11','C44','C55','C66','DET2','DET3']:
                all_data[scale][k].append(np.nan)

# ── 3×2 figure, shared X, hspace=0, bottom-row only ──
FW, FH = 20/2.54, 22/2.54
fig, axes = plt.subplots(3, 2, figsize=(FW, FH), sharex=True)
axes = axes.flatten()
fig.patch.set_facecolor('white')
fig.subplots_adjust(left=0.10, bottom=0.05, right=0.98, top=0.975,
                    hspace=0.0, wspace=0.22)

# All panels use data scaling + label multipliers
# C11/C44 ~10², C55 ~10², C66 ~10³, DET2 ~10⁵, DET3 ~10⁷
plot_items = [
    ('C11',  r'$C_{11}$', r'$C_{11}$ ($\times$10$^2$ GPa)', r'$C_{11}>0$',
     -0.5, 1.0, 1e-2),
    ('C44',  r'$C_{44}$', r'$C_{44}$ ($\times$10$^2$ GPa)', r'$C_{44}>0$',
     -0.5, 1.0, 1e-2),
    ('C55',  r'$C_{55}$', r'$C_{55}$ ($\times$10$^3$ GPa)', r'$C_{55}>0$',
     -0.1, 0.5, 1e-3),
    ('C66',  r'$C_{66}$', r'$C_{66}$ ($\times$10$^3$ GPa)', r'$C_{66}>0$',
     -0.1, 0.4, 1e-3),
    ('DET2', r'$C_{11}C_{22}-C_{12}^2$',
     r'Born criteria ($\times$10$^5$ GPa$^2$)', r'$C_{11}C_{22}-C_{12}^2>0$',
     -0.1, None, 1e-5),
    ('DET3', r'$C_{11}C_{22}C_{33}+2C_{12}C_{13}C_{23}$'+'\n'+r'$-C_{11}C_{23}^2-C_{22}C_{13}^2-C_{33}C_{12}^2$',
     r'Born criteria ($\times$10$^7$ GPa$^3$)',
     r'$C_{11}C_{22}C_{33}+2C_{12}C_{13}C_{23}$'+'\n'+r'$-C_{11}C_{23}^2-C_{22}C_{13}^2-C_{33}C_{12}^2>0$',
     None, None, 1e-7)
]

abc = ['(a)','(b)','(c)','(d)','(e)','(f)']
row_idx = [-1,-1, -1,-1, -1,-1]  # rows 0 to 5
for i in range(6): row_idx[i] = i // 2

for i, item in enumerate(plot_items):
    key, ylabel, criterion = item[0], item[2], item[3]
    ybot = item[4] if len(item) > 4 else None
    ystep = item[5] if len(item) > 5 else None
    scale = item[6] if len(item) > 6 else None

    ax = axes[i]
    ax.set_facecolor('white')
    ax.grid(False)

    # y=0 dashed line
    ax.axhline(0, color='#999999', linestyle='--', linewidth=1.0, zorder=1)

    for idx, scale_name in enumerate(SCALES):
        y = np.array(all_data[scale_name][key], dtype=float)
        if scale is not None:
            y = y * scale
        vm = ~np.isnan(y)
        if np.any(vm):
            ax.plot(np.array(TEMPS)[vm], y[vm],
                    color=COLORS[idx], marker=MARKERS[idx], ms=7,
                    markerfacecolor='none', lw=2.0,
                    label=f'SCRAPS-{scale_name}', zorder=2)

    # X label only on bottom row
    if row_idx[i] == 2:
        ax.set_xlabel("Temperature (K)", fontsize=FS_LABEL)
    else:
        ax.tick_params(axis='x', labelbottom=False, which='both')
    ax.set_ylabel(ylabel, fontsize=FS_LABEL)
    ax.set_xlim(0, 700)
    ax.set_xticks(np.arange(0,701,200))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    # Y-axis: 1 decimal, no scientific offset
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    try:
        ax.yaxis.get_offset_text().set_visible(False)
    except:
        pass

    if ybot is not None:
        ax.set_ylim(bottom=ybot)
    if ystep is not None:
        ax.yaxis.set_major_locator(MultipleLocator(ystep))
    # Panel (e) DET2: use 0.4 tick spacing to avoid label overlap
    if key == 'DET2':
        ax.yaxis.set_major_locator(MultipleLocator(0.4))

    ax.tick_params(axis='both', labelsize=FS_TICK, pad=3)

    # Legend top-right
    ax.legend(fontsize=FS_LEGEND, loc='upper right', frameon=False)

    # Sub-label: (a)(b) 0.98, (e)(f) 0.92, (c)(d) 0.90
    if i < 2:
        y_label = 0.98
    elif i >= 4:
        y_label = 0.94
    else:
        y_label = 0.90
    ax.text(-0.16, y_label, abc[i], transform=ax.transAxes,
            fontsize=FS_SUB, ha='left', va='bottom')

    # Criterion text center-bottom
    ax.text(0.50, 0.35, criterion, transform=ax.transAxes,
            fontsize=FS_INNER, ha='center', va='center')

# ── Save ──────────────────────────────────────────
out_dir = os.path.join(BASE_DIR, 'RDF_figures')
os.makedirs(out_dir, exist_ok=True)
for fmt in ['svg', 'pdf', 'png']:
    out = os.path.join(out_dir, f'Born_Stability_3x2.{fmt}')
    fig.savefig(out, format=fmt, dpi=600, bbox_inches='tight')
    print(f'[OK] {out}')
plt.close()
print('Done: Born 3×2 green→red')
