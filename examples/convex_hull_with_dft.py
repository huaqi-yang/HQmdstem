#!/usr/bin/env python3
"""Cu-Zn convex hull: NEP89 + NEPCu-Zn + DFT (excluding 3 bad structures)."""
import numpy as np, os, math
os.makedirs("figures", exist_ok=True)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from matplotlib.lines import Line2D

from plot_style import *
apply_style(font_size=11, axes_linewidth=1.0)

C89    = "#0F4D92"
C_CUZN = "#B64342"
C_DFT  = "#2ECC71"
C_RED  = "#B64342"
C_GRAY = "#767676"
C_DARK = "#272727"
C_GREEN= "#8BCF8B"
C_TEAL = "#42949E"
C_ORANGE="#E28E2C"

def form_energy(e, n_cu, n_zn, e_cu, e_zn):
    n = n_cu + n_zn
    return e - (n_cu/n)*e_cu - (n_zn/n)*e_zn

def hull_from_ordered(ordered_pts):
    """Convex hull from ordered phases only (no SRO)."""
    arr = np.array([(p[0], p[1]) for p in ordered_pts])
    hull = ConvexHull(arr)
    hv = arr[hull.vertices]
    hv = hv[np.argsort(hv[:, 0])]
    return hv

# ===== NEP89 =====
E89_Cu = -4.2760; E89_Zn = -1.4862
nep89_ord = [
    ("alpha-Cu FCC",     4,  0,  -4.2760),
    ("alpha-Cu3Zn L12",  3,  1,  -3.6498),
    ("beta-CuZn B2",     1,  1,  -2.9709),
    ("beta'-CuZn B2",    1,  1,  -2.9762),
    ("gamma-Cu5Zn8",    20, 32,  -2.6312),
    ("delta-CuZn3",     13, 39,  -2.2529),
    ("epsilon-CuZn4",    3, 13,  -2.0629),
    ("eta-Zn HCP",       0,  2,  -1.4862),
]
nep89_sro = [
    ("SRO-8",   4, 4,  -2.9691),
    ("SRO-16",  8, 8,  -2.9701),
    ("SRO-32", 16, 16, -2.9701),
]

# ===== NEPCu-Zn =====
ECu_zn = -3.5382; EZn_zn = -1.2394
nepcu_ord = [
    ("alpha-Cu FCC",     4,  0,  -3.5382),
    ("alpha-Cu3Zn L12",  3,  1,  -3.0385),
    ("beta-CuZn B2",     1,  1,  -2.5022),
    ("beta'-CuZn B2",    1,  1,  -2.5010),
    ("gamma-Cu5Zn8",    20, 32,  -2.2307),
    ("delta-CuZn3",     13, 39,  -1.8725),
    ("epsilon-CuZn4",    3, 13,  -1.7299),
    ("eta-Zn HCP",       0,  2,  -1.2394),
]
nepcu_sro = [
    ("SRO-8",   4, 4,  -2.4047),
    ("SRO-16",  8, 8,  -2.3952),
    ("SRO-32", 16, 16, -2.3960),
]

# ===== DFT (ABACUS) =====
E_CU_DFT = -3.5506; E_ZN_DFT = -1.2690
dft_ord = [
    ("alpha-Cu FCC",     4,  0,  -3.5506),
    ("alpha-Cu3Zn L12",  3,  1,  -3.0440),
    ("beta-CuZn B2",     1,  1,  -2.4964),
    ("beta'-CuZn B2",    1,  1,  -2.4949),
    ("epsilon-CuZn4",    3, 13,  -1.7501),
    ("eta-Zn HCP",       0,  2,  -1.2690),
]
dft_sro = [
    ("SRO-4",   1, 1,  -2.4169),
    ("SRO-8",   4, 4,  -2.4016),
    ("SRO-16",  8, 8,  -2.3927),
    ("SRO-32", 16, 16, -2.3937),
]

def make_pts(ordered, sro, e_cu, e_zn):
    """Return (ordered_pts, sro_pts) lists."""
    ord_pts = [(0.0, 0.0, 4, "Cu"), (1.0, 0.0, 2, "Zn")]
    for name, n_cu, n_zn, e in ordered:
        x = n_zn/(n_cu+n_zn)
        df = form_energy(e, n_cu, n_zn, e_cu, e_zn)
        if abs(df) > 0.5: continue
        ord_pts.append((x, df, n_cu+n_zn, name))
    sro_pts = []
    for name, n_cu, n_zn, e in sro:
        x = n_zn/(n_cu+n_zn)
        df = form_energy(e, n_cu, n_zn, e_cu, e_zn)
        sro_pts.append((x, df, n_cu+n_zn, name))
    return ord_pts, sro_pts

ord89, sro89 = make_pts(nep89_ord, nep89_sro, E89_Cu, E89_Zn)
ordcu, srocu = make_pts(nepcu_ord, nepcu_sro, ECu_zn, EZn_zn)
ord_dft, sro_dft = make_pts(dft_ord, dft_sro, E_CU_DFT, E_ZN_DFT)

h89 = hull_from_ordered(ord89)
hcu = hull_from_ordered(ordcu)
h_dft = hull_from_ordered(ord_dft)

# ===== Plot =====
fig, ax = plt.subplots(figsize=(8, 5.5), dpi=100)

ax.fill_between(h89[:,0], h89[:,1], 0.08, color=C89, alpha=0.04)
ax.fill_between(hcu[:,0], hcu[:,1], 0.08, color=C_CUZN, alpha=0.04)
ax.fill_between(h_dft[:,0], h_dft[:,1], 0.08, color=C_DFT, alpha=0.04)

ax.plot(h89[:,0], h89[:,1], color=C89, lw=2.0, ls='--', alpha=0.7, zorder=2)
ax.plot(hcu[:,0], hcu[:,1], color=C_CUZN, lw=2.0, ls='--', alpha=0.7, zorder=2)
ax.plot(h_dft[:,0], h_dft[:,1], color=C_DFT, lw=2.0, ls='--', alpha=0.7, zorder=2)

# Ordered phase points
for pts, color, marker in [(ord89, C89, 'o'), (ordcu, C_CUZN, 'D'), (ord_dft, C_DFT, '^')]:
    for x, df, n_atoms, name in pts:
        if name == "Cu" or name == "Zn": continue
        ax.scatter(x, df, s=60, facecolors='none', edgecolors=color,
                   linewidths=1.5, marker=marker, alpha=0.9, zorder=3)

# SRO points (all three models)
for pts, color, marker in [(sro89, C89, 's'), (srocu, C_CUZN, 's'), (sro_dft, C_DFT, 's')]:
    for x, df, n_atoms, name in pts:
        ax.scatter(x, df, s=50, facecolors='none', edgecolors=color,
                   linewidths=2.0, marker=marker, alpha=0.75, zorder=3)

ax.text(0.02, 0.002, "Cu\n(FCC)", fontsize=9, color=C_DARK, ha='left', va='bottom', fontweight='bold')
ax.text(0.98, 0.002, "Zn\n(HCP)", fontsize=9, color=C_DARK, ha='right', va='bottom', fontweight='bold')

legend = [
    Line2D([0],[0], marker='o', color=C89, mfc='none', mec=C89, mew=2, ms=8, lw=2, ls='--', alpha=0.7, label='NEP89'),
    Line2D([0],[0], marker='D', color=C_CUZN, mfc='none', mec=C_CUZN, mew=2, ms=7, lw=2, ls='--', alpha=0.7, label='NEPCu-Zn'),
    Line2D([0],[0], marker='^', color=C_DFT, mfc='none', mec=C_DFT, mew=2, ms=7, lw=2, ls='--', alpha=0.7, label='DFT (ABACUS)'),
    Line2D([0],[0], marker='s', color=C_ORANGE, mfc='none', mec=C_ORANGE, mew=2, ms=7, lw=0, label='SRO (all)'),
]
ax.legend(handles=legend, loc='upper center', bbox_to_anchor=(0.5, 0.88), fontsize=9, ncol=4)

ax.axhline(y=0, color=C_DARK, lw=0.6, alpha=0.25)
ax.set_xlabel("Zn fraction $x_{\\rm Zn}$")
ax.set_ylabel("Relative total Energy (eV/atom)")
ax.set_xlim(-0.03, 1.03)

ax_inset = ax.inset_axes([0.0283, 0.93, 0.9434, 0.04])
for x0, x1, lbl, c in [(0,0.375,"alpha",C_DARK),(0.375,0.56,"beta",C_RED),
                        (0.56,0.72,"gamma",C_GREEN),(0.72,0.88,"epsilon",C_GRAY),
                        (0.88,1.0,"eta",C_TEAL)]:
    ax_inset.fill_between([x0,x1],0,1,color=c,alpha=0.35)
    if x1-x0>0.05: ax_inset.text((x0+x1)/2,0.5,lbl,ha='center',va='center',fontsize=9,fontweight='bold')
ax_inset.set_xlim(0,1); ax_inset.set_ylim(0,1); ax_inset.set_yticks([])
ax_inset.set_xlabel("$x_{\\rm Zn}$",fontsize=10,labelpad=-2)

finalize_figure(fig,"figures/CuZn_convex_hull_with_dft",formats=["png","svg"],dpi=300,pad=0.5)
print("Done -> figures/CuZn_convex_hull_with_dft.png + .svg")
