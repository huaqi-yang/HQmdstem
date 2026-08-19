#!/usr/bin/env python3
"""
Experimental Cu-Zn phase diagram data (assessed, solid-solid boundaries).
Source: Massalski (ASM), Okamoto (JPG 2017), Gourdon et al. (Inorg. Chem. 2007)

Format: phase_boundary_name = [(T_C, x_Zn), ...]
T in Celsius, x_Zn in at.% (0-100)

Usage: python exp_phase_diagram.py  -- saves exp_CuZn_phases.png
"""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# ================ Phase boundary data ================
# α/(α+β) — α solvus
alpha_ab = [(25, 30.0), (200, 31.0), (400, 32.5), (454, 36.8),
            (500, 37.2), (558, 37.3), (600, 36.8), (700, 35.8),
            (800, 34.5), (902, 32.1)]

# β/(α+β) — β left (Zn-poor side)
beta_left = [(454, 44.8), (500, 44.5), (558, 44.2),
             (600, 43.8), (700, 43.2), (800, 42.5), (834, 41.5)]

# β/(β+γ) — β right (Zn-rich side)
beta_right = [(454, 48.5), (500, 51.0), (558, 55.0),
              (600, 56.0), (700, 57.0), (800, 57.5), (834, 57.8)]

# γ/(β+γ) — γ left
gamma_left = [(454, 48.5), (500, 55.0), (558, 57.6),
              (600, 58.5), (700, 58.8), (800, 60.0), (834, 60.5)]

# γ/(γ+ε) — γ right / ε left (peritectoid ~558°C)
gamma_right_eut = [(558, 68.0), (600, 67.5), (700, 66.5), (800, 65.0)]
epsilon_left_eut = [(558, 73.0), (600, 74.0), (700, 75.5), (800, 77.0)]

# ε/(ε+η) — ε right / η left (eutectic ~419°C)
epsilon_right = [(419, 85.0), (500, 83.5), (558, 83.0),
                 (600, 82.5), (700, 81.5)]
eta_left = [(419, 98.3), (500, 98.5), (558, 98.5),
            (600, 99.0), (700, 99.0)]

# β' ↔ β order-disorder transition (approximate)
beta_order = [(250, 44.0), (350, 44.5), (454, 44.8),
              (454, 48.5), (350, 49.0), (250, 49.5)]

# ================ Key horizontal lines (invariant reactions) ================
# Peritectoid: β + γ ↔ ε at 558°C
peritectoid_T = 558
# Eutectoid: β ↔ α + γ at 454°C  
eutectoid_T = 454
# Eutectic: ε ↔ γ + η at ~419°C? No, actually it's peritectic: L + ε ↔ η

def plot_phase_boundary(ax, data, **kwargs):
    pts = np.array(data)
    ax.plot(pts[:, 1], pts[:, 0], **kwargs)

# ================ Plot ================
fig, ax = plt.subplots(figsize=(7, 5.5), dpi=120)

# Fill phase fields (approximate)
# α phase
x_a = np.linspace(0, 32, 100)
t_a = np.interp(x_a, [p[1] for p in alpha_ab], [p[0] for p in alpha_ab])
ax.fill_between(x_a, 0, t_a, alpha=0.1, color="#4A90D9")

# Plot boundaries
plot_phase_boundary(ax, alpha_ab, color='#333333', lw=2, label='α/(α+β)')
plot_phase_boundary(ax, beta_left, color='#333333', lw=2, label='β/(α+β)')
plot_phase_boundary(ax, beta_right, color='#333333', lw=2, label='β/(β+γ)')
plot_phase_boundary(ax, gamma_left, color='#333333', lw=2, label='γ/(β+γ)')
plot_phase_boundary(ax, gamma_right_eut, color='#333333', lw=2, ls='--', label='γ/(γ+ε)')
plot_phase_boundary(ax, epsilon_left_eut, color='#333333', lw=2, ls='--', label='ε/(γ+ε)')
plot_phase_boundary(ax, epsilon_right, color='#333333', lw=2, label='ε/(ε+η)')
plot_phase_boundary(ax, eta_left, color='#333333', lw=2, label='η/(ε+η)')

# β order-disorder
b_ord = np.array(beta_order)
ax.plot(b_ord[:, 1], b_ord[:, 0], color='#666666', lw=1.5, ls=':', label="β'↔β")

# Horizontal invariant lines
ax.axhline(eutectoid_T, color='red', lw=1, ls='--', alpha=0.4)
ax.text(10, eutectoid_T+8, f'454°C', fontsize=8, color='red', alpha=0.6)

ax.axhline(peritectoid_T, color='red', lw=1, ls='--', alpha=0.4)
ax.text(10, peritectoid_T+8, f'558°C', fontsize=8, color='red', alpha=0.6)

# Phase labels
label_pos = {
    'α': (15, 150), 'β': (46, 600), "β'": (46, 350),
    'γ': (62, 650), 'ε': (78, 600), 'η': (95, 350)
}
for name, (x, y) in label_pos.items():
    ax.text(x, y, name, fontsize=14, fontweight='bold',
            ha='center', va='center',
            bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='#cccccc', alpha=0.8))

ax.set_xlim(0, 100)
ax.set_ylim(0, 1000)
ax.set_xlabel("Zn (at.%)", fontsize=12)
ax.set_ylabel("Temperature (°C)", fontsize=12)
ax.set_title("Cu–Zn Phase Diagram (Experimental, Assessed)", fontsize=13)
ax.grid(True, alpha=0.15)

plt.tight_layout()
plt.savefig("figures/exp_CuZn_phase_diagram.png", dpi=300)
plt.savefig("figures/exp_CuZn_phase_diagram.svg")
print("Saved: figures/exp_CuZn_phase_diagram.png + .svg")

# Also save the raw data as CSV
import csv
with open("exp_CuZn_boundaries.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["boundary", "T_C", "x_Zn_at_pct"])
    for name, data in [("alpha_ab", alpha_ab), ("beta_left", beta_left),
        ("beta_right", beta_right), ("gamma_left", gamma_left),
        ("gamma_right", gamma_right_eut), ("epsilon_left", epsilon_left_eut),
        ("epsilon_right", epsilon_right), ("eta_left", eta_left)]:
        for t, x in data:
            w.writerow([name, t, x])
print("Data saved: exp_CuZn_boundaries.csv")
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import csv, os

# ================ Digitized from Gourdon 2007 Figure 1 ================
# x: at.% Zn (30-80), y: Temperature (K, 873-1273)

# α/(α+β) — α-FCC solvus
alpha_ab = [(873, 32.0), (923, 34.0), (973, 35.5), (1023, 36.5),
            (1073, 36.0), (1123, 35.0), (1173, 33.5), (1223, 32.0)]

# β/(α+β) — β left (Zn-poor side of β)
beta_left = [(923, 40.0), (973, 41.0), (1023, 41.5),
             (1073, 42.0), (1123, 42.5), (1173, 43.0), (1223, 43.0)]

# β/(β+γ) — β right (Zn-rich side of β)
beta_right = [(973, 47.0), (1023, 48.5), (1073, 49.5),
              (1123, 52.0), (1173, 54.0), (1223, 55.0)]

# γ/(β+γ) — γ left (Zn-poor side of γ)
gamma_left = [(973, 52.0), (1023, 55.0), (1073, 56.0),
              (1123, 57.0), (1173, 58.0), (1223, 59.0)]

# γ/(γ+ε) — γ right (Zn-rich side of γ)
gamma_right = [(923, 69.0), (973, 68.5), (1023, 68.0),
               (1073, 67.0), (1123, 66.5), (1173, 66.0)]

# ε/(γ+ε) — ε left (Zn-poor side of ε)
epsilon_left = [(923, 72.0), (973, 72.5), (1023, 73.0),
                (1073, 74.0), (1123, 75.5), (1173, 77.0)]

# ε/(ε+L) — ε solidus (Zn-rich side of ε)
epsilon_solidus = [(1023, 77.0), (1073, 79.0), (1123, 80.0),
                   (1173, 81.5), (1223, 82.0)]

# Liquidus (approximate, from where ε+L ends)
liquidus = [(1023, 79.0), (1073, 81.0), (1123, 83.0),
            (1173, 84.0), (1223, 85.0)]

def plot_b(ax, data, **kw):
    pts = np.array(data)
    ax.plot(pts[:, 1], pts[:, 0], **kw)

fig, ax = plt.subplots(figsize=(7, 5.5), dpi=120)

# Fill phase fields
x_a = np.linspace(30, 38, 50)
t_a = np.interp(x_a, [p[1] for p in alpha_ab], [p[0] for p in alpha_ab])
ax.fill_between(x_a, 773, t_a, alpha=0.12, color="#4A90D9", label="_α")

# Plot boundaries
plot_b(ax, alpha_ab, color='#222222', lw=2, label='α/(α+β)')
plot_b(ax, beta_left, color='#222222', lw=2, label='β/(α+β)')
plot_b(ax, beta_right, color='#222222', lw=2, label='β/(β+γ)')
plot_b(ax, gamma_left, color='#222222', lw=2, label='γ/(β+γ)')
plot_b(ax, gamma_right, color='#222222', lw=2, ls='--', label='γ/(γ+ε)')
plot_b(ax, epsilon_left, color='#222222', lw=2, ls='--', label='ε/(γ+ε)')
plot_b(ax, epsilon_solidus, color='#222222', lw=2, label='ε/(ε+L)')
plot_b(ax, liquidus, color='#222222', lw=1.5, ls=':', label='Liquidus')

# Phase labels
label_pos = {'α': (33, 930), 'β': (45, 1120), 'γ': (62, 1100),
             'ε': (78, 1050), 'L': (83, 1150)}
for name, (x, y) in label_pos.items():
    ax.text(x, y, name, fontsize=14, fontweight='bold',
            ha='center', va='center',
            bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='#ccc', alpha=0.85))

ax.set_xlim(30, 80)
ax.set_ylim(873, 1273)
ax.set_xlabel("Zn (at.%)", fontsize=12)
ax.set_ylabel("Temperature (K)", fontsize=12)
ax.set_title("Cu–Zn Phase Diagram (Gourdon 2007, Fig.1)", fontsize=12)
ax.grid(True, alpha=0.15)

plt.tight_layout()
os.makedirs("figures", exist_ok=True)
plt.savefig("figures/exp_CuZn_phase_diagram.png", dpi=300)
plt.savefig("figures/exp_CuZn_phase_diagram.svg")
print("Saved: figures/exp_CuZn_phase_diagram.png + .svg")

# Save data
with open("exp_CuZn_boundaries.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["boundary", "T_K", "x_Zn_at_pct"])
    for name, data in [("alpha_ab", alpha_ab), ("beta_left", beta_left),
        ("beta_right", beta_right), ("gamma_left", gamma_left),
        ("gamma_right", gamma_right), ("epsilon_left", epsilon_left),
        ("epsilon_solidus", epsilon_solidus), ("liquidus", liquidus)]:
        for t, x in data:
            w.writerow([name, t, x])
print("Data saved: exp_CuZn_boundaries.csv")
