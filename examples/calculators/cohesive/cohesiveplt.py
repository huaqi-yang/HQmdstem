import numpy as np
import matplotlib.pyplot as plt
import os

# ================= 1. 统一字体设置（与示例图一致） =================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 14,          # 全局基础字号统一
    'axes.linewidth': 1.5,
    'axes.labelsize': 14,     # 坐标轴标签字号
    'xtick.major.size': 8,
    'xtick.major.width': 1.5,
    'xtick.direction': 'in',
    'ytick.major.size': 8,
    'ytick.major.width': 1.5,
    'ytick.direction': 'in',
    'lines.linewidth': 2.0,
    'lines.markersize': 8,
    'mathtext.fontset': 'stix',
    'legend.fontsize': 14,    # 图例字号
    'legend.frameon': False,  # 图例无边框
})

# ================= 2. 数据加载 =================
in_file = 'energy_vs_volume.txt'
if not os.path.exists(in_file):
    raise FileNotFoundError(f"数据文件未找到: {in_file}")

vol_energy = np.loadtxt(in_file, comments='#')
if vol_energy.ndim == 1:
    vol_energy = vol_energy.reshape(1, -1)

volume = vol_energy[:, 0]
energy = vol_energy[:, 1]

# ================= 3. 绘图 =================
plt.figure(figsize=(6, 5))
plt.plot(volume, energy, color='#1f77b4', linewidth=2.0, marker='o', markersize=4, label='FCC-Cu energy')

# 若需要最小值标注，可启用以下两行（自动找最小能量点）：
# min_idx = np.argmin(energy)
# plt.scatter(volume[min_idx], energy[min_idx], color='orange', s=80, zorder=5, label='min energy')

plt.xlabel('Volume (A^3)', fontsize=15)
plt.ylabel('Energy (eV)', fontsize=15)
plt.title('Cohesive Energy vs Volume for FCC-Cu', fontsize=16)
plt.grid(True, linestyle='--', alpha=0.5)

# 图例：若你只有一个系列，可以显示或隐藏图例
plt.legend(frameon=False, fontsize=12)

plt.tight_layout()
plt.savefig('energy_vs_volume_cohesive.png', dpi=300)
plt.close()
print("绘制完成，输出文件 energy_vs_volume_cohesive.png")
