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
# 读取 NEP 数据
in_file_nep = 'energy_vs_volume.txt'
if not os.path.exists(in_file_nep):
    raise FileNotFoundError(f"数据文件未找到: {in_file_nep}")

vol_energy_nep = np.loadtxt(in_file_nep, comments='#')
if vol_energy_nep.ndim == 1:
    vol_energy_nep = vol_energy_nep.reshape(1, -1)

volume_nep = vol_energy_nep[:, 0]
energy_nep = vol_energy_nep[:, 1]

# 读取 DFT 数据
in_file_dft = 'energy_vs_volumeDFT.txt'
if not os.path.exists(in_file_dft):
    raise FileNotFoundError(f"数据文件未找到: {in_file_dft}")

vol_energy_dft = np.loadtxt(in_file_dft, comments='#')
if vol_energy_dft.ndim == 1:
    vol_energy_dft = vol_energy_dft.reshape(1, -1)

volume_dft = vol_energy_dft[:, 0]
energy_dft = vol_energy_dft[:, 1]

# ================= 3. 绘图 =================
plt.figure(figsize=(6, 5))

# NEP 曲线（保持蓝色不变）
plt.plot(
    volume_nep, energy_nep,
    color='#1f77b4',
    linewidth=2.0,
    marker='o',
    markersize=4,
    label='NEP'
)

# DFT 曲线（红色虚线）
plt.plot(
    volume_dft, energy_dft,
    color='red',
    linestyle='--',
    linewidth=2.0,
    label='DFT'
)

plt.xlabel('Volume (A^3)', fontsize=15)
plt.ylabel('Energy (eV)', fontsize=15)
plt.title('Cohesive Energy vs Volume for FCC-Cu', fontsize=16)
plt.grid(True, linestyle='--', alpha=0.5)

plt.legend(frameon=False, fontsize=12)

plt.tight_layout()
plt.savefig('energy_vs_volume_cohesive.png', dpi=300)
plt.close()
print("绘制完成，输出文件 energy_vs_volume_cohesive.png")
