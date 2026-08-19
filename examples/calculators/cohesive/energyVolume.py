import numpy as np

# 假设你的初始体积 V0 是从 model.xyz 计算得到的（单位 A^3）
V0 = 25.949  # 请替换为你的实际数值

# 读取数据
data = np.loadtxt('cohesive.out')
factors = data[:, 0]
energies = data[:, 1]

# 转换为体积
volumes = V0 * (factors**3)

# 打印或保存结果
result = np.column_stack((volumes, energies))
np.savetxt('energy_vs_volume.txt', result, header='Volume(A^3) Energy(eV)')
