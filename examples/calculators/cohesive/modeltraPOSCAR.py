import numpy as np
from ase.io import read, write
import os

# 1. 加载初始模型
atoms = read('model.xyz')

# 2. 定义 11 个等间距的缩放系数，范围 0.9 到 1.2
scaling_factors = np.linspace(0.9, 1.2, 11)

# 输出目录
out_dir = 'vasp_validation'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

print(f"{'Index':<10} {'Factor':<15} {'a (Å)':<15} {'b (Å)':<15} {'c (Å)':<15}")

# 3. 循环生成 POSCAR
for i, s in enumerate(scaling_factors, 1):
    # 复制初始结构
    scaled_atoms = atoms.copy()
    
    # 等向性缩放：同时缩放晶格矢量和原子坐标
    scaled_atoms.set_cell(atoms.get_cell() * s, scale_atoms=True)
    
    # 输出的 POSCAR 文件名
    filename = os.path.join(out_dir, f"POSCAR_{i}.vasp")
    
    # 写入 VASP 格式（vasp5）
    write(filename, scaled_atoms, format='vasp', vasp5=True, direct=True)
    
    # 打印进度：提取晶格向量长度
    a_len, b_len, c_len = scaled_atoms.get_cell_lengths_and_angles()[:3]
    print(f"{i:<10} {s:<15.3f} {a_len:<15.6f} {b_len:<15.6f} {c_len:<15.6f}")

print("\n所有 11 个 POSCAR 文件已生成在 'vasp_validation' 文件夹中。")
