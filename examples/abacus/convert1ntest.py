from ase.io import read, write
from ase.calculators.singlepoint import SinglePointCalculator
import numpy as np

# 你的参考值 (确保单位为 eV)
E_Ti_ref = -1562.941684
E_Al_ref = -59.386039
E_V_ref = -2000.204008
E_C_ref= -155.0907296 
E_Zn_ref= -5444.916537
E_Cu_ref= -4965.389826

# 读取所有帧
try:
    frames = read('total_train.xyz', index=':')
except Exception as e:
    print(f"读取文件失败: {e}")
    exit()

processed_count = 0

for i, atoms in enumerate(frames):
    # 统计原子数量
    n_Ti = sum(atoms.symbols == 'Ti')
    n_Al = sum(atoms.symbols == 'Al')
    n_V = sum(atoms.symbols == 'V')
    n_C = sum(atoms.symbols == 'C')    
    n_Cu = sum(atoms.symbols == 'Cu')
    n_Zn = sum(atoms.symbols == 'Zn')
    # 1. 获取能量 (多重保险)
    old_energy = atoms.info.get('energy') or atoms.info.get('Energy')
    if old_energy is None and atoms.calc is not None:
        old_energy = atoms.calc.results.get('energy')

    if old_energy is None:
        print(f"警告: 第 {i} 帧缺失能量信息，已跳过。")
        continue

    # 2. 获取受力 (多重保险，决定坐标后的 forces 列)
    old_forces = atoms.arrays.get('forces') or atoms.arrays.get('force')
    if old_forces is None and atoms.calc is not None:
        old_forces = atoms.calc.results.get('forces')

    # 3. 关键修正：获取应力 (决定抬头中的 stress 字段)
    # 优先从 info 中提取原始的 stress 字符串或数组
    old_stress = atoms.info.get('stress') or atoms.info.get('Stress')
    if old_stress is None and atoms.calc is not None:
        old_stress = atoms.calc.results.get('stress')

    # 执行平移计算
    shifted_energy = old_energy - (n_Cu * E_Cu_ref + n_Zn * E_Zn_ref + n_C * E_C_ref + n_Ti * E_Ti_ref + n_Al * E_Al_ref + n_V * E_V_ref)

    # 清理 info 字典，防止旧标签与新计算器冲突
    # 注意：这里要清理 stress 相关的键，因为我们要把它通过 results 统一传入
    for key in ['energy', 'Energy', 'stress', 'Stress']:
        if key in atoms.info: del atoms.info[key]

    # 4. 统一打包：将所有物理量放入 results
    # 只有出现在这里的项，才会被 ASE 写回 Header 和原子行 [User Query]
    results = {
        'energy': shifted_energy, 
        'forces': old_forces
    }
    if old_stress is not None:
        results['stress'] = old_stress

    # 重新绑定计算器
    atoms.calc = SinglePointCalculator(atoms, **results)

    # 保持位力信息在 info 字典中 (如原文件有 Virial/virial)
    if 'Virial' in atoms.info:
        atoms.info['virial'] = atoms.info.pop('Virial')

    processed_count += 1

# 写入新文件
write('trainitest_manual_shifted.xyz', frames, format='extxyz')
print(f"处理完成！成功平移了 {processed_count} 帧，已生成 train_manual_shifted.xyz")
