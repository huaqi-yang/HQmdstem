import numpy as np
import os

def read_lammps_cfg(filename):
    """读取LAMMPS格式的.cfg文件"""
    with open(filename, 'r') as f:
        lines = f.readlines()

    atoms = []
    lattice_vectors = [[0,0,0], [0,0,0], [0,0,0]]
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        i += 1

        if not line:
            continue

        if 'H0(1,1)' in line:
            # 读取晶格矩阵
            h11 = float(line.split('=')[1].strip().replace('A', ''))
            h12 = float(lines[i].split('=')[1].strip().replace('A', '')); i += 1
            h13 = float(lines[i].split('=')[1].strip().replace('A', '')); i += 1
            h21 = float(lines[i].split('=')[1].strip().replace('A', '')); i += 1
            h22 = float(lines[i].split('=')[1].strip().replace('A', '')); i += 1
            h23 = float(lines[i].split('=')[1].strip().replace('A', '')); i += 1
            h31 = float(lines[i].split('=')[1].strip().replace('A', '')); i += 1
            h32 = float(lines[i].split('=')[1].strip().replace('A', '')); i += 1
            h33 = float(lines[i].split('=')[1].strip().replace('A', '')); i += 1

            lattice_vectors = [[h11, h12, h13],
                               [h21, h22, h23],
                               [h31, h32, h33]]

        elif line.isdigit():
            # 原子类型行
            current_atom_type = int(line)
            element_symbol = lines[i].strip() # 跳过元素符号行
            i += 1

            # 提取属于该类型的所有原子
            for j in range(current_atom_type):
                coord_line = lines[i].strip()
                i += 1
                parts = coord_line.split()
                if len(parts) >= 3:
                    x, y, z = map(float, parts[:3])
                    atoms.append([current_atom_type, x, y, z])

    return np.array(lattice_vectors), np.array(atoms)

def extract_and_write_unit_cell(filename, lattice, atoms, target_box, atom_types):
    """提取指定尺寸的单胞并写入.xyz文件"""
    extracted_atoms = []
    
    # 目标单胞的尺寸
    a, b, c = target_box
    
    for atom in atoms:
        atom_id, fx, fy, fz = atom
        
        # 将分数坐标转换为原超胞的笛卡尔坐标
        x = fx * lattice[0,0] + fy * lattice[1,0] + fz * lattice[2,0]
        y = fx * lattice[0,1] + fy * lattice[1,1] + fz * lattice[2,1]
        z = fx * lattice[0,2] + fy * lattice[1,2] + fz * lattice[2,2]

        # 考虑周期性边界引起的坐标微小越界 (将接近边界的原子移回原胞)
        # 例如原点处的原子可能因为热振动跑到 -0.01
        if x < -0.5: x += lattice[0,0]
        if y < -0.5: y += lattice[1,1]
        if z < -0.5: z += lattice[2,2]
        if x > lattice[0,0] - 0.5: x -= lattice[0,0]
        if y > lattice[1,1] - 0.5: y -= lattice[1,1]
        if z > lattice[2,2] - 0.5: z -= lattice[2,2]

        # 框选：只保留在目标单胞尺寸范围内的原子
        # 加上极小的容差(0.1A)防止处于边界上的原子被漏掉
        if (-0.1 <= x < a - 0.1) and (-0.1 <= y < b - 0.1) and (-0.1 <= z < c - 0.1):
            # 将原子映射回严格的 [0, a), [0, b), [0, c) 区间内
            x = x % a
            y = y % b
            z = z % c
            extracted_atoms.append([atom_id, x, y, z])

    # 写入XYZ
    with open(filename, 'w') as f:
        f.write(f"{len(extracted_atoms)}\n")
        # 注意这里输出的 Lattice 必须是截取后的单胞晶格！
        f.write(f"Lattice=\"{a:.6f} 0.000000 0.000000 "
                f"0.000000 {b:.6f} 0.000000 "
                f"0.000000 0.000000 {c:.6f}\" ")
        f.write(f"Properties=species:S:1:pos:R:3 pbc=\"T T T\"\n")

        for atom in extracted_atoms:
            atom_id, x, y, z = atom
            element = atom_types.get(int(atom_id), f'X{int(atom_id)}')
            f.write(f"{element} {x:.6f} {y:.6f} {z:.6f}\n")
            
    return len(extracted_atoms)

# 主程序
if __name__ == "__main__":
    filename = "POSCAR46zhengjiaolammps.cfg"
    if not os.path.exists(filename):
        print(f"错误: 文件 {filename} 不存在！")
        exit(1)

    print("="*60)
    print("LAMMPS CFG 单胞提取工具")
    print("="*60)

    # 1. 读取.cfg文件
    lattice, atoms = read_lammps_cfg(filename)
    
    # 根据您的截断需求，这里设置单胞的晶格常数 (单位: Å)
    # 请根据您的实际结构进行修改！
    # 如果是 B2 相有序 CuZn，晶格常数通常约为 2.95 A
    # 如果该 cfg 模型建模时经过了特定取向的旋转（比如 x, y, z 不是 [100], [010], [001]），这里需要改成对应的正交单胞长度。
    cell_a = 10
    cell_b = 10
    cell_c = 42.584
    
    target_box = (cell_a, cell_b, cell_c)

    # 原子类型映射 (cfg中有58, 如果有Zn请在此处添加，比如 30:'Zn')
    atom_mapping = {58: 'Cu', 30: 'Zn'} 

    # 2. 提取并写入XYZ格式
    xyz_file = "primitive_cell.xyz"
    extracted_count = extract_and_write_unit_cell(xyz_file, lattice, atoms, target_box, atom_mapping)
    
    print(f"\n成功截取单胞！")
    print(f"设定单胞尺寸: a={cell_a} Å, b={cell_b} Å, c={cell_c} Å")
    print(f"截取到的原子数: {extracted_count} 个")
    print(f"输出文件: {xyz_file}")
    print("="*60)
