import numpy as np
import os

def read_qstem_cfg(filename):
    """读取QSTEM的.cfg文件，处理包含单位的行"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    atoms = []
    lattice_vectors = []
    reading_atoms = False
    
    for i, line in enumerate(lines):
        line = line.strip()
        
        # 跳过空行
        if not line:
            continue
            
        if 'Number of particles' in line:
            num_atoms = int(line.split('=')[1].strip())
            print(f"原子数量: {num_atoms}")
        
        elif 'A =' in line:
            # 处理晶格矢量A，移除单位
            parts = line.split('=')[1].strip()
            # 移除可能的单位
            parts = parts.replace('Angstrom', '').replace('Å', '').strip()
            parts = parts.split()
            try:
                lattice_vectors.append([float(x) for x in parts])
            except ValueError as e:
                print(f"解析晶格矢量A时出错: {line}")
                print(f"清理后: {parts}")
                raise e
        
        elif 'B =' in line:
            # 处理晶格矢量B
            parts = line.split('=')[1].strip()
            parts = parts.replace('Angstrom', '').replace('Å', '').strip()
            parts = parts.split()
            try:
                lattice_vectors.append([float(x) for x in parts])
            except ValueError as e:
                print(f"解析晶格矢量B时出错: {line}")
                raise e
        
        elif 'C =' in line:
            # 处理晶格矢量C
            parts = line.split('=')[1].strip()
            parts = parts.replace('Angstrom', '').replace('Å', '').strip()
            parts = parts.split()
            try:
                lattice_vectors.append([float(x) for x in parts])
            except ValueError as e:
                print(f"解析晶格矢量C时出错: {line}")
                raise e
        
        elif 'H0(1,1)' in line:
            # 读取晶格矩阵格式
            h11 = float(line.split('=')[1].strip())
            for j in range(1, 9):
                if i+j < len(lines):
                    if 'H0(1,2)' in lines[i+j]:
                        h12 = float(lines[i+j].split('=')[1].strip())
                    elif 'H0(1,3)' in lines[i+j]:
                        h13 = float(lines[i+j].split('=')[1].strip())
                    elif 'H0(2,1)' in lines[i+j]:
                        h21 = float(lines[i+j].split('=')[1].strip())
                    elif 'H0(2,2)' in lines[i+j]:
                        h22 = float(lines[i+j].split('=')[1].strip())
                    elif 'H0(2,3)' in lines[i+j]:
                        h23 = float(lines[i+j].split('=')[1].strip())
                    elif 'H0(3,1)' in lines[i+j]:
                        h31 = float(lines[i+j].split('=')[1].strip())
                    elif 'H0(3,2)' in lines[i+j]:
                        h32 = float(lines[i+j].split('=')[1].strip())
                    elif 'H0(3,3)' in lines[i+j]:
                        h33 = float(lines[i+j].split('=')[1].strip())
            
            lattice_vectors = [[h11, h12, h13], [h21, h22, h23], [h31, h32, h33]]
        
        elif '.atom' in line:
            reading_atoms = True
            continue
        
        elif reading_atoms and len(line.split()) >= 4:
            parts = line.split()
            try:
                atom_type = int(parts[0])
                x, y, z = map(float, parts[1:4])
                atoms.append([atom_type, x, y, z])
            except ValueError as e:
                print(f"解析原子行时出错: {line}")
                print(f"行内容: {parts}")
                raise e
    
    return np.array(lattice_vectors), np.array(atoms)

def write_xyz(filename, lattice, atoms, atom_types={1:'Cu', 2:'Zn'}):
    """写入.xyz格式文件"""
    with open(filename, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"Lattice=\"{lattice[0,0]:.6f} {lattice[0,1]:.6f} {lattice[0,2]:.6f} "
                f"{lattice[1,0]:.6f} {lattice[1,1]:.6f} {lattice[1,2]:.6f} "
                f"{lattice[2,0]:.6f} {lattice[2,1]:.6f} {lattice[2,2]:.6f}\"\n")
        for atom in atoms:
            atom_id, x, y, z = atom
            element = atom_types.get(int(atom_id), f'Atom{int(atom_id)}')
            f.write(f"{element} {x:.6f} {y:.6f} {z:.6f}\n")

def write_cif(filename, lattice, atoms, atom_types={1:'Cu', 2:'Zn'}, title="CuZn Structure"):
    """写入.cif格式文件"""
    with open(filename, 'w') as f:
        f.write(f"data_{title.replace(' ', '_')}\n")
        f.write(f"_audit_creation_method 'QSTEM export'\n")
        f.write(f"_cell_length_a {lattice[0,0]:.6f}\n")
        f.write(f"_cell_length_b {lattice[1,1]:.6f}\n")
        f.write(f"_cell_length_c {lattice[2,2]:.6f}\n")
        f.write(f"_cell_angle_alpha 90.0000\n")
        f.write(f"_cell_angle_beta 90.0000\n")
        f.write(f"_cell_angle_gamma 90.0000\n")
        f.write(f"_symmetry_space_group_name_H-M 'P1'\n")
        f.write(f"loop_\n")
        f.write(f"_atom_site_label\n")
        f.write(f"_atom_site_type_symbol\n")
        f.write(f"_atom_site_fract_x\n")
        f.write(f"_atom_site_fract_y\n")
        f.write(f"_atom_site_fract_z\n")
        
        # 计算分数坐标
        lattice_inv = np.linalg.inv(lattice)
        for i, atom in enumerate(atoms):
            atom_id, x, y, z = atom
            element = atom_types.get(int(atom_id), f'Atom{int(atom_id)}')
            
            # 从直角坐标转换到分数坐标
            cart_coords = np.array([x, y, z])
            frac_coords = np.dot(lattice_inv, cart_coords)
            
            f.write(f"{element}{i+1} {element} {frac_coords[0]:.6f} {frac_coords[1]:.6f} {frac_coords[2]:.6f}\n")

def write_poscar(filename, lattice, atoms, atom_types={1:'Cu', 2:'Zn'}):
    """写入VASP POSCAR格式文件"""
    with open(filename, 'w') as f:
        f.write(f"CuZn from QSTEM export\n")
        f.write(f"1.0\n")
        
        # 写入晶格矢量
        for i in range(3):
            f.write(f"  {lattice[i,0]:.10f} {lattice[i,1]:.10f} {lattice[i,2]:.10f}\n")
        
        # 统计原子种类和数量
        atom_counts = {}
        for atom in atoms:
            atom_id = int(atom[0])
            element = atom_types.get(atom_id, f'X{atom_id}')
            atom_counts[element] = atom_counts.get(element, 0) + 1
        
        # 写入元素种类
        elements = list(atom_counts.keys())
        f.write("  " + " ".join(elements) + "\n")
        
        # 写入各元素原子数
        counts = [str(atom_counts[el]) for el in elements]
        f.write("  " + " ".join(counts) + "\n")
        
        f.write("Direct\n")
        
        # 计算分数坐标
        lattice_inv = np.linalg.inv(lattice)
        
        # 按元素类型写入原子
        for element in elements:
            for atom in atoms:
                atom_id, x, y, z = atom
                current_element = atom_types.get(int(atom_id), f'X{atom_id}')
                if current_element == element:
                    cart_coords = np.array([x, y, z])
                    frac_coords = np.dot(lattice_inv, cart_coords)
                    f.write(f"  {frac_coords[0]:.10f} {frac_coords[1]:.10f} {frac_coords[2]:.10f}\n")

def analyze_structure(lattice, atoms, atom_types={1:'Cu', 2:'Zn'}):
    """分析结构信息"""
    print("\n" + "="*60)
    print("结构分析结果")
    print("="*60)
    
    # 晶格信息
    print(f"\n晶格矢量 (Å):")
    for i, vec in enumerate(['a', 'b', 'c']):
        print(f"  {vec}: [{lattice[i,0]:.4f}, {lattice[i,1]:.4f}, {lattice[i,2]:.4f}]")
    
    # 晶格常数
    a_len = np.linalg.norm(lattice[0])
    b_len = np.linalg.norm(lattice[1])
    c_len = np.linalg.norm(lattice[2])
    print(f"\n晶格常数 (Å):")
    print(f"  a = {a_len:.4f}")
    print(f"  b = {b_len:.4f}")
    print(f"  c = {c_len:.4f}")
    
    # 晶胞体积
    volume = np.abs(np.linalg.det(lattice))
    print(f"  晶胞体积 = {volume:.4f} Å³")
    
    # 原子信息
    print(f"\n原子信息:")
    atom_counts = {}
    for atom in atoms:
        atom_id = int(atom[0])
        element = atom_types.get(atom_id, f'Atom{atom_id}')
        atom_counts[element] = atom_counts.get(element, 0) + 1
    
    for element, count in atom_counts.items():
        print(f"  {element}: {count} 个原子")
    
    # 计算分数坐标范围
    lattice_inv = np.linalg.inv(lattice)
    frac_coords = []
    for atom in atoms:
        cart_coords = np.array(atom[1:])
        frac = np.dot(lattice_inv, cart_coords)
        frac_coords.append(frac)
    
    frac_coords = np.array(frac_coords)
    print(f"\n分数坐标范围:")
    print(f"  x: [{frac_coords[:,0].min():.4f}, {frac_coords[:,0].max():.4f}]")
    print(f"  y: [{frac_coords[:,1].min():.4f}, {frac_coords[:,1].max():.4f}]")
    print(f"  z: [{frac_coords[:,2].min():.4f}, {frac_coords[:,2].max():.4f}]")

# 主程序
if __name__ == "__main__":
    # 检查文件是否存在
    filename = "POSCAR46zhengjiaolammps.cfg"
    if not os.path.exists(filename):
        print(f"错误: 文件 {filename} 不存在！")
        print("当前目录中的文件:")
        for f in os.listdir('.'):
            if '.cfg' in f or '.CFG' in f:
                print(f"  - {f}")
        exit(1)
    
    print("="*60)
    print("QSTEM CFG 文件转换工具")
    print("="*60)
    
    try:
        # 1. 读取.cfg文件
        print(f"\n读取文件: {filename}")
        lattice, atoms = read_qstem_cfg(filename)
        
        print(f"读取成功！")
        print(f"晶格矩阵形状: {lattice.shape}")
        print(f"原子数: {len(atoms)}")
        
        # 2. 分析结构
        analyze_structure(lattice, atoms)
        
        # 3. 写入多种格式
        print(f"\n导出结构文件:")
        
        # 导出为XYZ格式
        xyz_file = "single_cell.xyz"
        write_xyz(xyz_file, lattice, atoms)
        print(f"  ✓ XYZ格式: {xyz_file}")
        
        # 导出为CIF格式
        cif_file = "single_cell.cif"
        write_cif(cif_file, lattice, atoms)
        print(f"  ✓ CIF格式: {cif_file}")
        
        # 导出为POSCAR格式
        poscar_file = "POSCAR"
        write_poscar(poscar_file, lattice, atoms)
        print(f"  ✓ POSCAR格式: {poscar_file}")
        
        # 4. 验证输出
        print(f"\n验证输出文件:")
        for out_file in [xyz_file, cif_file, poscar_file]:
            if os.path.exists(out_file):
                with open(out_file, 'r') as f:
                    first_line = f.readline().strip()
                print(f"  {out_file}: {first_line[:50]}...")
            else:
                print(f"  ✗ {out_file}: 文件创建失败")
        
        print("\n" + "="*60)
        print("转换完成！")
        print("="*60)
        
    except Exception as e:
        print(f"\n错误: {e}")
        print("\n调试信息:")
        print("尝试直接查看文件内容...")
        
        # 显示文件前20行
        try:
            with open(filename, 'r') as f:
                for i in range(20):
                    line = f.readline()
                    if not line:
                        break
                    print(f"{i+1:3}: {line.rstrip()}")
        except Exception as e2:
            print(f"无法读取文件: {e2}")
