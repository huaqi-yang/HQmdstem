#!/usr/bin/env python3
"""
extract_cohesive_energy_multi.py
从 total_train.xyz 文件中提取每个结构的体积和总能量，
计算凝聚能：E_coh = E_tot - (n_Cu * E_Cu_bulk + n_Zn * E_Zn_bulk)
输出两列：Volume (Å^3) | Cohesive Energy (eV)
"""

import re
import numpy as np

# ================= 单质能量（从你提供的单原子文件获得） =================
# Cu 单质能量（单个原子，单位 eV）
E_Cu_bulk = -4965.3898265  # 来自 single-atom Cu total_train.xyz
# Zn 单质能量（单个原子，单位 eV）
E_Zn_bulk = -5444.9165374  # 来自 single-atom Zn total_train.xyz

def parse_structures(filename):
    """
    解析 total_train.xyz 格式的文件，返回每个结构的信息列表。
    每个结构是一个字典：{'volume': float, 'energy': float, 'n_Cu': int, 'n_Zn': int}
    """
    blocks = []
    current_block = {}
    in_block = False
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # 跳过空行
        if not line:
            i += 1
            continue
        
        # 检查是否为原子数行（数字开头）
        if line.isdigit():
            # 如果已经在处理一个块，先保存前一个块
            if in_block and 'volume' in current_block and 'energy' in current_block:
                blocks.append(current_block)
            
            # 开始新块
            current_block = {'n_Cu': 0, 'n_Zn': 0}
            in_block = True
            i += 1
            continue
        
        # 解析 Lattice 行（包含晶格和能量信息）
        if 'Lattice=' in line and in_block:
            # 提取晶格向量
            lattice_match = re.search(r'Lattice="([^"]+)"', line)
            if lattice_match:
                lattice_str = lattice_match.group(1)
                lattice_vals = list(map(float, lattice_str.split()))
                
                if len(lattice_vals) >= 9:
                    # 正交晶格：体积 = a * b * c
                    a = lattice_vals[0]
                    b = lattice_vals[4]
                    c = lattice_vals[8]
                    current_block['volume'] = a * b * c
            
            # 提取能量
            energy_match = re.search(r'energy=(-?\d+\.?\d*)', line)
            if energy_match:
                current_block['energy'] = float(energy_match.group(1))
        
        # 统计 Cu/Zn 原子
        if in_block and (line.startswith('Cu') or line.startswith('Zn')):
            parts = line.split()
            if parts:
                elem = parts[0]
                if elem == 'Cu':
                    current_block['n_Cu'] = current_block.get('n_Cu', 0) + 1
                elif elem == 'Zn':
                    current_block['n_Zn'] = current_block.get('n_Zn', 0) + 1
        
        i += 1
    
    # 添加最后一个块
    if in_block and 'volume' in current_block and 'energy' in current_block:
        blocks.append(current_block)
    
    return blocks

def main():
    # 输入文件名（根据你的实际路径修改）
    filename = 'total_train.xyz'
    
    # 解析所有结构
    blocks = parse_structures(filename)
    if not blocks:
        raise SystemExit("未解析到结构块，请检查输入格式。")
    
    # 输出文件和屏幕打印
    out_filename = 'cohesive_energy_vs_volume.txt'
    
    print(f"找到 {len(blocks)} 个结构")
    print(f"{'Volume (Å^3)':<20} {'Cohesive Energy (eV)':<25}")
    print("-" * 45)
    
    with open(out_filename, 'w') as fout:
        fout.write(f"# Volume (Å^3)    Cohesive Energy (eV)\n")
        
        for idx, b in enumerate(blocks, 1):
            V = b.get('volume')
            E_tot = b.get('energy')
            n_Cu = b.get('n_Cu', 0)
            n_Zn = b.get('n_Zn', 0)
            
            if V is None or E_tot is None:
                print(f"警告：第 {idx} 个结构缺少体积或能量信息，跳过")
                continue
            
            # 计算凝聚能
            E_coh = E_tot - (n_Cu * E_Cu_bulk + n_Zn * E_Zn_bulk)
            
            # 格式化输出
            vol_str = f"{V:.6f}"
            coh_str = f"{E_coh:.6f}"
            
            print(f"{vol_str:<20} {coh_str:<25}")
            fout.write(f"{vol_str}    {coh_str}\n")
    
    print(f"\n结果已保存到: {out_filename}")
    print(f"Cu 单质能量: {E_Cu_bulk} eV/atom")
    print(f"Zn 单质能量: {E_Zn_bulk} eV/atom")

if __name__ == "__main__":
    main()#!/usr/bin/env python3
"""
extract_cohesive_energy_multi.py
从 total_train.xyz 文件中提取每个结构的体积和总能量，
计算凝聚能：E_coh = E_tot - (n_Cu * E_Cu_bulk + n_Zn * E_Zn_bulk)
输出两列：Volume (Å^3) | Cohesive Energy (eV)
"""

import re
import numpy as np

# ================= 单质能量（从你提供的单原子文件获得） =================
# Cu 单质能量（单个原子，单位 eV）
E_Cu_bulk = -4965.3898265  # 来自 single-atom Cu total_train.xyz
# Zn 单质能量（单个原子，单位 eV）
E_Zn_bulk = -5444.9165374  # 来自 single-atom Zn total_train.xyz

def parse_structures(filename):
    """
    解析 total_train.xyz 格式的文件，返回每个结构的信息列表。
    每个结构是一个字典：{'volume': float, 'energy': float, 'n_Cu': int, 'n_Zn': int}
    """
    blocks = []
    current_block = {}
    in_block = False
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # 跳过空行
        if not line:
            i += 1
            continue
        
        # 检查是否为原子数行（数字开头）
        if line.isdigit():
            # 如果已经在处理一个块，先保存前一个块
            if in_block and 'volume' in current_block and 'energy' in current_block:
                blocks.append(current_block)
            
            # 开始新块
            current_block = {'n_Cu': 0, 'n_Zn': 0}
            in_block = True
            i += 1
            continue
        
        # 解析 Lattice 行（包含晶格和能量信息）
        if 'Lattice=' in line and in_block:
            # 提取晶格向量
            lattice_match = re.search(r'Lattice="([^"]+)"', line)
            if lattice_match:
                lattice_str = lattice_match.group(1)
                lattice_vals = list(map(float, lattice_str.split()))
                
                if len(lattice_vals) >= 9:
                    # 正交晶格：体积 = a * b * c
                    a = lattice_vals[0]
                    b = lattice_vals[4]
                    c = lattice_vals[8]
                    current_block['volume'] = a * b * c
            
            # 提取能量
            energy_match = re.search(r'energy=(-?\d+\.?\d*)', line)
            if energy_match:
                current_block['energy'] = float(energy_match.group(1))
        
        # 统计 Cu/Zn 原子
        if in_block and (line.startswith('Cu') or line.startswith('Zn')):
            parts = line.split()
            if parts:
                elem = parts[0]
                if elem == 'Cu':
                    current_block['n_Cu'] = current_block.get('n_Cu', 0) + 1
                elif elem == 'Zn':
                    current_block['n_Zn'] = current_block.get('n_Zn', 0) + 1
        
        i += 1
    
    # 添加最后一个块
    if in_block and 'volume' in current_block and 'energy' in current_block:
        blocks.append(current_block)
    
    return blocks

def main():
    # 输入文件名（根据你的实际路径修改）
    filename = 'total_train.xyz'
    
    # 解析所有结构
    blocks = parse_structures(filename)
    if not blocks:
        raise SystemExit("未解析到结构块，请检查输入格式。")
    
    # 输出文件和屏幕打印
    out_filename = 'cohesive_energy_vs_volume.txt'
    
    print(f"找到 {len(blocks)} 个结构")
    print(f"{'Volume (Å^3)':<20} {'Cohesive Energy (eV)':<25}")
    print("-" * 45)
    
    with open(out_filename, 'w') as fout:
        fout.write(f"# Volume (Å^3)    Cohesive Energy (eV)\n")
        
        for idx, b in enumerate(blocks, 1):
            V = b.get('volume')
            E_tot = b.get('energy')
            n_Cu = b.get('n_Cu', 0)
            n_Zn = b.get('n_Zn', 0)
            
            if V is None or E_tot is None:
                print(f"警告：第 {idx} 个结构缺少体积或能量信息，跳过")
                continue
            
            # 计算凝聚能
            E_coh = E_tot - (n_Cu * E_Cu_bulk + n_Zn * E_Zn_bulk)
            
            # 格式化输出
            vol_str = f"{V:.6f}"
            coh_str = f"{E_coh:.6f}"
            
            print(f"{vol_str:<20} {coh_str:<25}")
            fout.write(f"{vol_str}    {coh_str}\n")
    
    print(f"\n结果已保存到: {out_filename}")
    print(f"Cu 单质能量: {E_Cu_bulk} eV/atom")
    print(f"Zn 单质能量: {E_Zn_bulk} eV/atom")

if __name__ == "__main__":
    main()
