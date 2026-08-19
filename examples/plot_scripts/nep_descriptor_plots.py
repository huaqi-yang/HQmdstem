import os, glob, math
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import warnings
warnings.filterwarnings('ignore')

# ================= 1. 统一字体设置（与示例图一致） =================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 14,          # 全局基础字号统一
    'axes.linewidth': 1.5,
    'axes.labelsize': 14,    # 坐标轴标签字号
    'xtick.major.size': 6,
    'xtick.major.width': 1.5,
    'xtick.direction': 'in',
    'ytick.major.size': 6,
    'ytick.major.width': 1.5,
    'ytick.direction': 'in',
	'ytick.right': False,
	'xtick.top': False,
    'lines.linewidth': 2.0,
    'lines.markersize': 14,
    'mathtext.fontset': 'stix',
    'legend.fontsize': 14,   # 图例字号
    'legend.frameon': False, # 图例无边框
})

# ================= 2. 配色方案 =================
MAIN_COLORS = {
    'train': '#1f77b4',
    'test': '#d62728',
    'train_alt': '#2ca02c',
    'test_alt': '#ff7f0e',
    'diagonal': '#7f7f7f',
}

# Loss 图专用配色（严格对应各物理量）
LOSS_COLORS = [
    '#1f77b4',   # 0: E-train
    '#ff7f0e',   # 1: F-train
    '#2ca02c',   # 2: V-train
    '#d62728',   # 3: L_total / L1
    '#9467bd',   # 4: L2
    '#8c564b',   # 5: E-test
    '#e377c2',   # 6: F-test
    '#7f7f7f'    # 7: V-test
]

COMPONENT_COLORS = {
    'force': ['#e41a1c', '#377eb8', '#4daf4a'],
    'stress': ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33'],
    'virial': ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33'],
    'dipole': ['#e41a1c', '#377eb8', '#4daf4a'],
    'bec': ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', 
            '#a65628', '#f781bf', '#999999']
}

TRAIN_COLORS = ['#9BBBE1', '#EAB883', '#A9CA70', '#DD7C4F', '#F09BA0', '#B58C9A']
TEST_COLORS = ['#C44E52', '#55A868', '#8172B3', '#CCB974', '#64B5CD', '#DA8BC3']

def generate_colors(data):
    if three_six_component == 0 or data == 'energy':
        return MAIN_COLORS['train'], MAIN_COLORS['test']
    else:
        if data in COMPONENT_COLORS:
            return COMPONENT_COLORS[data][:3] if data in ['force', 'dipole'] else COMPONENT_COLORS[data]
        else:
            return TRAIN_COLORS, TEST_COLORS

# ================= 3. 全局参数设置 =================
three_six_component = 0
use_range = 0
charge_sign, charge_plot_method = 1, 'hist'

plot_range = {'energy': (-9, -8), 'force': (-20, 20), 'virial': (-10, 10),
       'stress': (-10, 10), 'dipole': (-10, 10), 'polarizability': (-10, 10)}

dipole_files, polar_files = glob.glob('dipole*'), glob.glob('polarizability*')
model_type = 'dipole' if dipole_files else 'polarizability' if polar_files else None

in_file = 'gnep.in' if os.path.exists('gnep.in') else 'nep.in'
lambda_1, lambda_v, charge_mode, lambda_q, lambda_z = 1.0, 0.1, 0, 0.1, 0.5
batch = 1 if os.path.exists('gnep.in') else 1000

if os.path.exists(in_file):
    with open(in_file) as file:
        for line in file:
            line = line.strip()
            if 'lambda_1' in line and not line.startswith('#'):
                lambda_1 = float(line.split()[1])
            if 'lambda_v' in line and not line.startswith('#'):
                lambda_v = float(line.split()[1])
            if 'batch' in line and not line.startswith('#'):
                batch = int(line.split()[1])
            if 'prediction' in line and not line.startswith('#'):
                batch = 1000000
            if 'charge_mode' in line and not line.startswith('#'):
                charge_mode = int(line.split()[1])
                if 'lambda_q' in line and not line.startswith('#'):
                    lambda_q = float(line.split()[1])
                if 'lambda_z' in line and not line.startswith('#'):
                    lambda_z = float(line.split()[1])

def get_indices(data, marker):
    def get_no_indices(path):
        idx = []
        with open(path) as f:
            for i, line in enumerate(f):
                parts = line.split()
                if not parts: continue
                last = parts[-1]
                if last == marker:
                    idx.append(i)
            total = i + 1 if 'i' in locals() else 0
        return idx, total
    train_nidx, test_nidx, train_idx, test_idx = [], [], [], []
    train_len, test_len = 0, 0
    if data == 'virial' and lambda_v == 0:
        return [], [], 0, [], [], 0
    tt = ['train', 'test'] if os.path.exists(f'{data}_test.out') else ['train']
    for t in tt:
        t_nidx, t_len = get_no_indices(f'{data}_{t}.out')
        t_idx = np.setdiff1d(np.arange(t_len), t_nidx, assume_unique=True)
        if len(t_nidx) > 0:
            np.savetxt(f'{t}_no_{data}_indices.txt', t_nidx, fmt='%d')
        if t == 'train':
            train_nidx, train_idx, train_len = t_nidx, t_idx, t_len
        else:
            test_nidx, test_idx, test_len = t_nidx, t_idx, t_len
    return train_nidx, train_idx, train_len, test_nidx, test_idx, test_len

if model_type is None:
    train_nv_idx, train_v_idx, train_v_len, test_nv_idx, test_v_idx, test_v_len = get_indices('virial', '-1e+06')
    if len(train_nv_idx) == train_v_len and train_v_len > 0:
        lambda_v = 0

def set_tick_params():
    tick_params(axis='x', which='both', direction='in', top=False, bottom=True, labelsize=12)
    tick_params(axis='y', which='both', direction='in', left=True, right=False, labelsize=12)

def get_counts2two(out_file):
    file_nums = int(out_file.shape[1]//2)
    new_nep, new_dft = out_file[:, :file_nums].flatten(), out_file[:, file_nums:].flatten()
    return np.column_stack((new_nep, new_dft))

def calc_r2_rmse(out_file, data_type=None):
    file_columns = int(out_file.shape[1]//2)
    numerator = np.sum((out_file[:, :file_columns] - out_file[:, file_columns:]) ** 2)
    denominator = np.sum((out_file[:, :file_columns] - np.mean(out_file[:, :file_columns])) ** 2)
    r2_data = 1.0 if denominator == 0 else 1 - (numerator / denominator)
    rmse_origin = np.sqrt(np.mean((out_file[:, :file_columns]-out_file[:, file_columns:])**2))
    if data_type == 'stress' or data_type == 'virial':
        rmse_data = rmse_origin
    else:
        rmse_data = rmse_origin * 1000 if rmse_origin < 1 else rmse_origin
    return rmse_origin, rmse_data, r2_data

def plot_value(values, colors, data):
    columns = int(values.shape[1]//2)
    if three_six_component == 0 or data == 'energy':
        plot(values[:, 1], values[:, 0], '.', color=colors)
    else:
        for i in range(columns):
            plot(values[:, i+columns], values[:, i], '.', color=colors[i])

units = {'force': r'eV/$\mathrm{\AA}$', 'stress': 'GPa', 'energy': 'eV/atom','virial': 'eV/atom', 'dipole': 'a.u./atom', 'polarizability': 'a.u./atom', 'bec' : 'e'}
munits = {'force': r'meV/$\mathrm{\AA}$', 'stress': 'MPa', 'energy': 'meV/atom','virial': 'meV/atom', 'dipole': 'ma.u./atom', 'polarizability': 'ma.u./atom', 'bec' : 'me'}

def get_unit(data, rmse_origin):
    if data == 'stress': return 'GPa'
    return munits.get(data, 'unknown unit') if rmse_origin < 1 else units.get(data, 'unknown unit')

def get_range(data, data_file):
    return np.floor(data_file.min()), np.ceil(data_file.max())

def check_loss(loss, label, idx):
    for i in sorted(idx, reverse=True):
        if loss[-1, i] == 0.0:
            loss = np.delete(loss, i, axis=1)
            label = np.delete(label, i)
    return loss, label

# ================= 4. 修改后的 Loss 绘图部分 =================
def plot_loss():
    print('plotting loss...')
    if not os.path.exists('loss.out'):
        print("Error: loss.out not found!")
        return

    if os.path.exists('gnep.in'):
        # GNEP: Columns are typically [Gen, L_total, E_train, F_train, V_train, E_test, F_test, V_test]
        loss_data = np.loadtxt('loss.out')
        loss = loss_data[:, 1:8]
        labels = [r'$L_{\text{total}}$', 'E-train', 'F-train', 'V-train', 'E-test', 'F-test', 'V-test']
        # 映射颜色: L_total(Red), E(Blue), F(Orange), V(Green)
        colors = [LOSS_COLORS[3], LOSS_COLORS[0], LOSS_COLORS[1], LOSS_COLORS[2], LOSS_COLORS[5], LOSS_COLORS[6], LOSS_COLORS[7]]
        
        if os.path.exists('test.xyz'):
            for i in range(loss.shape[1]):
                plt.loglog(loss[:, i], color=colors[i], label=labels[i], lw=1.5)
            plt.legend(ncol=2, frameon=False, fontsize=12, loc='upper right')
        else:
            for i in range(4): # 仅 Train 部分
                plt.loglog(loss[:, i], color=colors[i], label=labels[i], lw=1.5)
            plt.legend(ncol=1, frameon=False, fontsize=12, loc='lower left')
    else:
        # Standard NEP: [Gen, Time, L1, L2, E_tr, F_tr, V_tr, ..., E_te, F_te, V_te, ...]
        loss_data = np.loadtxt('loss.out')
        if charge_mode != 0:
            loss = loss_data[:, 2:14]
            label = [r'$L_1$', r'$L_2$', 'E-train', 'F-train', 'V-train', 'Q-train', 'Z-train', 'E-test', 'F-test', 'V-test', 'Q-test', 'Z-test']
        elif model_type in ['dipole', 'polarizability']:
            loss = loss_data[:, 2:6]
            label = [r'$L_1$', r'$L_2$', f'{model_type}-train', f'{model_type}-test']
        else:
            loss = loss_data[:, 2:10]
            label = [r'$L_1$', r'$L_2$', 'E-train', 'F-train', 'V-train', 'E-test', 'F-test', 'V-test']

        # 这里的颜色逻辑：L1/L2 用紫色系，其余物理量用定义的 LOSS_COLORS
        standard_colors = [LOSS_COLORS[3], LOSS_COLORS[4], LOSS_COLORS[0], LOSS_COLORS[1], LOSS_COLORS[2], LOSS_COLORS[5], LOSS_COLORS[6], LOSS_COLORS[7]]
        
        # 动态清理无用的 loss 列
        if lambda_v == 0 and 'V-train' in label:
            v_idx = [i for i, l in enumerate(label) if 'V-' in l]
            loss, label = check_loss(loss, np.array(label), [0, 1] + v_idx)
        else:
            loss, label = check_loss(loss, np.array(label), [0, 1])
        
        # 开始绘图
        has_test = os.path.exists('test.xyz')
        num_to_plot = loss.shape[1] if has_test else (3 if charge_mode==0 else 5)
        
        for i in range(num_to_plot):
            # 自动匹配颜色，如果超出范围则循环
            plt.loglog(loss[:, i], color=standard_colors[i % len(standard_colors)], label=label[i], lw=1.5)
        
        plt.legend(ncol=2 if has_test else 1, frameon=False, fontsize=12, loc='lower left')

    set_tick_params()
    plt.xlabel('Epoch' if os.path.exists('gnep.in') else 'Generation / 100', fontsize=15)
    plt.ylabel('Loss', fontsize=15)
    plt.tight_layout()

# ================= 5. 其他函数保持不动 =================
def plot_diagonal(data):
    color_train, color_test = generate_colors(data)
    label_unit = units.get(data, 'unknown unit')
    
    if os.path.exists(f'{data}_train.out'):
        raw_train = np.loadtxt(f'{data}_train.out')
        if data in ['virial', 'stress'] and 'train_v_idx' in globals() and len(train_nv_idx) < train_v_len:
            raw_train = raw_train[train_v_idx]
        data_train = raw_train if data == 'energy' else get_counts2two(raw_train)
    else: return
    data_test = None
    if os.path.exists(f'{data}_test.out'):
        raw_test = np.loadtxt(f'{data}_test.out')
        if data in ['virial', 'stress'] and 'test_v_idx' in globals() and len(test_nv_idx) < test_v_len:
            raw_test = raw_test[test_v_idx]
        data_test = raw_test if data == 'energy' else get_counts2two(raw_test)
    
    plot_value(data_train, color_train, data)
    origin_rmse_train, rmse_data_train, r2_data_train = calc_r2_rmse(data_train, data)
    unit_label = get_unit(data, origin_rmse_train)
    
    if data_test is not None:
        plot_value(data_test, color_test, data)
        origin_rmse_test, rmse_data_test, r2_data_test = calc_r2_rmse(data_test, data)
        legend([f'train RMSE: {rmse_data_train:.3f} {unit_label}', f'test RMSE: {rmse_data_test:.3f} {unit_label}'], loc='upper left', frameon=False, fontsize=12, handletextpad=0.1)
        annotate(f'train R$^2$: {r2_data_train:.5f}', xy=(0.55, 0.18), fontsize=14, xycoords='axes fraction')
        annotate(f'test R$^2$: {r2_data_test:.5f}', xy=(0.55, 0.10), fontsize=14, xycoords='axes fraction')
    else:
        legend([f'train RMSE: {rmse_data_train:.3f} {unit_label}'], loc='upper left', frameon=False, fontsize=14, handletextpad=0.1)
        annotate(f'train R$^2$: {r2_data_train:.5f}', xy=(0.6, 0.10), fontsize=14, xycoords='axes fraction')
    ax = gca()
    xlabel(f"DFT {data} ({label_unit})", fontsize=14, labelpad=5)
    ylabel(f"NEP {data} ({label_unit})", fontsize=14)
    
    range_min, range_max = get_range(data, data_train)
    if data_test is not None:
        t_min, t_max = get_range(data, data_test)
        range_min, range_max = min(range_min, t_min), max(range_max, t_max)
    
    padding = (range_max - range_min) * 0.1
    adj_min, adj_max = range_min - padding, range_max + padding
    xlim(adj_min, adj_max); ylim(adj_min, adj_max)
    plot([adj_min, adj_max], [adj_min, adj_max], 'k--', lw=1.5, zorder=0)
    ax.set_aspect('equal', adjustable='box')
    set_tick_params()
    xticks(fontsize=12); yticks(fontsize=12)
    ax.yaxis.set_label_coords(-0.12, 0.5)

def plot_descriptor():
    """绘制PCA描述符图：修正测试集为空心圆圈，并移除图例边框"""
    if not os.path.exists('descriptor.out'): 
        return
    
    print('plotting PCA descriptor (Train & Test)...')
    
    # 1. 加载数据
    desc_train = np.loadtxt('descriptor.out')
    eng_train = np.loadtxt('energy_train.out')
    min_l_tr = min(len(desc_train), len(eng_train))
    desc_train, eng_train = desc_train[:min_l_tr], eng_train[:min_l_tr]

    has_test = os.path.exists('descriptor_test.out') and os.path.exists('energy_test.out')
    
    # 2. PCA 降维
    from sklearn.decomposition import PCA
    reducer = PCA(n_components=2)
    proj_train = reducer.fit_transform(desc_train)
    
    plt.figure(figsize=(7, 6))
    ax = plt.gca()

    # 3. 绘制训练集 (实心圆点)
    sc = plt.scatter(proj_train[:, 0], proj_train[:, 1], c=eng_train[:, 1], 
                     cmap='plasma', marker='o', edgecolor='none', alpha=0.5, s=70, label='Train Set')

    # 4. 绘制测试集 (关键修改：实现空心颜色圈)
    if has_test:
        desc_test = np.loadtxt('descriptor_test.out')
        eng_test = np.loadtxt('energy_test.out')
        min_l_te = min(len(desc_test), len(eng_test))
        desc_test, eng_test = desc_test[:min_l_te], eng_test[:min_l_te]
        
        proj_test = reducer.transform(desc_test)
        
        # 核心设置：facecolors='none' 确保不填充，edgecolors 映射颜色
        plt.scatter(proj_test[:, 0], proj_test[:, 1], 
                    facecolors='none', 
                    edgecolors=plt.cm.plasma((eng_test[:, 1] - eng_train[:, 1].min()) / (eng_train[:, 1].max() - eng_train[:, 1].min())), 
                    marker='o', linewidths=1.5, s=80, label='Test Set', alpha=0.9)

    # 5. 添加色标
    cbar = plt.colorbar(sc)
    cbar.set_label('Energy (eV/atom)', fontsize=14)
    
    # 6. 标注文本 (位置微调)
    ax.text(1.0, 0.55, 'Zn', color='#1C068D', fontsize=20, fontweight='bold', ha='center')
    ax.text(-0.65, -0.45, 'Cu', color='#F79B56', fontsize=20, fontweight='bold', ha='center')
    ax.text(-0.1, 0.15, 'Cu-Zn', color='#9F2697', fontsize=20, fontweight='bold', ha='center')

    # 7. 格式化图表
    plt.xlabel('PC1', fontsize=15)
    plt.ylabel('PC2', fontsize=15)
    
    # 8. 修正图例：去掉边框 (frameon=False)
    leg = plt.legend(loc='upper left', fontsize=12, frameon=False) 
    if has_test:
        # 让图例里的图标也显示为空心黑圈
        leg.legend_handles[1].set_facecolor('none')
        leg.legend_handles[1].set_edgecolor('black')
        leg.legend_handles[1].set_linewidth(1.5)

    set_tick_params()
    plt.savefig('nep-descriptor-final.png', dpi=300, bbox_inches='tight')
    print("描述符图已更新：Test Set 已改为真正的空心颜色圈，图例边框已移除。")

def plot_base_picture():
    abc = ['(a)', '(b)', '(c)', '(d)']
    base_diag_types = ['energy', 'force', 'virial', 'stress']
    
    if model_type is not None:
        # 单图模式
        plt.figure(figsize=(6, 6))
        ax = plt.gca()
        plot_diagonal(f'{model_type}')
        # 取消 fontweight='bold'，使用 normal
        ax.text(-0.1, 1.02, abc[0], transform=ax.transAxes, 
                fontsize=16, fontweight='normal', va='bottom', ha='right')
        plt.savefig(f'nep-{model_type}-diagonal.png', dpi=300, bbox_inches='tight')
    else:
        # 四图紧凑模式
        plt.figure(figsize=(11, 10)) # 调整画布比例
        for i, dtype in enumerate(base_diag_types):
            ax = plt.subplot(2, 2, i + 1)
            plot_diagonal(dtype)
            
            # (a)(b)(c)(d) 标签：不加粗，位置紧凑
            ax.text(-0.1, 1.02, abc[i], transform=ax.transAxes, 
                    fontsize=16, fontweight='normal', va='bottom', ha='right')
        
        # 核心修改：wspace 和 hspace 设为较小值 (0.2 左右)
        # top 留出空间给标签，bottom 留出给 x 轴
        plt.subplots_adjust(left=0.1, right=0.98, bottom=0.08, top=0.95, wspace=0.22, hspace=0.25)
        plt.savefig('nep-efvs-diagonals.png', dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    if os.path.exists('energy_train.out'): globals()['energy_train'] = np.loadtxt('energy_train.out')
    if os.path.exists('force_train.out'): globals()['force_train'] = np.loadtxt('force_train.out')
    
    plot_base_picture()
    plt.close('all')
    
    if os.path.exists('loss.out'):
        plt.figure(figsize=(5.5, 5))
        plot_loss()
        plt.savefig('nep-loss.png', dpi=300)
    
    plot_descriptor()
    print("所有图片绘制完成！")

