# -*- coding: utf-8 -*-
"""
UMAP闄嶇淮 + FPS绛涢€夎缁冮泦鑴氭湰
鐢ㄦ硶:
  python umap_select.py <nep.txt> <鍏冪礌> <姹犲瓙.xyz> [宸叉湁璁粌闆?xyz ...] [--exclude 宸查€夊簭鍙?txt] [--n 300]
杈撳嚭:
  umap_descriptors.png         鈥?UMAP闄嶇淮鏁ｇ偣鍥?  *_umap.selected.xyz           鈥?FPS绛涢€夊嚭鐨勫疄闄呯粨鏋?  selected_indices.txt         鈥?瀵瑰簲鍘熷搴忓彿
"""
import sys, os, warnings, argparse, numpy as np
warnings.filterwarnings('ignore')
from ase.io import read, write
from calorine.nep import get_descriptors
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import umap

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 14,
    "axes.linewidth": 1.5,
    "axes.labelsize": 14,
    "axes.facecolor": "white",
    "axes.edgecolor": "black",
    "axes.labelcolor": "black",
    "axes.titlecolor": "black",
    "figure.facecolor": "white",
    "text.color": "black",
    "legend.facecolor": "white",
    "legend.edgecolor": "black",
    "legend.labelcolor": "black",
    "legend.fontsize": 12,
    "legend.frameon": False,
    "xtick.color": "black",
    "ytick.color": "black",
    "xtick.major.size": 6,
    "xtick.major.width": 1.5,
    "xtick.direction": "in",
    "xtick.top": False,
    "ytick.major.size": 6,
    "ytick.major.width": 1.5,
    "ytick.direction": "in",
    "ytick.right": False,
    "lines.linewidth": 2.0,
    "lines.markersize": 14,
    "mathtext.fontset": "stix",
    "grid.color": "#cccccc",
    "grid.alpha": 0.5,
})


parser = argparse.ArgumentParser(description='UMAP+FPS selection for training set')
parser.add_argument('nep', help='NEP model file nep.txt')
parser.add_argument('element', help='element symbol e.g. Li, or total')
parser.add_argument('pool', help='total pool XYZ file')
parser.add_argument('n_fps', nargs='?', type=int, default=300, help='FPS selection count (default 300)')
parser.add_argument('--mode', type=int, default=1, help='selection mode: 1=high-dim FPS (default), 2=UMAP + FPS')
parser.add_argument('train', nargs='*', help='existing train set XYZ files')
parser.add_argument('--already', help='already selected structures XYZ file to exclude', default=None)
try:
    args = parser.parse_args()
except SystemExit:
    print()
    print("="*60)
    print("Usage example (correct parameter order):")
    print("="*60)
    print()
    print("  python umap_select.py nep.txt Li pool.xyz 500 train1.xyz train2.xyz")
    print()
    print("  # For your data:")
    print("  python umap_select.py nep.txt Li trainr0sum.xyz 500 trainround2+neptrain.xyz")
    print()
    print("  python umap_select.py nep.txt Li pool.xyz 300 train.xyz --already selected.xyz")
    print()
    print("  python umap_select.py nep.txt Li pool.xyz 300 --mode 2")
    print()
    print("="*60)
    print("Note: n_fps is the 4th positional argument (after pool), NOT --n")
    print("="*60)
    print()
    sys.exit(1)

model_file = args.nep; element = args.element
cand_file = args.pool; train_files = args.train
N_FPS = args.n_fps
fps_mode = args.mode
already_xyz = args.already
if already_xyz:
    train_files = list(train_files) + [already_xyz]

print(f"NEP model: {model_file}")
print(f"Element: {element}")
print(f"Candidate pool: {cand_file}")
print(f"Existing train sets: {train_files}")

# ============ 1. 瀵规墍鏈夋枃浠舵彁鍙?NEP descriptor ============
all_descriptors = []
all_labels = []
all_indices = []   # 鍏ㄥ眬缁撴瀯绱㈠紩

def compute_descriptors(xyz_file, label):
    print(f"\nProcessing {label} ({xyz_file})...")
    atoms_list = read(xyz_file, index=':')
    n = len(atoms_list)
    descs = []
    for i, atoms in enumerate(atoms_list):
        if (i+1) % 50 == 0 or i == 0:
            print(f"  {i+1}/{n}")
        d = get_descriptors(atoms, model_filename=model_file)
        if element.lower() == 'total':
            sel = d
        else:
            sym = np.array(atoms.get_chemical_symbols())
            idx = np.where(sym == element)[0]
            if len(idx) == 0:
                continue
            sel = d[idx, :]
        descs.append(np.mean(sel, axis=0))
    descs = np.array(descs)
    all_descriptors.append(descs)
    all_labels.append(label)
    print(f"  {label}: {len(descs)} structures, shape={descs.shape}")
    return descs

# 鍊欓€夋睜
cand_desc = compute_descriptors(cand_file, os.path.basename(cand_file))
all_indices.append(np.arange(len(cand_desc)))

# Existing train sets
for tf in train_files:
    d = compute_descriptors(tf, os.path.basename(tf))
    all_indices.append(np.arange(len(d)))

# ============ 2. UMAP 闄嶇淮 ============
print("\nRunning UMAP on combined descriptors...")
combined = np.concatenate(all_descriptors, axis=0)
scaler = StandardScaler()
combined_s = scaler.fit_transform(combined)
reducer = umap.UMAP(n_components=2, random_state=42)
reduced = reducer.fit_transform(combined_s)

# 鎷嗗垎鍥炲悇鏂囦欢
reduced_splits = []
start = 0
for d in all_descriptors:
    end = start + d.shape[0]
    reduced_splits.append(reduced[start:end])
    start = end

# 淇濆瓨
np.save('umap_descriptors_2d.npy', reduced)
np.save('descriptors_raw.npy', combined)
print(f"  Reduced dim: {reduced.shape}, saved to umap_descriptors_2d.npy")
print(f"  Raw descriptors: {combined.shape}, saved to descriptors_raw.npy")

# ============ 3. FPS 绛涢€?(mode 1=楂樼淮, mode 2=UMAP) ============
mode_name = "High-dim (raw descriptor)" if fps_mode == 1 else "UMAP 2D"
print(f"\nRunning FPS - Mode {fps_mode}: {mode_name} (target {N_FPS} points)...")
cand_desc_raw = all_descriptors[0]
cand_reduced = reduced_splits[0]
n_cand = len(cand_desc_raw)

fps_space = cand_desc_raw if fps_mode == 1 else cand_reduced

exclude_set = set()
if already_xyz:
    print(f"  Matching already-selected from {already_xyz}...")
    al_atoms = read(already_xyz, index=':')
    al_descs = []
    for atoms in al_atoms:
        d = get_descriptors(atoms, model_filename=model_file)
        if element.lower() == 'total': sel = d
        else:
            sym = np.array(atoms.get_chemical_symbols())
            idx = np.where(sym == element)[0]
            if len(idx) == 0: continue
            sel = d[idx, :]
        al_descs.append(np.mean(sel, axis=0))
    al_descs = np.array(al_descs)
    for ad in al_descs:
        dists = np.sqrt(((cand_desc_raw - ad)**2).sum(axis=1))
        exclude_set.add(dists.argmin())
    print(f"  Matched {len(exclude_set)} to exclude")

keep_mask = np.ones(n_cand, dtype=bool)
if exclude_set:
    keep_mask = np.array([i not in exclude_set for i in range(n_cand)])
    fps_space = fps_space[keep_mask]
    cand_reduced = cand_reduced[keep_mask]
    cand_desc_raw = cand_desc_raw[keep_mask]
    n_cand = len(fps_space)
    print(f"  Remaining: {n_cand}")

selected = []
for k in range(N_FPS):
    if k == 0:
        dists = np.sqrt((fps_space**2).sum(axis=1))
    else:
        min_dists = np.full(n_cand, np.inf)
        for s in selected:
            d = np.sqrt(((fps_space - fps_space[s])**2).sum(axis=1))
            min_dists = np.minimum(min_dists, d)
        dists = min_dists
    best = dists.argmax()
    selected.append(best)
    if (k+1) % 50 == 0:
        print(f"  FPS {k+1}/{N_FPS}")

# --- 浠庡€欓€夋睜 XYZ 涓彁鍙栬閫変腑鐨勫抚锛屽啓鍑轰负鏂?XYZ ---
print(f"\nExtracting {len(selected)} selected structures from {cand_file}...")
cand_atoms = read(cand_file, index=':')
from ase.io import write
# If --exclude was used, map filtered indices back to original
if exclude_set:
    orig_indices = np.where(keep_mask)[0]  # filtered -> original
    sel_orig = orig_indices[selected]
else:
    sel_orig = selected
out_fps = cand_file.replace('.xyz', '_umap.selected.xyz')
sel_atoms = [cand_atoms[i] for i in sel_orig]
write(out_fps, sel_atoms)
np.savetxt('selected_indices.txt', np.array(sel_orig), fmt='%d')
print(f"  Saved {len(sel_atoms)} structures to {out_fps}")
print(f"  Selected {len(selected)} structures, indices saved to selected_indices.txt")

        
# ============ 4. 画图 ============
# --- 辅助函数：计算每个结构的颜色值（基于 descriptor 向量的模） ---
def compute_color_values(descriptors):
    vals = np.linalg.norm(descriptors, axis=1).ravel()
    vals = (vals - vals.min()) / (vals.max() - vals.min() + 1e-10)
    return vals

# 计算颜色值
combined_color_vals = compute_color_values(combined)
cand_color_vals = compute_color_values(cand_desc_raw)

cmap_name = "viridis"

# --- 图1：双面板总览 - 所有结构和 FPS 选择 ---
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), facecolor="white")
ax1, ax2 = axes

# 左侧：所有结构散点图
sc1 = ax1.scatter(reduced[:, 0], reduced[:, 1],
                  c=combined_color_vals, cmap=cmap_name, alpha=0.7, s=18,
                  edgecolors='none')
ax1.set_xlabel("UMAP 1")
ax1.set_ylabel("UMAP 2")
ax1.set_title("All structures")
ax1.set_facecolor("white")
ax1.grid(True, alpha=0.15)
cbar1 = fig.colorbar(sc1, ax=ax1, shrink=0.7, pad=0.02)
cbar1.set_label("Descriptor magnitude", color='black')
cbar1.ax.yaxis.set_tick_params(color='black')
plt.setp(cbar1.ax.get_yticklabels(), color='black')

# 右侧：FPS 选择散点图
sc2 = ax2.scatter(cand_reduced[:, 0], cand_reduced[:, 1],
                  c=cand_color_vals, cmap=cmap_name, alpha=0.5, s=15,
                  edgecolors='none', label='Pool')
# 显示已有训练集（如果有）
if train_files:
    for ri in range(1, len(all_labels)):
        rd = reduced_splits[ri]
        ax2.scatter(rd[:, 0], rd[:, 1], c='#3498DB', alpha=0.8, s=20,
                    edgecolors='black', linewidths=0.3, label=all_labels[ri])
# 显示 FPS 选中点（星形标记）
ax2.scatter(cand_reduced[selected, 0], cand_reduced[selected, 1],
            c='#2ECC71', marker='*', s=120, alpha=0.95,
            edgecolors='black', linewidths=0.5,
            label=f'FPS selected ({len(selected)})', zorder=10)
ax2.set_xlabel("UMAP 1")
ax2.set_ylabel("UMAP 2")
ax2.set_title(f"FPS selection ({len(selected)} / {n_cand})")
ax2.set_facecolor("white")
ax2.grid(True, alpha=0.15)
cbar2 = fig.colorbar(sc2, ax=ax2, shrink=0.7, pad=0.02)
cbar2.set_label("Descriptor magnitude", color='black')
cbar2.ax.yaxis.set_tick_params(color='black')
plt.setp(cbar2.ax.get_yticklabels(), color='black')

for ax in axes:
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['right'].set_color('black')
    ax.tick_params(colors='black')

ax1.legend(fontsize=8, loc='best')
ax2.legend(fontsize=8, loc='best')
plt.tight_layout()
plt.savefig('umap_descriptors.png', dpi=300, facecolor="white")
plt.close()
print("  umap_descriptors.png saved (dual-panel, white theme)")

# --- 单图1：全部结构散点图（深色主题） ---
fig, ax = plt.subplots(1, 1, figsize=(7, 5.5), facecolor="white")
sc = ax.scatter(reduced[:, 0], reduced[:, 1],
                c=combined_color_vals, cmap=cmap_name, alpha=0.7, s=20,
                edgecolors='none')
ax.set_xlabel("UMAP 1")
ax.set_ylabel("UMAP 2")
ax.set_title("UMAP projection of all structures")
ax.set_facecolor("white")
ax.grid(True, alpha=0.15)
for spine in ax.spines.values():
    spine.set_color('black')
ax.tick_params(colors='black')
cbar = fig.colorbar(sc, ax=ax, shrink=0.75, pad=0.02)
cbar.set_label("Descriptor magnitude", color='black')
cbar.ax.yaxis.set_tick_params(color='black')
plt.setp(cbar.ax.get_yticklabels(), color='black')
plt.tight_layout()
plt.savefig('umap_all.png', dpi=300, facecolor="white")
plt.close()
print("  umap_all.png saved (single, white theme)")

# --- 单图2：FPS 选择结果散点图（深色主题） ---
fig, ax = plt.subplots(1, 1, figsize=(7, 5.5), facecolor="white")
sc = ax.scatter(cand_reduced[:, 0], cand_reduced[:, 1],
                c=cand_color_vals, cmap=cmap_name, alpha=0.5, s=15,
                edgecolors='none', label='Pool')
if train_files:
    for ri in range(1, len(all_labels)):
        rd = reduced_splits[ri]
        ax.scatter(rd[:, 0], rd[:, 1], c='#3498DB', alpha=0.8, s=20,
                   edgecolors='black', linewidths=0.3, label=all_labels[ri])
ax.scatter(cand_reduced[selected, 0], cand_reduced[selected, 1],
           c='#2ECC71', marker='*', s=120, alpha=0.95,
           edgecolors='black', linewidths=0.5,
           label=f'FPS selected ({len(selected)})', zorder=10)
ax.set_xlabel("UMAP 1")
ax.set_ylabel("UMAP 2")
ax.set_title("UMAP + FPS selection")
ax.set_facecolor("white")
ax.grid(True, alpha=0.15)
for spine in ax.spines.values():
    spine.set_color('black')
ax.tick_params(colors='black')
cbar = fig.colorbar(sc, ax=ax, shrink=0.75, pad=0.02)
cbar.set_label("Descriptor magnitude", color='black')
cbar.ax.yaxis.set_tick_params(color='black')
plt.setp(cbar.ax.get_yticklabels(), color='black')
ax.legend(fontsize=8, loc='best')
plt.tight_layout()
plt.savefig('umap_fps.png', dpi=300, facecolor="white")
plt.close()
print("  umap_fps.png saved (single, white theme)")


# ============ UMAP vs NEPTRAIN(PCA) descriptor comparison ============
from sklearn.decomposition import PCA
pca = PCA(n_components=2)
pca_proj = pca.fit_transform(combined_s)

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), facecolor="white")
sc1 = axes[0].scatter(reduced[:, 0], reduced[:, 1], c=combined_color_vals,
                      cmap="viridis", alpha=0.7, s=18, edgecolors="none")
axes[0].set_xlabel("UMAP 1")
axes[0].set_ylabel("UMAP 2")
axes[0].set_title("UMAP (Viridis)")

sc2 = axes[1].scatter(pca_proj[:, 0], pca_proj[:, 1], c=combined_color_vals,
                      cmap="magma", alpha=0.7, s=18, edgecolors="none")
axes[1].set_xlabel("PC1")
axes[1].set_ylabel("PC2")
axes[1].set_title("NEPTRAIN/PCA (Magma)")

for ax in axes:
    ax.set_facecolor("white")
    ax.grid(True, alpha=0.15)
    for sp in ax.spines.values():
        sp.set_color("black")
    ax.tick_params(colors="black")

cbar1 = fig.colorbar(sc1, ax=axes[0], shrink=0.7, pad=0.02)
cbar1.set_label("Descriptor magnitude", color="black")
cbar1.ax.yaxis.set_tick_params(color="black")
cbar2 = fig.colorbar(sc2, ax=axes[1], shrink=0.7, pad=0.02)
cbar2.set_label("Descriptor magnitude", color="black")
cbar2.ax.yaxis.set_tick_params(color="black")

plt.tight_layout()
plt.savefig("umap_vs_neptrain.png", dpi=300, facecolor="white")
plt.close()
print("  umap_vs_neptrain.png saved (UMAP=viridis, NEPTRAIN/PCA=magma)")
