# HQmdstemkit 环境依赖与迁移指南（换电脑必读）

本文说明：换一台新电脑后，如何快速把 HQmdstemkit 装好并跑通。
核心思想：**一键脚本自动建环境，机器相关的路径全部放进 `config.json`**，
换电脑只需要改 `config.json` 一个文件（或重新跑一次 `install.sh`）。

## 1. 一键安装（推荐）

```bash
cd /path/to/HQmdstem   # 或你拷贝到的任何目录
bash install.sh
```

`install.sh` 会自动完成：

1. 定位 conda（PATH 或常见安装路径）；
2. 选择 Python：优先复用已有 `nepkit` 环境（不克隆 55505 个文件、不下载），其次 `HQmdstemkit`，最后系统 `python3`；
3. 只缺少依赖时才自动 `pip install numpy matplotlib scipy pillow pandas ase`，已有则跳过；
4. 检测外部程序：`gpumd`、`qstem`、`gnep`、`nep`；
5. 写入 `config.json`（环境 python 路径 + 各程序路径 + 全部内置示例/势文件路径）；
6. 把 HQmdstemkit 目录写入 `~/.bashrc`，并把 `HQmdstemkit.sh` 软链到 conda base `bin/`，当前终端立即可用；
7. 运行 `hq_env.py check` 自检。

安装完成后直接使用：

```bash
HQmdstemkit.sh
```

新终端也无需手动 `source ~/.bashrc`。

## 2. Python 环境策略（不复制大环境）

```text
优先级：nepkit 环境 > HQmdstemkit 环境 > 系统 python3
核心包：numpy, scipy, matplotlib, pillow, pandas, ase
用途：所有 HQmdstemkit 脚本的 Python 运行环境
```

install.sh 直接复用已有 `nepkit`，不执行 `conda create --clone`，因此不占额外磁盘、不用下载。

若新机器没有 nepkit，建议手动创建一个轻量环境：

```bash
conda create -n HQmdstemkit python=3.10 -y
conda activate HQmdstemkit
pip install numpy matplotlib scipy pillow pandas ase
```

之后重跑 `bash install.sh` 即可自动选用。

## 3. 外部可执行程序（不在 conda 里）

| 程序 | 用途 | 安装方式 |
|------|------|----------|
| `gpumd` | GPUMD 分子动力学 | 从 GPUMD 源码编译，放到 PATH |
| `qstem` | 高分辨电镜模拟 | 从 QSTEM 官网/源码安装，放到 PATH |
| `gnep`/`nep` | NEP 训练/推理 | nepkit 环境自带，或用 GPUMD 编译产物 |

自检方法：

```bash
./HQmdstemkit.sh 7 706
```

对应直接命令：

```bash
python3 scripts/hq_env.py check
```

## 4. config.json（机器相关路径）

所有机器相关路径集中在这一个文件里：

```json
{
  "env": "nepkit",
  "conda": "/path/to/miniconda3/bin/conda",
  "python": "/path/to/miniconda3/envs/nepkit/bin/python3",
  "gpumd": "/usr/local/bin/gpumd",
  "qstem": "qstem",
  "gnep": "gnep",
  "nep": "nep",
  "gpumd_example_dir": "examples/gpumd",
  "nep_example_dir": "examples/nep",
  "elastic_example_dir": "examples/elastic",
  "eam_input_dir": "examples/elastic/eam",
  "potential_cuzn_specific": "examples/gpumd/nep.txt",
  "potential_cuzn_generic": "examples/gpumd/nep_gen159000.txt",
  "potential_nep89_universal": "examples/nep/nep89_20250409.txt",
  "potential_nep89_restart": "examples/nep/nep89_20250409.restart",
  "potential_nep89_nepin": "examples/nep/nep89_20250409.nep.in",
  "potential_eam_cuzn": "examples/elastic/eam/CuZn.eam.alloy",
  "train_xyz_example": "examples/nep/train.xyz",
  "select_xyz_example": "examples/nep/selectsum23n16.xyz",
  "rdf_base": "examples/rdf_data",
  "chain_analysis_script": "scripts/hq_tem_core.py",
  "qstem_example_qsc": "examples/qstem/46qstem.qsc",
  "elastic_raw_csv": "examples/elastic_constants_raw_data.csv"
}
```

换电脑后只需要把这里的路径改成新电脑上的实际路径即可；
也可以用 `bash install.sh` 自动重新生成。

查看当前配置：

```bash
./HQmdstemkit.sh 7 707   # 或 python3 scripts/hq_env.py config
```

## 5. 迁移步骤清单（不跑 install.sh 时）

1. 安装 WSL + Ubuntu，安装 Miniconda；
2. 拷贝整个 `HQmdstem` 文件夹到新电脑（含 scripts/、examples/、README.md）；
3. 直接 `bash install.sh` 自动复用 nepkit；没有 nepkit 时手动创建轻量环境（见第 2 节）；
4. 安装 GPUMD、QSTEM，确认 `gpumd`/`qstem` 在 PATH 中；
5. 修改 `config.json` 中的路径（或直接 `bash install.sh`）；
6. 运行自检：
   ```bash
   ./HQmdstemkit.sh 7 706
   ```
   看到 `MISSING: none` 即环境就绪。

## 6. Windows / WSL 路径对照

```text
Windows:  C:\path\to\HQmdstem
WSL:      /path/to/HQmdstem

数据盘:   <盘符>:\...        ->  /mnt/<盘符>/...
```

HQmdstemkit 建议在 WSL/bash 里运行（和 GPUMDkit 一致）；
纯 Python 功能（结构生成、平移、筛选、弹性分析、绘图）在 Windows 下也能直接跑。