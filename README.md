# HQmdstem

> 一个便捷的 Cu-Zn 体系一体化计算工作流工具包

**版本：1.0 beta（测试版本）**

## 简介

HQmdstem 是一个面向 Cu-Zn 体系的命令行工作流工具包，主要用于快速生成和处理 GPUMD、NEP、QSTEM、CP2K、ABACUS 等计算软件所需的输入文件，并对计算结果进行分析与可视化。对于经常进行分子动力学模拟、机器学习势训练和高分辨电镜模拟的科研人员和学生来说，它可以在结构建模、势函数训练、弹性常数计算以及微结构分析等方面提供一定便利。

HQmdstem 是一个基于命令行的菜单式工具包，风格参考 VASPKIT 与 GPUMDkit，旨在帮助用户更方便地完成 Cu-Zn 合金从结构建模、机器学习势（NEP）训练、GPUMD 分子动力学模拟，到 QSTEM 高分辨电镜模拟以及 ABACUS / CP2K 第一性原理计算的一体化流程。用户既可以通过编号菜单点选功能，也可以通过命令行方式直接调用，减少手动编辑文件的工作量。

该工具主程序为 `HQmdstemkit.sh`，启动后显示编号主菜单，包括格式转换、示例结构、工作流、计算器、分析器、可视化、工具等模块；输入对应数字即可进入二级功能。涉及 GPUMD、NEP、QSTEM 等实际计算的功能，需要预先安装对应的外部程序。

## 工具特色

- **一体化流程**：从固溶体结构建模、机器学习势训练、分子动力学模拟，到弹性常数、相图、RDF 与电镜图像分析，上下游步骤一条命令即可衔接。
- **双模式操作**：既支持编号菜单交互式点选，也支持命令行直达，便于脚本化与高通量批处理。
- **路径无关、一键迁移**：所有机器相关路径集中在 `config.json` 中，换电脑只需重跑 `install.sh` 或修改一个文件。
- **内置势函数与示例**：CuZn 专用 NEP 势、通用 NEP89 势、EAM 势及各类输入示例随包分发，开箱即用。

## 主要功能

- **结构格式转换**：支持 xyz 与 cfg 等结构文件格式互转，以及 NEP 训练集能量平移（shift）与坏点筛选（remove）。
- **固溶体结构生成**：支持有序（L12、B2、L10）与无序（Cu80Zn20）Cu-Zn 固溶体结构的自动生成。
- **GPUMD 分子动力学模拟**：自动生成 run.in 输入文件与升温脚本，用于升温观察相变过程。
- **NEP 势训练与微调**：支持神经演化势（NEP）的从零训练（train）与基于通用 NEP89 势的微调（finetune）。
- **弹性常数计算**：支持 0 K 应力-应变法与有限温度应变涨落法计算弹性常数，并进行 Born 稳定性判据分析。
- **凸包与实验相图**：支持凸包（convex hull）构建，以及实验相图数据点的自动识别与绘图。
- **QSTEM 高分辨电镜模拟**：支持 xyz → cfg → qsc 结构转换与多核并行 QSTEM 模拟。
- **CP2K 单点能计算**：支持 LDA / PBE / PBE-D3 泛函的单点能输入文件生成，以及能量、应力提取。
- **ABACUS 前处理**：支持 SCF 前处理、能量提取与能量平移。
- **微结构自动分析**：支持电镜白色链长统计、团簇尺寸分布、成分偏析、孪晶检测、晶粒尺寸与晶向统计。
- **数据可视化**：支持径向分布函数（RDF）、Born 稳定性、弹性常数、凸包与实验相图的绘图。

## 功能菜单一览

主菜单共 8 个分类：

| 编号 | 菜单 | 说明 |
|------|------|------|
| 1 | Format Conversion | 格式转换 |
| 2 | Sample Structures | 示例结构 |
| 3 | Workflow | 工作流 |
| 4 | Calculators | 计算器 |
| 5 | Analyzer | 分析器 |
| 6 | Visualization | 可视化 |
| 7 | Utilities | 工具 |
| 8 | Developing... | 开发中 |
| 0 | Exit | 退出 |

各分类下的二级功能（编号或关键字均可直达）：

| 分类 | 编号 | 功能 |
|------|------|------|
| 格式转换 1xx | 101–104 | xyz2cfg / cfg2xyz / shift 能量平移 / remove 坏点筛选 |
| 结构生成 2xx | 201–205 | L12 / B2 / L10 / Cu80Zn20 / disordered 无序固溶体 |
| 工作流 3xx | 301–313 | GPUMD prepare/run、NEP train/finetune、QSTEM、弹性 0K、ABACUS/CP2K 前处理与提取 |
| 计算 4xx | 401–407 | 弹性常数 auto/0K、凸包 hull/hull-dft、实验相图 exp、Born、strain |
| 分析 5xx | 501–508 | 链长 / 团簇 / 偏析 / 孪晶 / 晶粒 / 晶向 / 热力学 / 相变 |
| 可视化 6xx | 601–606 | RDF 4x1 / RDF / Born 图 / 弹性图 / 凸包 / 实验相图 |
| 工具 7xx | 701–707 | QSTEM prepare/run/list、说明、结构、环境自检、结果汇总 |

## 内置势函数与示例资源

| 类型 | 内容 |
|------|------|
| CuZn 专用 NEP 势 | `nep.txt`（`nep_best3n16v3.txt`，用于训练 / 微调 / 模拟） |
| 通用 NEP89 势 | `nep89_20250409`（`.txt` / `.restart` / `.nep.in`） |
| EAM 势 | `CuZn.eam.alloy`（LAMMPS 弹性常数计算） |
| GPUMD 输入示例 | `model.xyz`、`nep.txt`、`run.in` |
| NEP 训练 / 微调集 | `train.xyz`、`selectsum23n16.xyz` |

> **注意**：ABACUS 赝势 / 轨道（`apns-pseudopotentials`、`apns-orbitals`）与 CP2K 基组 / 势 / 色散参数（`BASIS_MOLOPT_UZH`、`POTENTIAL_UZH`、`dftd3.dat`）为第三方资源，本仓库**不随包分发**。请从官方渠道自行下载后放入对应目录：ABACUS 用 `--pp DIR --orb DIR` 指定路径；CP2K 下载地址见 `examples/cp2k/README.txt`。

## 适用场景

- 研究 Cu-Zn 等二元合金相变与力学性能的科研人员；
- 需要训练或微调机器学习势（NEP）并进行分子动力学模拟的学生与研究者；
- 需要进行高分辨透射电镜（QSTEM）模拟与图像分析的用户；
- 需要快速准备 ABACUS / CP2K 第一性原理计算输入文件的用户；
- 希望将结构建模、计算、分析与绘图整合为一体化流程的计算材料研究者。

## 安装与依赖

- **一键安装**：运行 `bash install.sh`，自动复用 nepkit 环境、检测外部程序、写入 `config.json` 并配置 PATH。
- **Python 依赖**：`numpy`、`scipy`、`matplotlib`、`pandas`、`ase`、`pillow`。
- **外部程序**：`gpumd`（分子动力学）、`nep` / `gnep`（NEP 训练）、`qstem`（电镜）；CP2K、ABACUS 仅需在计算集群端安装。
- **运行平台**：建议在 WSL / Linux 的 bash 环境运行（与 GPUMDkit 一致）；纯 Python 功能在 Windows 下亦可直接运行。
- **环境自检**：运行 `./HQmdstemkit.sh 7 706`，出现 `MISSING: none` 即环境就绪。

## 使用方式

首先运行安装脚本完成环境配置：

```bash
bash install.sh
```

随后可通过编号菜单点选功能，或使用命令行方式直达，例如：

```bash
./HQmdstemkit.sh                                                   # 启动主菜单，输入数字进入二级功能
./HQmdstemkit.sh 1 101 model.xyz model.cfg                         # 命令行直达：格式转换
./HQmdstemkit.sh gpumd prepare examples/gpumd --T 300 --dT 100 --Tmax 600 --steps 1100000
```

## 引用

如果您在科研工作中使用了本工具，请引用以下 NEP 相关的原始工作：

- Z. Fan, Z. Zeng, C. Zhang, Y. Wang, K. Song, H. Dong, Y. Chen, and T. Ala-Nissila, "Neuroevolution machine learning potentials: Combining high accuracy and low cost in atomistic simulations and application to heat transport," *Physical Review B* **104**, 104309 (2021). DOI: [10.1103/PhysRevB.104.104309](https://doi.org/10.1103/PhysRevB.104.104309)
- C. Chen, Y. Li, R. Zhao, Z. Liu, Z. Fan, G. Tang, and Z. Wang, "NepTrain and NepTrainKit: Automated active learning and visualization toolkit for neuroevolution potentials," *Computer Physics Communications* **317**, 109859 (2025). DOI: [10.1016/j.cpc.2025.109859](https://doi.org/10.1016/j.cpc.2025.109859)
- Z. Fan et al., "GPUMD: A package for constructing accurate machine-learned potentials and performing highly efficient atomistic simulations," *The Journal of Chemical Physics* **157**, 114801 (2022). DOI: [10.1063/5.0106617](https://doi.org/10.1063/5.0106617)
