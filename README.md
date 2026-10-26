# floquet_soc — 通用 Wannier–Floquet–SOC Engineering Framework

**状态：Stage A / B / C / D 已实现并通过 PASS/FAIL 测试；全部四个阶段都有可运行
CLI + 报告。** 框架对任意 `wannier90.win/hr/chk/eig/spn` 输入通用，不硬编码
轨道数、电子数或 SOC 类型；PbNV 只是 Benchmark 6。

## Benchmark 矩阵（真实材料 + 解析模型）

| Benchmark | 体系 | 轨道数 | SOC / 对称性 | 输入格式 | 状态 |
|---|---|---|---|---|---|
| B2 | Rashba 2 带（解析） | 2 | 面内 Rashba | 解析 | PASS |
| B3 | Rashba+Dresselhaus（解析） | 2 | PST | 解析 | PASS |
| B4 | 驱动 Rashba（解析） | 2 | 周期驱动 | 解析 | PASS |
| B6 | **PbNV** | **156** | 非中心、Ising 型 SOC | win+hr+chk+eig+spn | validate **26 PASS / 2 WARN / 0 FAIL** |
| B7 | **fcc Pt** | **18** | 中心、强 SOC 金属 | win+chk+eig+spn（无 hr.dat） | validate **24 PASS / 4 WARN / 0 FAIL** |

- **B7 Pt 数据来源**：公开仓库
  [`wannier-developers/wannier-tutorials`](https://github.com/wannier-developers/wannier-tutorials)
  的 `2022_05_Trieste/DAY2_WBERRI/advanced/spin_hall/data_Pt`；副本与授权说明见
  `benchmarks_data/Pt/README.md`（6.7 MB）。
- Pt 数据集不含 `wannier90_hr.dat`，框架用新增的 `chk+eig` 路径重建 H(R)：
  在 DFT 网格上精确（2.2e-14 eV），网格点之间是 Fourier 插值
  （与 `Pt_band.dat` 参考的 index-aligned 偏差：中位 54.1 meV、最大 1.14 eV；
  参考 run 可能使用不同 Wannier 设置，故最终 Benchmark 7 将该比较记为
  WARN——已在报告中作为已知限制说明）。

## 当前实现（四个阶段）

- **Stage A**：通用 win/hr/centres 解析；`chk` 解析、`U/v` 矩阵验证；二进制/
  格式化 `.spn` 自动识别；`S_W(k)=v†Sv` 旋转与实空间插值；静态自旋纹理；
  `validate` 输出 PASS/FAIL/WARN/BLOCKED 报告。
- **Stage B**：U(T) 演化（中点指数化 + 中点求积 + Schur）、Sambe 扩展
  Hamiltonian、高频 van Vleck 展开；A(k,E) 未跟踪谱函数（含光子权重与别名
  检查）；周期平均自旋（含 micromotion）、自旋分辨谱；Floquet 态按 overlap
  匹配到静态带（报告匹配纯度），不做“单带标签”比较。
- **Stage C**：SOC 相关能带分裂（相邻带 + 自旋配对，明确称为
  spin-related spectral splitting）；通用有效项线性拟合（用户给基函数/算符，
  报告残差与条件数）；PST 共线性/主轴/旋转分析；SOC on/off 缺同基组时显式
  返回 NEEDS_INPUT。
- **Stage D**：`scan`（ħω/光强扫描）、`converge`（n_step/m_sp 收敛，报告
  不变量与生产设置误差）、Benchmark 7（Pt，18 轨道）；跨材料统一 CLI/配置。

## 快速开始

```bash
# 输入验证（PbNV: 26 PASS/2 WARN；Pt: 24 PASS/4 WARN）
python floquet.py validate  --config configs/pbnv.json
python floquet.py validate  --config configs/pt.json

# 全部 benchmark（解析 + PbNV + Pt）
python floquet.py benchmark --which all

# 静态/Floquet 自旋纹理
python floquet.py spin_static  --config configs/pbnv.json
python floquet.py spin_floquet --config configs/pbnv.json --kmesh-n 5 --kmax 0.12

# Stage C：SOC 分裂 / PST
python floquet.py soc --config configs/pbnv.json --nk 20
python floquet.py pst --config configs/pbnv.json
python floquet.py soc --config configs/pt.json --nk 10

# Stage C：线性 Rashba 拟合（带边劈裂对 |k| 线性时给出 alpha_R，单位 eV·Å）
python floquet.py rashba --config configs/pbnv.json

# Stage D：扫描 / 收敛
python floquet.py scan     --config configs/pt.json
python floquet.py converge --config configs/pt.json

# 功能状态表；测试（约 10 分钟）
python floquet.py status
python -m unittest discover -s tests -v
```

依赖：Python ≥3.10，numpy，scipy，matplotlib。

## Linux / 服务器部署

```bash
python -m pip install .              # 或 pip install -e ".[dev]" 做开发
python -m floquet_soc --version
python -m floquet_soc validate --config configs/pt.json
python -m unittest discover -s tests -v
```

- **无显示环境**：绘图后端在所有绘图模块里强制 `Agg`，不需要 X11 / DISPLAY。
- **输出与缓存**：产物写 `$FLOQUET_SOC_OUT`（否则 `./runs`）；模型缓存在
  `$FLOQUET_SOC_CACHE`（否则系统临时目录）。
- **路径**：配置里的 `root` 与文件路径支持 `${VAR}` 和 `~` 展开，例如
  `"root": "${FLOQUET_SOC_DATA}/PbNV"`，因此同一份配置可在 Windows / Linux 共用。
- **并行**：`--jobs N` 通过 `python -m floquet_soc _kchunk` 子进程切分 k 点，
  Linux 下与 Windows 行为一致。
- **字体**：无 Arial 时自动回退 Liberation Sans / DejaVu Sans（只影响出图，不影响数值）。

## 目录结构

完整布局见 [`DIRECTORY.md`](DIRECTORY.md)。摘要：

```
floquet.py                     # 源码树 CLI 入口（安装后为 `floquet-soc` 或 `python -m floquet_soc`）
pyproject.toml / LICENSE       # 打包元数据
configs/                       # 材料配置：pbnv / pt / rashba / rashba_dresselhaus
benchmarks_data/Pt/            # Benchmark 7 公开数据副本 + README（来源/授权）
floquet_soc/                   # 包：core / drives / solvers / observables /
                               #      benchmarks / workflows / validation
tests/                         # 140 项测试（4 项数据依赖 skip）
docs/
  release/                     # 发布说明与测试结果存档（v0.1.1 … v0.1.7）
  reports/                     # Stage A–D 与 Pt / PbNV 专题报告
  STATUS_REGISTRY.md / NEW_MATERIAL_GUIDE.md / REQUIREMENTS_GAP_MATRIX.md …
outputs/                       # 运行产物（git 忽略；每个计算一个子目录）
```

## 重要修正（相对早期报告）

Stage A/B 时期的 TRS 质量判据曾错误地写成
`max|eig S(-k)+eig S(k)|`；正确条件是
`spec S(-k) = -spec S(k)`，即 `sorted(S(-k)) = -reverse(sorted(S(k)))`。
修正后：

- PbNV（新 run，VASP LWRITE_SPN）：TRS 偏差 **0.002**（旧判据误报 2.000），
  `validate` 从 24 PASS/4 WARN 变为 **26 PASS / 2 WARN**（其中
  `spin_wannier_hermiticity_offgrid` 为后加的 off-grid 检查，自带
  NOT_VALIDATED 注记）；
- Pt（QE pw2wannier90）：TRS 偏差 **0.0000**；
- 早期“PbNV `.spn` ~1% TRS 残差”的说法部分源于该判据 bug，已在
  [`docs/reports/STAGE_A_REPORT.md`](docs/reports/STAGE_A_REPORT.md) §11 与
  [`docs/reports/STAGE_D_REPORT.md`](docs/reports/STAGE_D_REPORT.md) §4 更正。

## 已知限制

1. hr-less（Pt）重建的带间插值精度有限（与 `Pt_band.dat` 的 index-aligned
   偏差中位 54.1 meV / 最大 1.14 eV，且参考 run 的设置可能不同）；定量带结构
   分析需补 `hr.dat`。
2. Pt 为金属，`static_gap` 报告为 WARN；“SOC 分裂”一律称为
   spin-related spectral splitting，不冒充 SOC 常数。
3. 自旋纹理/PST 分析对近简并态使用子空间/去简并判据；Floquet 态必须按
   overlap 匹配（已实现并报告纯度）。
4. 可选：同一 tutorial 仓库的 bcc Fe（18 spinor Wannier、磁性 SOC，
   `mp_grid=8×8×8`）可作为磁性 Benchmark 8 的候选数据。

## 最终验收状态（2026-10-23）

详见 [`docs/reports/PT_ACCEPTANCE_REPORT.md`](docs/reports/PT_ACCEPTANCE_REPORT.md)、
[`docs/reports/PT_INTERPOLATION_REPORT.md`](docs/reports/PT_INTERPOLATION_REPORT.md)、
[`docs/reports/PBNV_PST_REPORT.md`](docs/reports/PBNV_PST_REPORT.md)、
[`docs/reports/PHOTON_TRUNCATION_REPORT.md`](docs/reports/PHOTON_TRUNCATION_REPORT.md)。

- **软件实现**：PASS（WS/ndegen Hermiticity、自旋算符；全套 unittest
  **140 项 OK，4 项数据依赖 skip，0 failure**，见
  [`docs/release/RELEASE_TEST_RESULTS.txt`](docs/release/RELEASE_TEST_RESULTS.txt)）。
- **静态对称性**：PbNV PASS（T-Kramers 仅限 TRIM，残差 0.31–0.49 meV）；
  Pt 部分 WARN（DFT 网格 Kramers 精确，离网格残差 0.43 meV 中位是
  数值插值误差，`pt_kramers_offgrid=WARN`）。
- **Floquet 数值收敛**：PASS（n_step 收敛；m_sp 需 >=32 才满足
  光子权重完整，m_sp=16 有 2 个态权重缺失 ~100%）。
- **真实材料定量精度**：尚未达到。Pt 官方 Wannier90 hr.dat 未生成，
  `Pt_band.dat` 无法复现 → Pt 仅作软件接口/错误检测 benchmark；
  PbNV PST 为近似 Ising-like（WARN，轴稳定 0.998，共线性规范涨落 0.038，
  bootstrap 0.942–0.973，近 Γ 区共线性 0.865 vs 远区 0.964）。
