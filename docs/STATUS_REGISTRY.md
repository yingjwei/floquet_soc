# STATUS_REGISTRY — floquet_soc v0.1.0-rc1 非 PASS 状态登记

本文件是冻结快照 `floquet_soc v0.1.0-rc1` 的**非 PASS 状态登记表**。所有
`WARN` / `NOT_APPLICABLE` / `NOT_VALIDATED` / `UNCONVERGED` / `BLOCKED` /
`NEEDS_INPUT` / `FAIL` / `SKIP` 状态在此逐条保留：来源文件、含义、限制类别
（物理 / 数值 / 数据-环境 / 接口约定）、解决条件。**任何后续研究引用本框架的
数字或结论前，必须先读本表；不得把 WARN 写成 PASS，也不得引用已撤回的旧 PASS。**

冻结基线：

- 版本：`VERSION` = `0.1.0-rc1`；`floquet_soc/version.py::__version__` 与
  `floquet_soc/__init__.py::__version__` 均导出 `0.1.0-rc1`（冻结时已统一）。
- git：annotated tag `v0.1.0-rc1`（提交信息 `floquet_soc v0.1.0-rc1 frozen
  baseline`）；精确 commit hash 见快照目录外的 `FROZEN_SNAPSHOT.md`。
- 测试：`python -m unittest discover -s tests -v` → **Ran 74 tests, OK
  (skipped=2), failures=0, errors=0**，输出冻结为 `RELEASE_TEST_RESULTS.txt`
  （UTF-8）。2 个 skip 见 §7 的 SKIP-01/02。
- 测试运行环境注记：冻结环境（DSH Windows 沙箱）禁止命名管道，因此
  `ProcessPoolExecutor` 不可用；验收运行设置 `FLOQUET_SOC_JOBS=1` 强制串行。
  该变量是框架已支持的环境覆盖（见 `floquet_soc/cli.py::_jobs_for`），
  不改变任何物理算法或数值结果。

限制类别：

| 类别 | 含义 |
|---|---|
| 物理 | 由体系对称性/物理定义决定，换代码也消不掉（如金属无带隙、PT 体系 SOC 分裂不适用） |
| 数值 | 截断/插值/收敛误差，可通过更密网格、更大窗口、官方 hr.dat 降低 |
| 数据-环境 | 输入数据缺失、外部可执行文件缺失、历史 run 不可复现 |
| 接口约定 | 框架显式拒绝/如实返回的状态（NEEDS_INPUT、BLOCKED、NOT_VALIDATED 等） |

证据层级：**当前证据** = `outputs/PbNV/`、`outputs/Pt/`、`outputs/benchmarks/`
（冻结版本的 validate/benchmark 输出）；**历史证据** = `*_verify/`、
`outputs/_tests*/validate_pre_freeze_*`、`pt_audit.json`、日期化报告
（`PT_*`、`STAGE_*`）。历史证据保留原样，不与当前结论混用。

---

## 1. 当前 validate 状态

### 1.1 PbNV — `outputs/PbNV/validate.md` / `.json`

汇总：**PASS=26, WARN=2, FAIL=0, BLOCKED=0, SKIP=0**。

| ID | 检查 | 状态 | 证据 | 类别 | 含义 / 影响 | 解决条件 |
|---|---|---|---|---|---|---|
| W-01 | `spn_spin_sum_rule` | **WARN** | `eig(ΣS²)∈[0.719,3.000]` | 物理+数据 | 输入 VASP LWRITE_SPN 的赝波函数近似 + 156 维 Wannier 窗未覆盖全部自旋权重（文件为 192 带基）。自旋纹理/PST 只能是半定量 | 提供自旋完备基/更大 Wannier 窗，或明确接受半定量；不能靠放宽容差解决 |
| W-02 | `spin_wannier_sum_rule` | **WARN** | `eig(ΣS²)∈[2.786,3.000]` | 物理+数据 | 旋转到 Wannier 基后更接近 3，但仍缺窗口外自旋权重；与 W-01 同源 | 同 W-01 |
| N-01 | `spin_wannier_hermiticity_offgrid` | **PASS（携带 NOT_VALIDATED）** | detail = `NOT_VALIDATED: shifted-mesh alias fallback (off-grid \|S-S†\| = 4.26e-02); grid values are exact but off-grid spin physics is not validated` | 数值/后端 | **PASS 仅指 DFT 网格上**；离网格的自旋算符物理未验证。禁止用离网格自旋做定量结论 | 换无位移网格 + WS/ndegen 后端，或提供官方变换文件；在此之前按未验证处理 |

### 1.2 Pt — `outputs/Pt/validate.md` / `.json`

汇总（2026-10-23 冻结复跑，当前代码）：**PASS=24, WARN=4, FAIL=0, BLOCKED=0,
SKIP=0**。旧文件记录 23 PASS/4 WARN（phase-2 之前的 validate 输出，缺
`spin_wannier_hermiticity_offgrid` 行），已原样保留为
`outputs/Pt/validate_pre_freeze_2026-10-23.md/.json`（见 H-F-02）。

| ID | 检查 | 状态 | 证据 | 类别 | 含义 / 影响 | 解决条件 |
|---|---|---|---|---|---|---|
| W-03 | `wannier_centres` | **WARN** | `tau source = wannier90.chk centres` | 数据来源 | 无 `wannier90_centres.xyz`，中心从 chk 读取（等价、仅来源不同）；不影响已通过的 H/自旋检查 | 提供 `wannier90_centres.xyz` 可消除（可选） |
| W-04 | `static_gap` | **WARN** | `gap = -2.6952 eV`（路径上无隙） | 物理 | fcc Pt 是金属，无带隙；如实标注，不是数值错误 | 不适用（换半导体/绝缘体材料才可能有 gap） |
| W-05 | `spn_spin_sum_rule` | **WARN** | `eig(ΣS²)∈[0.333,3.000]` | 物理+数据 | QE `pw2wannier90` 赝波函数/PAW 近似；wannier 窗不完整。自旋纹理半定量 | 提供自旋完备/更高质量 spn；接受半定量 |
| W-06 | `spin_wannier_sum_rule` | **WARN** | `eig(ΣS²)∈[2.918,3.000]` | 物理+数据 | 旋转后接近 3，仍非严格 3；与 W-05 同源 | 同 W-05 |

> 冻结复跑新增（PASS，非登记项的 NOT_VALIDATED）：
> `spin_wannier_hermiticity_offgrid` = PASS，`max|S(k)-S(k)†| off-grid = 6.66e-16`
> （WS/ndegen 后端，已做离网格验证）；这正是 Pt validate 从旧 23 PASS 变为
> **24 PASS** 的唯一原因。

---

## 2. 当前 benchmark 状态

### 2.1 Benchmark 7：fcc Pt — `outputs/benchmarks/pt/report.md` / `.json`

汇总：**PASS=13, WARN=5, FAIL=0, BLOCKED=0, SKIP=0**（`benchmarks_verify/pt/`
的 9 PASS / 1 WARN 是更早版本的副本，见 §2.5）。

| ID | 检查 | 状态 | 证据 | 类别 | 含义 / 影响 | 解决条件 |
|---|---|---|---|---|---|---|
| W-07 | `pt_hr_reconstruction` | **WARN** | 后端标签：grid-consistent WS/ndegen reconstruction；`official-hr comparison pending` | 数据-环境 | 公开数据无官方 `wannier90_hr.dat`，H(R) 由 chk+eig 重建（网格上精确、网格间为 Fourier 插值） | 在 Linux/HPC 用 `restart=plot` + `write_hr=.true.` 生成官方 hr.dat（recipe 见 `PT_INTERPOLATION_REPORT.md`）后重新比较 |
| W-08 | `pt_kramers_offgrid` | **WARN** | 离网格 δ_PT：median **0.433 meV**，max **4.54 meV**（k=[-0.1027,-0.3645,0.0031]） | 数值 | 严格 PT 对称晶体的离网格 Kramers 残差应恒为 0；这些 meV 值是 4×4×4 网格 Fourier 插值误差，**不是物理 PT 破缺**。不得据此声称简并 PASS | 官方 hr.dat/更密网格 + 离网格对比；在此之前保持 WARN |
| W-09 | `pt_quantitative_status` | **WARN** | 官方 hr.dat 缺失、`Pt_band.dat` 不能复现、离网格 Kramers 为数值残差 | 数据-环境 | 发布级别结论：Pt 目前**只能作软件接口/错误检测 benchmark**，不得用于定量 SOC/Floquet 材料结论 | 官方 hr.dat + 带结构复现 + 收敛证据齐备后才能升级 |
| W-10 | `pt_pst_gauge_invariance` | **WARN**（语义 = **NOT_APPLICABLE**） | 随机 U(2) 测试：PST 轴 dot 最低 **0.21**，collinearity 0.368–0.413 | 物理 | PT 对称非磁体系的 Kramers 简并使单带 PST 规范依赖：PST 轴/共线性随子空间旋转任意改变。Pt 没有可定义的持久自旋织构 | 对非简并/非中心体系才定义单带 PST；Pt 上不适用（代码已走 NOT_APPLICABLE 分支） |
| W-11 | `pt_band_interpolation_vs_wannier90` | **WARN** | index-aligned 593 点：median **54.1 meV**，max **1141.8 meV** | 数据-环境 | 与 `Pt_band.dat` 的差异无法用插值实现解释；参考 run 可能使用不同 Wannier 设置（如 `use_ws_distance`、解纠缠窗口）。旧 nearest-reference 指标（median 33.7 / max 773.8 meV）曾记为 PASS，已撤回（见 H-PASS-01） | 官方 hr.dat 或同设置参考带；在此之前保持 WARN |
| NA-01 | `pt_soc_splitting_not_applicable` | **PASS（语义 = NOT_APPLICABLE）** | "PT-symmetric non-magnetic static Pt: SOC spin splitting is NOT APPLICABLE (Kramers pairs must be degenerate). `soc` reports delta_PT and inter-doublet gaps, no spin-paired splitting." | 物理 | Pt 的 SOC 分裂量在该对称性下无定义；CLI `soc` 只输出 δ_PT 与双带间距（`outputs/Pt/soc_splitting.npz` 的 `note` 字段再次声明）。旧的 "48.8/62/126 meV SOC splitting" 是插值假象，已废弃 | 不适用；换非中心体系才讨论 SOC 自旋分裂 |

### 2.2 Benchmark 6：PbNV — `outputs/benchmarks/pbnv_spin/report.md` / `.json`

汇总：**PASS=10, WARN=0, FAIL=0, BLOCKED=0, SKIP=0**。无非 PASS 项。
注意其中两条是 informational（`pbnv_kramers_degeneracy`: PbNV 无空间反演，Kramers
简并本就不期待；`pbnv_adjacent_spin_antisymmetry`: 无 Kramers 配对），不要当成
"Kramers 简并已通过"的证据。

### 2.3 解析模型 benchmark（`outputs/benchmarks/`）

`rashba`（6 PASS）、`rashba_dresselhaus`（5 PASS）、`floquet_rashba`
（8 PASS）全部 PASS，无非 PASS 状态。

### 2.4 WCCl（外部数据，归档证据）— `outputs/benchmarks/wccl/report.md` / `.json`

汇总：**PASS=6, WARN=2, FAIL=0**。**注意**：`wccl` 不在
`floquet_soc/benchmarks/__init__.py::BENCHMARKS` 注册表中，仓库也没有 WCCl 的
config/数据，因此这是**历史运行证据，无法用当前仓库 CLI 复现**，只作记录。

| ID | 检查 | 状态 | 证据 | 类别 | 含义 / 影响 | 解决条件 |
|---|---|---|---|---|---|---|
| W-12 | `wccl_spn_quality` | **WARN** | spn eig range `[-1.431,+1.431]`，sum-rule `[0.636,6.132]` | 数据 | 赝波函数 .spn（非正交、非自旋完备）→ 自旋纹理仅半定量 | 提供正确的 Wannier 基 .spn |
| W-13 | `wccl_band_polarization_bound` | **WARN** | max per-band `\|<σ>\|=1.226 > 1` | 数据 | 同上，Pauli 界被违反 | 同 W-12 |

### 2.5 历史 / 验证副本（不得与当前证据混用）

| ID | 来源 | 状态 | 证据 | 说明 |
|---|---|---|---|---|
| H-W-01 | `outputs/PbNV_verify/validate.md/.json` | **WARN**（历史） | `spin_wannier_trs = 2.000` | 旧 TRS 判据 `max|eig S(-k)+eig S(k)|` 的错误结果；修正后当前值为 **0.0017 PASS**。该副本保留作审计轨迹，不得引用 2.000 |
| H-W-02 | `outputs/PbNV_verify/` | WARN（历史，重复项） | `spn_spin_sum_rule ∈ [0.719,3.000]`（旧文本） | 与 W-01 同源；汇总 24 PASS / 3 WARN |
| H-W-03 | `outputs/PbNV_verify/` | WARN（历史，重复项） | `spin_wannier_sum_rule ∈ [2.786,3.000]` | 与 W-02 同源 |
| H-W-04 | `outputs/Pt_verify/validate.md/.json` | WARN ×4（历史，重复项） | 同 W-03…W-06 | phase-2 中间验证副本（TRS 判据已修正，但 WS/ndegen 与 off-grid 检查尚未包含）：23 PASS / 4 WARN；冻结复跑为 **24 PASS / 4 WARN**（见 §1.2） |
| H-W-05 | `outputs/benchmarks_verify/pt/report.md/.json` | WARN（历史） | `pt_hr_reconstruction`（旧措辞）；该版本整体 9 PASS / 1 WARN | 更早版本的 benchmark；band 比较当时记为 PASS（见 H-PASS-01），已被最终 13 PASS / 5 WARN 取代 |
| H-PASS-01 | `outputs/benchmarks_verify/pt/report.md` | **撤回的 PASS**（历史） | `pt_band_interpolation_vs_wannier90` = PASS，nearest-reference median 33.7 meV / max 773.8 meV | 最终验收改用 index-aligned 指标：median 54.1 meV / max 1141.8 meV → **WARN（W-11）**。旧 PASS 已撤回，不得引用 |
| H-W-06 | `outputs/_tests/validate_pre_freeze_2026-10-23.md/.json`（冻结前副本） | WARN（历史） | `spn_pauli_eigenvalue_bound`：eig range `[-1.632,+1.632]`（违反 Pauli 界） | 旧（2026-01）run 的 .spn 数据质量问题；当前 PbNV 新 run 该检查已是 PASS |
| H-W-07 | 同上 | WARN（历史） | `spn_spin_sum_rule ∈ [0.541,7.962]` | 同上（旧 run 的赝波函数数据） |
| H-B-01 | 同上 | **BLOCKED**（历史） | `spin_operator_wannier_basis`：spn 为 192 带基而 H 为 156 带基，缺 chk/u.mat/u_dis.mat/amn/mmn | 2026-10-22 用户提供同 run 的 `wannier90.chk/eig` 后解除（见 STAGE_A_REPORT §10） |
| H-B-02 | 同上 | **BLOCKED**（历史） | `spin_texture_wannier`：缺失变换文件时拒绝猜测带序 | 同 H-B-01，已解除 |
| H-F-01 | `outputs/pt_audit.json`（phase-1 审计，无 token 字段） | 历史数值 | `kramers_offgrid` median **67 meV** / max **0.89 eV**；`kramers_mesh` max **7.9e-8 eV**；`inter_doublet_gap_meV` max 22.37 eV | 插值修复前的数值；现已被 phase-2 的 W-08 取代（median 0.43 meV / max 4.54 meV）。归档证据，不得引用 |
| H-F-02 | `outputs/Pt/validate_pre_freeze_2026-10-23.md/.json`、`outputs/Pt/pst_pre_freeze_2026-10-23.npz`（冻结前副本） | 历史文件 | Pt validate 旧文件：**23 PASS / 4 WARN**、旧 TRS 判据文本、`ndegen [1,1]`、`H_W-H_hr = 2.19e-14 eV`；Pt pst.npz 旧 run：base collinearity 0.3965、axis dot min 0.345 | 冻结前磁盘上的旧证据（validate 生成于 phase-2 之前；pst 为更早 run）。冻结复跑后当前文件为 24 PASS / 4 WARN（validate）与 base 0.4003 / axis dot 0.285（pst）；**结论不变**（Pt PST 仍 NOT_APPLICABLE）。旧文件保留供审计，不得引用旧计数 |

---

## 3. 收敛 / 未量化状态

| ID | 状态 | 来源 | 证据 | 类别 | 含义 / 影响 | 解决条件 |
|---|---|---|---|---|---|---|
| U-01 | **UNCONVERGED** | `PHOTON_TRUNCATION_REPORT.md`；`outputs/Pt/converge.npz`；`RELEASE_VALIDATION.md` §5/§6；CLI `converge` 输出 | `m_sp=16`：max per-state deficit **1.00**，2 个态 `D_i>1e-3`；`m_sp=32/48`：max deficit `3.9e-14`，0 个态；`accepted_m_sp = 32`；`n_deficit_gt_1e3_m = [2, 0, 0]` | 数值 | 光子窗截断：m_sp=16 时目标能窗内 2 个态几乎丢失全部光子权重。任何 `D_i>1e-3` 的态**不得**用于光子分辨谱/自旋谱权重/PST 类求和；`m_sp=32` 只是当前 Pt/PbNV 设置的收敛选择，不是普适默认 | 每个材料重跑 deficit 检查；用 `accepted_m_sp` 或更大窗口；受影响态排除 |
| U-02 | 未量化（光子分辨） | `outputs/PbNV/converge.npz`（旧格式） | 仅含 `status=['PASS']`，**没有** `w_deficit_*` / `accepted_m_sp` 字段；`outputs/PbNV/` 无 `spectrum.npz` | 数据/数值 | PbNV 的"收敛 PASS"只覆盖 n_step 与旧版 m_sp 全局量；PbNV 的光子分辨光谱未做过 deficit 检查 → 记为**未量化**，不等于已收敛 | 用当前 `converge` 重新跑 PbNV 生成 deficit 字段，并按 U-01 判据处理 |
| U-03 | 近似（数值/定义） | `README.md` 已知限制；`outputs/PbNV/spin_floquet.npz`；`STAGE_B_REPORT.md` §5 | overlap 匹配纯度 median **0.996**，`>0.8` 占 **88.9%**（即 11.1% 的 k 点纯度 <0.8） | 数值 | 近简并 Floquet 态在 Schur 分解下会混合；低纯度点的单态结论不可靠 | 做简并子空间投影，或剔除低纯度点后再统计 |
| U-04 | **WARN（近似 Ising-like）** | `PBNV_PST_REPORT.md`；`README.md`；`RELEASE_VALIDATION.md` §5；`outputs/PbNV/pst.npz` | M=(0.9570,0.0401,0.0029)，collinearity 0.9570；gauge axis dot min 0.998、collinearity spread **0.038**；bootstrap 0.9421–0.9730；near-Γ 0.865 vs far 0.964 | 物理/规范 | PbNV PST 轴稳定但共线性有 3.8% 规范涨落（近 TRIM 简并）；只能作为**近似** 2D Ising-like 参考，不是精确 PST 常数 | 更大 k 网格、去简并子空间定义、明确规范后重新评估 |

---

## 4. 接口 / 输入状态

| ID | 状态 | 来源 | 证据 | 类别 | 含义 / 影响 | 解决条件 |
|---|---|---|---|---|---|---|
| NI-01 | **NEEDS_INPUT** | `STAGE_C_REPORT.md` §2/§4；`README.md`（Stage C）；`RELEASE_VALIDATION.md` §7；`NEW_MATERIAL_GUIDE.md` §5；`observables/soc_splitting.py::soc_off_difference` | SOC on/off 两套 run 的 win/hr/centres 哈希不一致或文件缺失时返回 NEEDS_INPUT，拒绝"伪减法" | 接口约定 | `H_SOC = H(λ=1)-H(λ=0)` 只在同基组/同规范下有意义；缺输入时框架显式拒绝而不是给错误数字 | 提供同基组、同规范的 SOC-off run（win/hr/centres 一致） |
| B-01 | **BLOCKED** | `PT_HAMILTONIAN_REPORT.md`；`PT_INTERPOLATION_REPORT.md`；`PT_PHASE2_REPORT.md`（Remaining 1）；`PT_ACCEPTANCE_REPORT.md`（Blocked/external）；`RELEASE_VALIDATION.md` §11.1；`README.md` | 官方 Pt `wannier90_hr.dat` 不存在：本机无 wannier90 可执行文件，公开 tutorial 数据也没有 hr.dat | 数据-环境 | 官方插值对比、离网格 Kramers/带结构定量验证全部依赖该文件；Pt 因此停留在 interface benchmark（W-09） | 在 Linux/HPC 上执行 `restart = plot` + `write_hr = .true.` 生成 hr.dat（recipe 见 `PT_INTERPOLATION_REPORT.md`），再跑 Hamiltonian/带结构/Kramers/PST 对比 |
| B-02 | **BLOCKED（代码路径，未触发）** | `floquet_soc/cli.py::_pst_not_applicable` | Pt `pst` 分支在 DFT k-mesh 缺失时打印 `BLOCKED: no DFT k-mesh available for the gauge test` 并返回 2 | 接口约定 | 保留显式拒绝语义；当前 Pt 数据含 win kpoints，故未触发 | 提供含 k 网格的 win/chk |
| B-03 | **BLOCKED（历史，已解除）** | `STAGE_A_REPORT.md` §4/§8 | `.spn`（192 带基）→ Wannier（156 带基）旋转缺 `wannier90.chk`/`u.mat`/`u_dis.mat`/`amn`/`mmn` | 数据-环境 | 2026-10-22 同 run 的 chk/eig 到位后解除；当前 validate 中 `spin_operator_wannier_basis` 为 PASS | 已解决，保留历史记录 |
| NV-01 | **NOT_VALIDATED** | `NEW_MATERIAL_GUIDE.md` §1；`RELEASE_VALIDATION.md` §6；§1.1 N-01 | 有位移 k 网格的 `chk+eig` alias 回退：网格值精确、离网格物理 `NOT_VALIDATED`（自旋算符离网格 Hermiticity 检查直接标 NOT_VALIDATED） | 数值/后端 | 任何 shifted-mesh 材料都不得用离网格结果做定量结论 | 换无位移网格（WS/ndegen 重建）、提供官方 hr.dat，或对具体材料做离网格实测验证 |

---

## 5. 报告层 audit trail（`PT_*` / `STAGE_*` / `PBNV_*` / `PHOTON_*` / 根目录 md）

> 这些是日期化历史报告，保留原文与原始数字；其中若干计数已被最终验收取代
> （本表逐条注明）。最终结论以 §1–§4 与 `RELEASE_VALIDATION.md` 为准。

| 报告 | 状态 / 结论 | 说明 |
|---|---|---|
| `PT_ACCEPTANCE_REPORT.md` | Level 1 **PASS**；Level 2 **PARTIAL**（PbNV PASS / Pt WARN）；Level 3 **PASS**；Level 4 **NOT YET**；Pt Benchmark **13 PASS / 5 WARN**；官方 hr.dat **BLOCKED** | 发布级最终结论；Level 4 "NOT YET" 是发布定性（Pt 定量验证未完成），不是测试失败 |
| `PT_PHASE2_REPORT.md` | 历史 "PASS=14 / WARN=3"；**BLOCKED** hr.dat；Remaining limitations 1–5 | 数字为 phase-2 中间结果；冻结版最终 13 PASS / 5 WARN（报告内已加冻结说明）。限制 1–5 仍然有效（hr.dat、参考带设置、离网格 Kramers、PbNV PST 规范、m_sp 每材料重验） |
| `PT_PHYSICS_AUDIT.md` | 历史 "PASS=14 / WARN=3"；**BLOCKED** hr.dat；SOC/PST **NOT_APPLICABLE**；旧 "48.8/62/126 meV" 已废弃 | 报告内已加冻结说明（最终 13/5）；审计问答 1–6 保留 |
| `PT_KRAMERS_REPORT.md` | 修复前 **FAIL**（median 67 meV / max 0.89 eV）→ 修复后离网格 **WARN 级**（median 0.43 meV / max 4.54 meV）；DFT 网格 PASS | FAIL 是历史数值；当前状态见 W-08。容差 rationale（DFT 1e-6 eV；离网格 median<1 meV、max<10 meV 因 4×4×4 截断）保留 |
| `PT_GAUGE_INVARIANCE_REPORT.md` | 单带 PST 判定 **FAIL**（历史）→ 现在 **NOT_APPLICABLE**；PbNV PST 规范涨落 | Pt 单带 PST 轴 dot 最低 0.345、invariant spin mean 0.0024；PbNV 轴稳定 0.998 / 共线性 spread 0.038（见 U-04） |
| `PT_HAMILTONIAN_REPORT.md` | WS/ndegen 修复说明；官方 hr.dat **BLOCKED**；后端标签 "Grid-consistent WS/ndegen reconstruction (Hermitian); official-hr comparison pending" | 见 W-07；修复前非 Hermiticity 1.27 eV 为历史 |
| `PT_INTERPOLATION_REPORT.md` | 官方 hr.dat **BLOCKED**；各插值约定对比表中所有重构都无法复现 `Pt_band.dat`；旧 alias 的"更接近参考"是假象（非 Hermitian/破坏 PT） | 见 W-07/W-11；包含 Linux/HPC 生成 hr.dat 的 recipe |
| `PBNV_PST_REPORT.md` | PbNV PST **WARN（近似 Ising-like）** | 见 U-04；数据文件路径已修正为 `outputs/PbNV/pst.npz` |
| `PHOTON_TRUNCATION_REPORT.md` | `m_sp=16` → 2 个态 **UNCONVERGED**（D_i≈1.0）；`m_sp=32/48` 收敛；production `[n_step=160, m_sp=32]` | 见 U-01；包含 npz 字段清单（`w_deficit_max_m`、`n_deficit_gt_1e3`、`accepted_m_sp`） |
| `STAGE_A_REPORT.md` | 历史 **BLOCKED**（自旋算符/纹理，§4/§8）已被 §10 取代；`.spn` 数据质量 **WARN**；历史测试 18 项 | BLOCKED 记录保留；当前状态 PASS（H-B-01/H-B-02 已解除）。历史 "25 PASS" 是加 off-grid 检查之前的计数，现为 26 |
| `STAGE_B_REPORT.md` | 全部解析检查 PASS；数据近似 ~1e-3；micromotion 差异 8.9e-3（如实报告）；匹配纯度 caveat；25 k 演示网格 | 自旋结论按半定量使用（W-01/W-02）；低纯度点见 U-03；25 k 仅为演示规模 |
| `STAGE_C_REPORT.md` | 解析验证 PASS；SOC on/off **NEEDS_INPUT**；PbNV 无 Kramers 简并（相关检查视为 N/A）；PbNV 分裂/PST 半定量 | 见 NI-01、NA（Kramers）；历史 PST 数字 0.945 与最终 0.9570 因网格不同，以最终输出为准 |
| `STAGE_D_REPORT.md` | Pt validate 23/4；Pt benchmark 历史 9/1（nearest 指标）；PbNV 历史 25 PASS；`static_gap` WARN；m_sp 收敛说明 | 历史计数已被冻结复跑的 **Pt validate 24/4**、**PbNV validate 26/2**、**Pt benchmark 13/5** 取代（报告内已加冻结说明）；其余限制仍有效 |
| `RELEASE_VALIDATION.md` | §5 WARN 清单、§6 未定量验证（含 `D_i>1e-3` UNCONVERGED 规则）、§7 token policy、§11 已知未解决限制 1–5 | 发布主文档；本登记表与它互补（本表逐条给出来源与解决条件） |
| `NEW_MATERIAL_GUIDE.md` | 新材料的硬性 gate：0 FAIL、WARN 必读、shifted-mesh `NOT_VALIDATED`、PST 规范测试（axis dot<0.99 或 spread>1e-2 → WARN/N/A）、光子 deficit `D_i<1e-3` 否则 UNCONVERGED、`NEEDS_INPUT` 时不得伪减法 | 后续研究必须遵守；所有非 PASS token 必须继续保留 |
| `README.md` | Benchmark 矩阵 + 已知限制（Pt 插值/金属、近简并处理、bcc Fe 候选） | 口径与 `RELEASE_VALIDATION.md` 一致 |
| `benchmarks_data/Pt/README.md` | Pt 数据来源/授权 + hr.dat 缺失说明；参考带偏差已更新为最终 index-aligned 数字 | 数据再分发仅用于框架验证，授权见上游仓库 |

---

## 6. npz / 元数据中内嵌的状态

| 文件 | 字段 | 值 / 含义 |
|---|---|---|
| `outputs/Pt/pst.npz` | `note` | `PST NOT_APPLICABLE: gauge-dependent single-band spins in Kramers-degenerate PT system` |
| `outputs/Pt/pst.npz` | `system_type` | `pt_nonmagnetic`；含 `gauge_collinearity`、`gauge_axis`、`gauge_invariant_spin` 证据 |
| `outputs/Pt/pst_pre_freeze_2026-10-23.npz` | 旧 run 证据 | 冻结前 `pst.npz` 原样副本（base collinearity 0.3965 / axis dot 0.345）；结论同为 NOT_APPLICABLE，见 H-F-02 |
| `outputs/Pt/soc_splitting.npz` | `note` | `Kramers residual / inter-doublet gaps; SOC spin splitting is not defined for a PT-symmetric non-magnetic system`；含 `kramers_delta`、`inter_doublet_gaps` |
| `outputs/Pt/converge.npz` | `n_deficit_gt_1e3_m`, `w_deficit_max_m`, `accepted_m_sp`, `status` | `[2,0,0]`；`[1.0, 3.9e-14, 3.9e-14]`；`32`；`PASS`（PASS 指整体判据，m_sp=16 的 2 个态仍为 UNCONVERGED，见 U-01） |
| `outputs/PbNV/converge.npz` | `status` | `PASS`（旧格式，无 deficit 字段 → 光子分辨未量化，见 U-02） |
| `outputs/PbNV/spin_floquet.npz` | `purity` | overlap 匹配纯度（median 0.996；低纯度点见 U-03） |
| `outputs/PbNV/pst.npz` | `gauge_coll_spread`, `gauge_axis_dot_min`, `gauge_invariant_mean_norm` | 0.038 / 0.998 / 0.035（近似 Ising-like，见 U-04） |
| `outputs/PbNV/validate.json`、`outputs/Pt/validate.json` | `checks[].status` | 与 §1 一致；`summary` 保留 WARN 计数 |

---

## 7. 状态统计与测试 skip

冻结验收（`RELEASE_TEST_RESULTS.txt`）：**Ran 74 tests, OK (skipped=2)**，
0 failure、0 error。skip 项：

| ID | 状态 | 测试 | 原因 | 类别 | 解决条件 |
|---|---|---|---|---|---|
| SKIP-01 | **SKIP** | `test_spin_wannier.TestCheckpoint.test_read_wannier_spn_dispatch` | 依赖 `C:\Users\15756\Desktop\wannier90.spn`（formatted fixture），不在快照/CI 环境 | 数据-环境 | 把该 fixture 纳入仓库或生成脚本；注意该测试对 config `.spn` 的旧断言期待 `vasp-binary-spn`，而当前 PbNV spn 为 formatted，恢复 fixture 前需同步更新断言 |
| SKIP-02 | **SKIP** | `test_stage_b.TestPbNVRegression`（setUpClass） | 跨引擎参考 `outputs/_tests_b/old_engine_ref.npz` 缺失（需 `gen_old_engine_ref.py` 生成，脚本未随仓库提供） | 数据-环境 | 在具备旧引擎参考的机器上重新生成参考 npz 后解 skip |

登记条目统计（不重复计 §5 报告中与 §1–§4 同源的引用）：

| token | 当前证据 | 历史/副本证据 | 合计登记条目 |
|---|---|---|---|
| WARN | 14（W-01…W-13 + U-04） | 7（H-W-01…H-W-07） | 21 |
| NOT_APPLICABLE | 2（NA-01、W-10 语义） | 0 | 2 |
| NOT_VALIDATED | 2（N-01、NV-01） | 0 | 2 |
| UNCONVERGED | 1（U-01） | 0 | 1 |
| BLOCKED | 2（B-01、B-02） | 1（B-03） | 3 |
| NEEDS_INPUT | 1（NI-01） | 0 | 1 |
| FAIL | 0 | 2（PT_KRAMERS 修复前、PT_GAUGE 单带 PST） | 2 |
| SKIP | 2（SKIP-01/02） | 0 | 2 |
| 撤回的 PASS（单独记录） | 0 | 1（H-PASS-01） | 1 |
| **合计** | **24** | **11** | **35** |

说明：W-10 作为一行 WARN 记录，同时承载 NOT_APPLICABLE 语义（在两个 token 列
各计一次）；`outputs/Pt_verify` 的 4 个 WARN 与当前 Pt validate 完全同源，计入
H-W-04；`benchmarks_verify/pt` 的 1 个 WARN 计入 H-W-05。

---

## 8. 使用规则（冻结声明）

1. 引用任何数值前先读本表；WARN 不得表述为 PASS，"未量化/近似"不得表述为"已验证"。
2. 历史证据（`*_verify/`、`_tests*/validate_pre_freeze_*`、`pt_audit.json`、
   日期化 `PT_*`/`STAGE_*` 报告）只作审计轨迹；与当前证据冲突时以 §1–§4 为准。
3. 已被取代/撤回的旧结论（H-W-01 的 TRS=2.000、H-PASS-01 的带比较 PASS、
   `pt_audit.json` 的 67 meV 离网格 Kramers、STAGE_A 的 BLOCKED）**禁止再引用**。
4. 新增材料必须按 `NEW_MATERIAL_GUIDE.md` 执行；产生的非 PASS 状态必须继续以
   WARN/NOT_APPLICABLE/UNCONVERGED/BLOCKED/NEEDS_INPUT token 保留，不得用全局
   PASS 计数掩盖。
5. 本快照冻结后如需修改代码/数据，必须新开分支并打新 tag（如
   `v0.1.1-rc1`），不得在 `v0.1.0-rc1` 上就地修改。
