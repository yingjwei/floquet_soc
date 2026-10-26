# Stage D 报告 — 多材料 benchmark、统一 CLI/配置、扫描与收敛

**状态：PASS。** Stage D 的交付内容：

1. `scan` / `converge` 工作流（CLI：`python floquet.py scan|converge`）；
2. **Benchmark 7：fcc Pt**（18 个 Wannier 函数、SOC、真实 DFT 数据；来自公开
   Wannier90 tutorial 仓库）；
3. 新输入路径 `chk+eig`（无需 `wannier90_hr.dat`），以及该路径的验证与
   文档化限制；
4. 跨材料矩阵：2 个解析模型 + PbNV（B6）+ Pt（B7），轨道数 2/156/18、
   三种输入格式、不同 SOC 对称性。

## 1. Benchmark 7 数据来源（GitHub 公开数据）

- 仓库：`https://github.com/wannier-developers/wannier-tutorials`
  （branch `master`）
- 目录：`2022_05_Trieste/DAY2_WBERRI/advanced/spin_hall/data_Pt`
- 本地副本：`benchmarks_data/Pt/`（约 6.7 MB），来源与授权说明见
  `benchmarks_data/Pt/README.md`
- 文件：`Pt.win / Pt.chk / Pt.eig / Pt.spn / Pt.wout / Pt.nnkp /
  Pt_band.dat / Pt_band.kpt / scf.out`
- 体系：fcc Pt，非共线 spinor Wannier 函数（SOC，spin Hall tutorial），
  `num_bands=24`、`num_wann=18`、`mp_grid=4×4×4`（64 个 k 点），
  `num_electrons=10`（从 `scf.out` 读取，写在配置中，**框架不硬编码**）

> 该数据集**不含** `wannier90_hr.dat`。框架新增 `WannierModel.from_chk_eig`
> 输入路径：由 `H_W(k)=v†diag(eig)v` 在 DFT k 网格上做 Fourier 变换得到
> H(R)。该路径在网格点上**精确**（见 §2），网格点之间是 Fourier 插值；
> 与 wannier90 官方 Wannier 插值的差异在 §2 中量化并作为已知限制报告。

## 2. Benchmark 7 验证结果（PASS/WARN）

`python floquet.py validate --config configs/pt.json`：
**23 PASS / 4 WARN / 0 FAIL，exit 0**。关键数字：

| 检查 | 结果 |
|---|---|
| chk/eig 重建在 64 个网格 k 上 vs `v†diag(eig)v` | **2.19e-14 eV** |
| `U` 幺正性 / `v` 等距 / 晶格 | 1.35e-13 / 1.35e-13 / 9.1e-18 |
| 自旋矩阵（QE `pw2wannier90` 产物）本征值 | **[-1.0000, +1.0000]** |
| 自旋求和规则 `eig(ΣS²)` | **[0.333, 3.000]** |
| TRS：`spec S(-k) = -reverse spec S(k)` | **0.0000** |
| `spin_wannier_interpolation` | 1.3e-14 |
| 与 wannier90 插值带（`Pt_band.dat`，593 点） | index-aligned 中位 **54.1 meV**，最大 **1141.8 meV**（hr-less 重建的 Fourier 插值限制；参考 run 设置可能不同，WARN 说明） |

`python floquet.py benchmark --which pt`：**9 PASS / 1 WARN / 0 FAIL**（该数字为
当时 nearest-reference 比较；最终冻结版本的 index-aligned 比较把它降级为
**13 PASS / 5 WARN**，见 `RELEASE_VALIDATION.md` §7 与
`outputs/benchmarks/pt/report.md`），
报告与图：`outputs/benchmarks/pt/`。

## 3. Stage C/D 在 Pt 上的运行结果（第二材料全链路）

- `scan`（ħω=0.5/1.0/1.5 eV，γ=0/0.02/0.04；n_step=160，m_sp=48）：
  光子权重归一化误差 ≤2.2e-14；准能量展宽随 ħω 增大（0.376→1.304 eV 的
  折叠区间内变化）；周期平均自旋 |S̄|≈0.744，随 γ 轻微下降；数据：
  `outputs/Pt/scan.npz`、`scan.png`。
- `converge`：PASS。
  - 生产设置（n_step=160, m_sp=32）对参考（n_step=240, m_sp=48）：
    `dĒbar=7.3e-7 eV`，`dS̄=2.7e-6`；
  - 不变量（ΣĒbar、ΣS̄、Σw）偏差 **3.4e-12**；
  - 数据：`outputs/Pt/converge.npz`、`converge.png`。
- `soc`（**2026-10-23 修正**）：Pt 是 PT 对称非磁体系，SOC splitting
  **NOT_APPLICABLE**；旧的 48.8/62.0/126.1 meV 是插值假劈裂与双重态间距
  （见 `PT_PHYSICS_AUDIT.md`）。现在 `soc` 只报 delta_PT 与 inter-doublet
  gaps：路径上 delta_PT 最大 4.37 meV、中位 0.27 meV（WS/ndegen 插值后；
  修正前离网格中位 67 meV、最大 0.89 eV）。输出`soc_splitting.npz`/
  `soc_splitting.png`。## 4. 一个重要的判据修正：TRS 本征值匹配

Stage A/B 的 TRS 质量检查原来写成
`max|eig S(-k) + eig S(k)|`，这是**错误的**：正确条件是
`spec S(-k) = -spec S(k)`，即排序后的本征值应满足
`sorted(S(-k)) = -reverse(sorted(S(k)))`。旧公式在“每个 k 的本征谱本身
关于 0 对称”（例如中心对称体系）时会把正确数据误判为破坏 TRS。

修正后：

| 材料 | 旧（错误）TRS 判据 | 修正后 TRS | validate 汇总 |
|---|---|---|---|
| PbNV（新 run，VASP LWRITE_SPN） | 2.000（误报） | **0.002** | **25 PASS / 2 WARN / 0 FAIL** |
| Pt（QE pw2wannier90） | 2.000（误报） | **0.0000** | 23 PASS / 4 WARN / 0 FAIL |

因此之前“PbNV `.spn` 有 ~1% TRS 残差”的结论有一部分来自判据 bug；
新 run 的 `.spn`（Pauli 界 100%、求和规则 [2.786,3.000]、TRS 2e-3）是
质量合格的 VASP LWRITE_SPN 数据。PbNV 剩余的 WARN 只有求和规则下界
（2.786）与深能带/窗口相关项，需在 PbNV 具体问题的讨论中按需处理。

## 5. 跨材料矩阵（Stage D 验收）

| 材料 | 轨道数 | SOC/对称性 | 输入 | 链路 |
|---|---|---|---|---|
| Rashba 解析（B2） | 2 | 面内 Rashba | 解析 | A/B/C 全 PASS |
| Rashba+Dresselhaus（B3） | 2 | PST 面内 | 解析 | A/B/C 全 PASS |
| 驱动 Rashba（B4，解析） | 2 | 周期驱动 | 解析 | Stage B 8/8 PASS |
| **PbNV（B6）** | **156** | 非中心、Ising 型 SOC | win+hr+chk+eig+spn | A/B/C PASS，validate 25/2/0 |
| **Pt（B7）** | **18** | 中心、强 SOC 金属 | win+chk+eig+spn（无 hr） | A/C/D PASS，validate 23/4/0 |

## 6. 已知限制与下一步

1. Pt 数据无 `hr.dat`：网格点之间是 Fourier 插值（与 `Pt_band.dat` 的
   index-aligned 偏差中位 54.1 meV / 最大 1.14 eV，参考 run 的设置可能不同）。
   对“框架验证”足够；若要做 Pt 的定量能带物理，应补 `hr.dat` 或 wannier90
   输出。
2. Pt 是金属：`static_gap` 如实报告为 WARN（无隙），SOC 窗口在 `soc`
   命令中自动切换为费米面附近窗口；这些都在报告中标注。
3. 高对称点路径写在配置里（`kpath`），框架不硬编码材料。
4. 可选扩展：同一 tutorial 仓库还有 **bcc Fe**（18 个 spinor Wannier，
   `.chk + .eig + .spn`，磁性 SOC，`mp_grid=8×8×8`）可作磁性 benchmark 8；
   Si/GeTe 等无 SOC 数据也可用于纯 Hamiltonian 回归。

## 7. 复现

```bash
python floquet.py validate  --config configs/pt.json
python floquet.py benchmark --which pt
python floquet.py scan      --config configs/pt.json
python floquet.py converge  --config configs/pt.json
python floquet.py soc       --config configs/pt.json --nk 10
python floquet.py pst       --config configs/pt.json
python -m unittest tests.test_stage_d -v
```

## 8. 冻结说明（2026-10-23，v0.1.0-rc1）

本报告的计数为 Stage D 当时状态；冻结前用 README/RELEASE_VALIDATION 中的命令
**复跑**后，当前代码的实际输出为：

| 项 | 本报告历史计数 | 冻结复跑（当前代码） | 差异原因 |
|---|---|---|---|
| PbNV validate | 25 PASS / 2 WARN | **26 PASS / 2 WARN** | 后加 `spin_wannier_hermiticity_offgrid` 检查 |
| Pt validate | 23 PASS / 4 WARN | **24 PASS / 4 WARN** | 同上（WS/ndegen 后端 off-grid 自旋厄米性 PASS） |
| Pt benchmark 7 | 9 PASS / 1 WARN（nearest 参考带指标） | **13 PASS / 5 WARN**（index-aligned） | 最终验收改用 index-aligned 比较，并把 `pt_band_interpolation_vs_wannier90`、`pt_kramers_offgrid`、`pt_quantitative_status` 记为 WARN |

历史计数保留不改；以 `outputs/*/validate.*`、`outputs/benchmarks/pt/report.*`
（冻结复跑/最终验收生成）、`RELEASE_VALIDATION.md` 与 `STATUS_REGISTRY.md` 为准。
当时全套测试为 56 项规模；冻结版为 **74 项（2 项数据依赖 skip）**，见
`RELEASE_TEST_RESULTS.txt`。
