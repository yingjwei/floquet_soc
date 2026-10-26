# Stage C 报告 — SOC 相关能带分裂、有效项拟合、PST

**状态：PASS（解析 Rashba/RD 全项验证 + PbNV 实际数据应用）。**
实现：`floquet_soc/observables/soc_splitting.py`、`observables/pst.py`（扩展），
CLI `soc` / `pst`；测试 `tests/test_stage_c.py`（10 项，1 skip）。

## 1. 物理定义（避免误判）

- **相邻能带分裂** `Δ_n(k)=E_{n+1}(k)−E_n(k)`：直接、无模型假设；
- **自旋配对分裂**：在给定能量窗口内，把“能量最近 + 自旋反平行”的相邻能带配对，
  规则显式记录（相邻带 j=i+1、能量差 ≤ `--max-gap`、反平行度作为输出）；
- 报告一律称为 **spin-related spectral splitting**，不直接叫“SOC 强度”，
  因为分裂可能来自 SOC、交换场或其它对称性破缺；
- `Δsplit^F`（驱动下的准能量分裂）需要与 Stage B 的 overlap matching 联用，
  属于 Stage C/D 的下一步（当前框架已具备 `S_bar` 匹配数据）。

## 2. 解析验证（PASS/FAIL）

| 检查 | 基准 | 结果 |
|---|---|---|
| 相邻分裂 = 2αk | Rashba α=0.1 | max dev **<1e-12** |
| 自旋配对分裂 gap = 2α|k|，反平行度 | Rashba | max dev **<1e-10**，min 反平行 **>0.99** |
| 有效项拟合 | `H = c k²I + α(kyσx − kxσy)`，基 `[k², ky, −kx]` | 系数 **0.5000, 0.1000, 0.1000**，rms **<1e-12** |
| PST 本征值 | Rashba | **(0.5, 0.5, 0)**，dev <1e-8 |
| PST 共线性/主轴 | RD α=β | **(1,0,0)**，主轴 **(1,−1,0)/√2**，角误差 <1e-8 |
| 非共线判据 | RD α≠β | collinearity **0.65 < 1** |
| SOC on/off 门控 | 两套不同基组 | 返回 **NEEDS_INPUT**（拒绝直接相减） |

## 3. PbNV 应用（Benchmark 6，156 轨道，112 电子）

窗口 `[−1.08, 1.01] eV`（VBM−0.3 到 CBM+0.3），路径 Γ–X–S–Y–Γ（81 点）：

**自旋配对分裂（相邻能带，`--max-gap 0.2 eV`；本表为 v0.1.0-rc1 贪心配对口径的历史值）：**

| 能带对 | median | max |
|---|---|---|
| 110–111（VBM） | 0.53 meV | **29.98 meV** |
| 112–113（CBM） | 0.38 meV | **103.28 meV** |
| 114–115 | 0.45 meV | **71.98 meV** |

> **注（v0.1.3-rc1 起）**：`soc` 命令只报 **gap-edge 两对**的直接能差
> （VBM-1/VBM = 110-111、CBM/CBM+1 = 112-113），不再对 SOC 窗口内所有相邻带做
> 贪心配对。当前口径下 PbNV 的值：110-111 max **29.98** / median **0.51** meV；
> 112-113 max **107.66** / median **0.58** meV（81 k）。上表 103.28/71.98 为历史
> 口径，保留不改；当前值以 `outputs/PbNV/soc_splitting.npz` 的
> `edge_pair_rows`/`spin_pair_rows` 为准。

空间分布高度各向异性：Γ–X–S 段接近 0，分裂集中在 **Y–Γ 附近**（与静态自旋纹理的
面内过渡带、Floquet 自旋改变的热点位置一致）。

**PST（能带 110–113，k 网格 11×11）：**

- `M` 本征值 = **(0.9450, 0.0547, 0.0003)**，共线性 0.945；
- 主轴 = **(0.0164, 0.0183, 0.9997) ≈ z**，平均极化 0.916，484 个态；
- 解释：PbNV 近带边是 **面外（Ising 型）近共线自旋纹理**，不是 Rashba 面内织构；
  这与静态 `Sz≈±1` 分区图一致。共线性 0.945 而非 1，说明 ky≈0 过渡带
  存在明显的面内自旋分量。

## 4. SOC on/off 策略

直接 `H_SOC=H(λ=1)−H(λ=0)` 只在两套 run 共享 Wannier 基组/规范时才有意义。
`floquet_soc.observables.soc_splitting.soc_off_difference` 对 win/hr/centres
做哈希比对；不同基组或文件缺失时返回 **NEEDS_INPUT** 并说明原因，不做“伪减法”。

## 5. 注意事项

1. PbNV 无空间反演：Kramers 简并相关检查不适用（输出应视为 N/A），
   带分辨自旋纹理本身是物理的（DFT 能级无 Kramers 简并，最近邻带间距中位 31 meV）。
2. 新 `.spn` 的 TRS/求和规则残差在 1e-2 量级（VASP 赝波函数近似），
   因此分裂/PST 结果按半定量使用；若要定量 SOC 常数，建议：
   - 用 VASP 原生 `LWRITE_SPN`（当前文件即此路径）+ 更大 k 网格；或
   - 对给定低能模型做 `fit_linear_model` 拟合并报告残差（框架已提供）。
3. 有效项拟合是**用户指定基函数与算符**的最小二乘，不预设 Rashba/Dresselhaus
   形式；对任意材料需先给低能模型，再报告拟合残差与条件数。

## 6. 复现

```bash
python floquet.py soc --config configs/pbnv.json --nk 20
python floquet.py pst --config configs/pbnv.json
python -m unittest tests.test_stage_c -v
```

输出：`outputs/PbNV/soc_splitting.npz`、`soc_splitting.png`、
`soc_spin_splitting.png`、`pst.npz`、`pst.png`。
