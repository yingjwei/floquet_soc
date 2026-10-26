# 阶段 B 补充验收（子空间自旋 + Benchmark 1/5）

## 1. PbNV / Pt 的简并子空间自旋

单带自旋在近简并点依赖 gauge，不能直接用于 PbNV/Pt。改用投影权重平均：

\[
\bar{\mathbf S}_{\rm sub}=
\frac{\sum_i w_i \mathbf S_i}{\sum_i w_i},\qquad
w_i=\|P_{\rm sub}u_i\|^2
\]

| 材料/点 | 子空间 | `Ebar` (eV) | `S_bar` (Pauli) | overlap | `w_main` | 子空间 gap |
|---|---|---|---|---|---|---|
| PbNV X | valence pair | −1.249 | (0.004, 0.006, 0.004) | 1.00 | 1.000 | 0.39 meV |
| PbNV X | conduction pair | 1.077 | (0.045, 0.136, 0.005) | 1.00 | 1.000 | 0.13 meV |
| Pt Γ | valence pair | 16.549 | (−0.023, 0.029, −0.054) | 0.89 | 0.000 | 0 meV |
| Pt Γ | conduction pair | 16.549 | (0.010, −0.010, 0.195) | 0.89 | 0.000 | 0 meV |

解释：

- PbNV 的价带对平均自旋几乎为零，是因为近简并对内 `S_z` 反号、相互抵消；单带看到的 ±`S_z` 是规范选择，不能直接当“Ising PST 的带自旋”报告。
- Pt 的 `w_main=0` 说明所选 Floquet 态在目标子空间没有主光子权重（金属/简并），只能用于 solver/unitarity 测试，不能作 SOC 物理结论。
- 这两项仍为 **WARN**，但已从“无法处理”变为“有定义、可复现的子空间量”。

## 2. Benchmark 1：无 SOC 自旋简并

`python floquet.py benchmark --which spin_degenerate`：

| 检查 | 结果 |
|---|---|
| 静态自旋简并 | 0.00e+00 eV |
| 驱动后准能级简并 | 0.00e+00 eV |
| 无假自旋信号（子空间平均） | 0.00e+00 |

**PASS 3 / FAIL 0**。

## 3. Benchmark 5：Rashba ferromagnet

`H = c k^2 + α(k_yσ_x−k_xσ_y) + Mσ_z`

| 检查 | 结果 |
|---|---|
| Γ 点交换劈裂 = 2M | 0.400000 eV = 0.4 eV |
| Γ 点自旋沿 ±z | \|⟨S_z⟩\|=0.5, max\|S_xy\|=0 |
| 有限 k 出现 SOC 贡献 | ΔE(0.1)=0.4005 eV，偏离 2M 5e−4 eV |
| γ=0 恢复静态（模 ħω） | 2.2e−16 eV |

**PASS 4 / FAIL 0**。

## 4. 当前判定

- Benchmark 1/5：**PASS**；
- PbNV/Pt 带边自旋：**WARN（子空间量已可算，但单带 Ising 标签不可用）**；
- MoSSe：**PASS**（带边劈裂足够，单带/子空间一致）。

复现：

```powershell
python floquet.py benchmark --which spin_degenerate
python floquet.py benchmark --which magnetic_soc
python floquet.py floquet_points --config configs\pbnv.json --points "0.5,0,0" --gamma 0 0.02 --pol circ+ --n-step 40 --m-sp 8 --subspace
python floquet.py floquet_points --config configs\pt.json --points "0,0,0" --gamma 0 0.02 --pol circ+ --n-step 80 --m-sp 16 --subspace
```
