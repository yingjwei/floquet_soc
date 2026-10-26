# Stage B 报告 — Floquet 自旋纹理与自旋分辨谱

**状态：PASS（8/8 解析恒等式/收敛性检查 + PbNV 实际数据应用）。**
实现文件：`floquet_soc/solvers/{evolution,sambe,effective}.py`、
`observables/{spectral,floquet_spin}.py`；测试 `tests/test_stage_b.py`（13 项 OK）。

## 1. 求解器约定（已验证）

- `U(T)|ψ⟩ = e^{-iεT}|ψ⟩`；物理态展开 `ψ(t)=Σ_m e^{-i(ε-mħω)t}u^{(m)}`；
- 因此频谱能量 `E_{αm}=ε_α-mħω`，光子权重 `w_{αm}=‖u^{(m)}‖²`；
- 周期平均自旋 `S̄ = (1/T)∫₀ᵀ⟨ψ(t)|S|ψ(t)⟩dt`（不是 t=0 快照）；
- 所有实现都与以下恒等式对拍：`Σ_m w=1`、`Ē=ε-ħωΣ_m m w`、
  `Σ_m s_{αm}=S̄_α`、`∫A^S dE=Σ_α S̄_α`。

## 2. 数值验证（PASS/FAIL）

| 编号 | 检查 | 结果 |
|---|---|---|
| B1 | 零驱动准能量 vs 静态能带（mod ħω） | dev **6.5e-16 eV** |
| B1 | 零驱动周期平均自旋 vs 静态自旋 | dev **1.2e-13** |
| B2 | 光子权重归一化（m_sp=12） | max\|Σw−1\| = **3.8e-14** |
| B3 | Sambe 扩展哈密顿量 vs 时间演化准能量（n_harm=3） | **8.2e-12 eV** |
| B4 | 高频 van Vleck 收敛（ω=1,2,4,8 eV） | 误差 1.27e-8 → 3.16e-9 → 7.91e-10 → **1.98e-10 eV**；ratio(8/1)=**1.56e-2** |
| B5 | 自旋谱权重恒等式 | **6.1e-16** |
| B5 | 自旋谱函数能量积分 | **6.2e-16** |
| B6 | micromotion 检验（S̄ vs t=0 快照） | **8.9e-3**（报告差异，不忽略） |

跨实现回归：在同一份 PbNV 数据上，新引擎与已验证 PbNV 引擎的准能量差
**<1e-8 eV**（机器精度）；数据集更新后用 `gen_old_engine_ref.py`
重新生成参考值（不放宽容差）。

## 3. 自旋分辨谱与别名控制

- `A(k,E)` 使用未跟踪谱重构（不做能带标签跟踪），只依赖 `w_{αm}` 与 `E=ε−mω` 样本；
- bin 宽度显式写入输出，谱函数按 bin 宽度归一；
- 深能带别名检查：γ=0 时全路径 0 个窗口内别名；γ=0.2 时 21 个 k 出现 84 个
  混叠实例，但最小 `|E_alias|=4.83 eV`，远在窗口边界 2.21 eV 之外 → 窗口干净；
- 采样 Nyquist：`|m|<n_step/2=75`（|E−ε|<30 eV），默认 m_sp=16 足以覆盖谱窗。

## 4. PbNV 应用（Stage B）

数据：新 run（`wannier90.chk/eig/hr/spn`，validate 24 PASS / 4 WARN，0 FAIL）。

- 静态自旋纹理（带 111 = VBM）：`Sz ≈ +1 (ky>0)`、`Sz ≈ −1 (ky<0)`，
  在 ky≈0 附近经过面内 (1,1) 方向过渡带；自旋几乎饱和（|S|≈0.96）。
- Floquet 自旋纹理（ħω=0.4 eV, γ=0.1, 25 个 k）：
  - Floquet 态按 overlap 对到静态带，匹配纯度中位 **0.995**，90.6% 的点 >0.8；
  - 光致自旋改变高度局域在 avoided crossing / Y–Γ 附近：
    带 110–111（VBM 对）max `|ΔS|` = **0.351**；
    带 112 max **1.031**，带 113 max **1.221**；
    中位变化 1e-4–0.08；
  - 不做单带标签比较（直接按准能量排序会得到虚假的大变化，已修正）。

## 5. 注意事项（写进报告，不隐藏）

1. 新 `.spn` 为 VASP LWRITE_SPN 数据；修正 TRS 判据后残差仅 2e-3。
   Stage B/C 的自旋结论仍有数据近似（~1e-3），不给出“精确 SOC 常数”。
2. Floquet 网格 25 个 k 为演示规模；出图/定量拟合时可增大到 81 个 k（约 5–6 min）。
3. 简并/近简并 Floquet 态必须按子空间/overlap 分析；当前实现给出匹配纯度，
   纯度低的点应在后续分析中剔除或做简并子空间投影。

## 6. 复现

```bash
python floquet.py validate --config configs/pbnv.json
python floquet.py benchmark --which floquet_rashba
python floquet.py spectrum --config configs/pbnv.json --nk 5
python floquet.py spin_floquet --config configs/pbnv.json --kmesh-n 5 --kmax 0.12
python -m unittest tests.test_stage_b -v
```
