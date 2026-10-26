# 阶段 B 验收报告（Floquet spin texture / solver）

对象：解析 Rashba benchmark + PbNV / Pt / MoSSe 三套 Wannier 数据（只含 SOC-on Wannier）。
核心：v0.1.7-rc1 工作区。

## 1. 解析 benchmark（独立正确性证据）

`python floquet.py benchmark --which floquet_rashba`：

| 检查 | 结果 |
|---|---|
| γ=0 准能级恢复静态 | 6.5e−16 eV |
| γ=0 自旋纹理恢复 | 1.2e−13 |
| 光子权重归一化 | 3.8e−14 |
| Sambe vs 时间演化 | 8.2e−12 eV |
| 高频有效收敛 | 1.98e−10 eV（ω=8 eV） |
| 自旋谱权重恒等式 | 6.1e−16 |
| **总判定** | **PASS 8 / FAIL 0** |

## 2. 三种材料

命令：`floquet_points --points ... --gamma 0 0.02 --pol circ+`；MoSSe 用 K 点（n_step=150, m_sp=16），PbNV/Pt 用较小设置（40–80 步，m_sp=8）。

| 材料 | 体系 | validate | Floquet solver | 带边自旋 |
|---|---|---|---|---|
| MoSSe | 2D，34 Wannier，gapped | 28 PASS / 0 FAIL / 1 WARN | unitarity 3.6e−14；γ=0 `dS=0`；γ=0.02 max\|dS\|=1.54e−4；w_main≈0.9996 | **PASS**（4 带全部 reliable） |
| PbNV | 2D，156 Wannier，gapped | 27 PASS / 0 FAIL / 2 WARN | unitarity 1.2e−14；w_main≈0.9997 | **WARN**：Γ/X 的目标带间 gap < 5 meV，单带匹配全部 unreliable |
| Pt | 3D fcc，18 Wannier，金属 | 25 PASS / 0 FAIL / 4 WARN | unitarity 1.3–1.9e−14 | **WARN**：金属无带边；目标态简并，单带光子权重无意义 |

## 3. 判定与解释

- **阶段 B 的核心 solver/自旋平均逻辑：PASS**（解析 Rashba 独立验证 + 三材料 unitarity）。
- **MoSSe：PASS**（这是目前唯一同时满足“gapped + 带边劈裂足够 + 自旋算符可靠”的材料）。
- **PbNV：WARN**。PbNV 的 SOC 劈裂很小（<5 meV），当前 `floquet_match_tol=5 meV` 把所有目标带标为不可靠。这不是程序错误，而是**单带标签不适用**；需要用简并子空间投影/自旋谱权重分析。
- **Pt：WARN**。金属体系没有“带边”；目标带简并且光子权重对简并基矢敏感。只能用于 solver/unitarity 测试，不能用于 SOC 物理结论。

> 因此：**不能**用 PbNV/Pt 的这组结果声称“Floquet 调控了 SOC”；只能写“solver 在真实多带数据上可运行且 unitarity 正常”。

## 4. 只提供 SOC-on Wannier 是否够？

**对 Floquet 调控研究：够。** 你可以用 SOC-only 数据做能带/谱函数、自旋纹理、PST、Rashba/Dresselhaus 有效趋势、频率/强度/偏振扫描。

**对以下项目：不够。**

- SOC-on/off 分解 `H = H_scalar + λH_SOC`；
- “纯 SOC 强度变化”的定量声明；
- 与 SOC-off 的差分对照。

没有匹配的 SOC-off 数据时，结果只能标为 **spin-related spectral splitting**，不能标为 SOC strength；这是规格里明确要求的语义区分。若以后需要 SOC 分解，需要同一 Wannier 子空间/规范下的 SOC-off 计算，或解析/对称性分解，不能直接把两组不同计算的矩阵相减。

## 5. 复现

```powershell
Set-Location D:\Floquet\floquet_soc
python floquet.py benchmark --which floquet_rashba
python floquet.py validate --config configs\pbnv.json
python floquet.py validate --config configs\pt.json
python floquet.py validate --config D:\Floquet\materials\MoSSe\configs\mosse.json
python floquet.py floquet_points --config D:\Floquet\materials\MoSSe\configs\mosse.json --points K --gamma 0 0.02 --pol circ+ --n-step 150 --m-sp 16
```
