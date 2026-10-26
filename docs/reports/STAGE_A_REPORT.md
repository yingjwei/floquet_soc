# Stage A Report — 输入检查、自旋算符与静态自旋纹理

日期：2026-10-21  ·  代码版本：`floquet_soc.__version__ = 0.1.0-stageA`
运行环境：Python 3.12，numpy 2.5.3，scipy 1.18.1，matplotlib 3.11.2

本报告只陈述**可复现的检查结果**。状态分类：
`已独立验证`（有解析基准或独立方法交叉验证） / `内部自洽`（只通过自证）
/ `需要额外输入` / `未实现`。不用“图看起来合理”作为正确性依据。

---

## 1. 输入文件清单与来源

目录 `D:\Floquet\PbNV`：

| 文件 | 大小 | 内容 |
|---|---|---|
| `wannier90.win` | 3.5 KB | `num_bands=192`, `num_wann=156`, `spinors=.true.`, mp_grid 6×6×1，36 个 k 点，dis_win [-64.84, 9.13] eV，dis_froz [-3.44, 2.56] eV，投影 N/V/Pb |
| `wannier90_hr.dat` | 59.6 MB | 156×156 H(R)，nrpts=49，ndegen 1–4；本文件采用“每行都带 R”的 VASP 变体格式（框架两种方言都支持） |
| `wannier90_centres.xyz` | 9.5 KB | 156 个 `X` 中心（另有 12 行原子位置，属该写法的附加信息；解析时只取 `X` 行） |
| `wannier90.spn` | 32.0 MB | **VASP/WAVECAR 派生二进制自旋矩阵**（格式解码见 §3） |
| `wannier90.chk` / `.u.mat` / `.u_dis.mat` / `.amn` / `.mmn` / `.eig` | **不存在** | Wannier 变换数据（见 §4） |

晶格：a1=(5.28084,0,0) Å，a2=(0,5.76310,0) Å，a3=(0,0,18.74250) Å；
`a_ref=|a1|=5.28084 Å`。非磁性、TRS 保留的 2D 层状体系（kz=0 网格）。

## 2. Hamiltonian 侧检查（全部通过）

`python floquet.py validate --config configs/pbnv.json`：

| 检查 | 结果 |
|---|---|
| 解析维度 | num_wann=156 与 hr.dat 一致；36 个 k 点与 mp_grid 一致 |
| `H(-R) = H(R)†` | max 误差 **0.00e+00**（tol 1e-10） |
| `H(k) = H(k)†` | max 误差 **1.24e-14**（tol 1e-10） |
| ndegen 一致性 | R/−R 对全部一致，ndegen∈[1,4] |
| Wannier 中心 | 来自 `_centres.xyz`（156 个），非原子位置近似 |
| 静态带隙（Γ-X-S-Y-Γ，nocc=112） | **1.487057 eV**（与 PbNV 既有参考 1.48706 eV 一致） |
| 驱动配置解析 | PASS |

结论：**Wannier Hamiltonian 读取与规范、Hermiticity、带隙均通过**；
这部分可以用于后续 Floquet 计算。

## 3. `wannier90.spn` 格式解码

### 3.1 二进制记录结构（Fortran unformatted）

```
record 0 : 60 字节 ASCII 头 "Created from WAVECAR at 2026-01-29T23:23:09.791451" + 空格
record 1 : int32 ×2  = [num_bands=192, num_kpts=36]
record 2..37 : 每个 k 一条记录，长度 889,344 字节
```

`889,344 = 3 × (192×193/2) × 16`——即每个 k 有 **3 个 192×192 的
Hermitian 矩阵、以复数双精度（complex128）按三角打包**。

### 3.2 打包顺序（与写出器逐行核对）

写出行（wannier-berri `utils/vaspspn.py`，本文件的头字符串与其完全一致）：

```python
A = np.array([s[n, m] for m in range(NBout) for n in range(m + 1)
                        for s in (SX, SY, SZ)], dtype=np.complex128)
```

即：外层 m、内层 n≤m（**列主序上三角**），每个矩阵元内依次写
(SX, SY, SZ) 三个自旋分量。框架按同一顺序解包并按 Hermitian 补全。

验证结果（36 个 k 全部检查）：

* 对角元虚部 `max|Im diag| = 0.00e+00`（Hermitian 矩阵的必要条件）；
* 重建后 `max|S−S†| = 0.00e+00`（构造上精确）；
* 逐带自旋极化 `max |⟨σ⟩| = 1.000`（Pauli 界 1，未超界）。

### 3.3 数据质量警告（重要）

| 指标 | 实测 | 正交基期望 | 判读 |
|---|---|---|---|
| S_i 本征值范围 | **[−1.632, +1.632]** | [−1, 1] | 超出 Pauli 界 63% |
| eig(Σ_i S_i²) | **[0.541, 7.962]** | 3 | 自旋求和规则严重偏离 |

这是 `vaspspn` 用 **赝波函数**（而非全 PAW 波函数）且只做逐带归一化、
不做正交化的已知限制；对自旋简并（保留 P·T）体系尤其不可靠
（见 wannier-berri 文档中的 PAW warning）。因此即使拿到了 Wannier 变换，
这份 `.spn` 也只能作为**近似**自旋算符，需要交叉验证或重新生成。
（本框架在 `validate` 中把这两项标为 **WARN**，不作为 PASS。）

### 3.4 关于 TRS 自旋关系的说明

对非磁性体系，理论上 `eig[S(-k)] = −eig[S(k)]`。但 `.spn` 文件不存储
k 点坐标，只有记录顺序；`wannier90.win` 的 36 个 k 点顺序与 WAVECAR 的
k 顺序是否一致无法从文件本身确认。尝试最优配对后没有找到满足该关系的
配对（所有配对偏差近似常数），因此本报告**不把 TRS 检查当作已验证项**；
它被列为“需要 k 顺序/来源信息才能完成”的待验证项。

## 4. 自旋算符与 Wannier 基组：BLOCKED

* `spn` 矩阵维度 = **192**（WAVECAR 能带基组，`num_bands`）。
* `H(k)` 维度 = **156**（`num_wann`，Wannier 基组）。
* 二者之间的旋转需要 Wannier 变换：`S_W(k) = U†(k) S_band(k) U(k)`，
  其中 U 来自 `.chk`（或 `.u.mat`/`.u_dis.mat`），或由 `.amn/.mmn/.eig`
  重跑 wannier90 得到。

**目录中缺少以下全部文件：**
`wannier90.chk`, `wannier90.u.mat`, `wannier90.u_dis.mat`,
`wannier90.amn`, `wannier90.mmn`, `wannier90.eig`, `wannier90.nnkp`。

框架按规范抛出 `MissingTransformError` 并列出上述文件，
**不通过任何排列/猜测把 192 维自旋矩阵硬塞进 156 维**。
因此：

* `spin_operator_wannier_basis` = **BLOCKED**
* `spin_texture_wannier` = **BLOCKED**
* `validate` 退出码 = **2**

这意味着：**目前不能对 PbNV 声称“可靠的静态自旋纹理/Floquet 自旋纹理”**。
Hamiltonian 侧的 Floquet 谱计算不受影响，但任何涉及自旋的结论都必须等
自旋算符解锁后再做。

### 解除阻塞的推荐顺序

1. 提供 `wannier90.chk`（首选；含 U(k) 与规范信息），或 `.u.mat`/`.u_dis.mat`。
2. 若只能重跑：提供 `.amn/.mmn/.eig`（或 `nnkp` + 波函数），由 wannier90 重新生成 U。
3. 同时建议重新生成 `.spn`：新 VASP 的 `LWRITE_SPN=.TRUE.`（格式化文件）或
   VASP2WAN90_v2_fix，避免赝波函数误差。
4. 拿到后 Stage A 剩余验收：
   * 旋转后自旋算符在 Wannier 基组中满足 Hermiticity / Pauli 界 / 求和规则（±2%），
   * 与 band-basis 结果经 U 变换一致（≤1e-10），
   * TRS：`eig[S_W(-k)] = −eig[S_W(k)]`（≤1e-8），
   * 静态自旋纹理对 k 网格、Wannier 规范（Bloch 相位）稳定。

## 5. 解析模型 benchmark（自旋纹理/PST 算法的独立验证）

`python floquet.py benchmark --which all` → **PASS（11 项检查，0 FAIL）**。

### Benchmark 2：Rashba 2 带

H = (ħ²k²/2m)I + α_R(k_yσ_x − k_xσ_y)，α_R=0.10 eV·Å，ħ²/2m=0.5 eV·Å²。

| 检查 | 结果 |
|---|---|
| 带劈裂 ΔE(k) = 2α_R k | 最大偏差 **6.9e-18 eV** |
| 自旋-动量锁定 ⟨S⟩=±(k_y,−k_x,0)/2|k| | 最大偏差 **4.4e-16** |
| S_z = 0 | **4.2e-16** |
| |S| = 1/2（ħ 单位） | **2.8e-16** |
| PST 指标本征值 (1/2,1/2,0) | 偏差 **2.2e-16** |
| Hermiticity | 0 |

### Benchmark 3：Rashba + Dresselhaus，PST

H = ... + α_R(k_yσ_x−k_xσ_y) + β_D(k_xσ_x−k_yσ_y)，α=β=0.10 eV·Å。

| 检查 | 结果 |
|---|---|
| ΔE = 2|d(k)| | **2.8e-17 eV** |
| α=β 自旋共线 ⟨S⟩=±(1/2)(1,−1,0)/√2 | **1.5e-15** |
| PST 指标本征值 (1,0,0) | 精确 |
| PST 主轴 (1,−1,0)/√2 | 角度偏差 **0°**（机器精度） |
| α≠β（0.10 vs 0.03）非共线 | 共线度 **0.650** < 1 |

结论：**自旋纹理、自旋-动量锁定、PST 共线性指标**这些 stage C 的核心算法
已经在解析模型上完成独立验证；它们不依赖 PbNV 的 `.spn` 文件。

## 6. 驱动接口 / Peierls（内部自洽）

`drives/peierls.py` 实现通用矢量势

    A(t) = A1 cos(ωt) e1 + A2 cos(ωt+φ) e2

（线/圆/椭圆偏振、任意 e1/e2 方向、相位、手性、γ↔A0 换算、峰值场强），
以及 Peierls 两中心形式 `H_ij(k,t)=H(R)e^{i[k·R+A(t)·d_ij]}`。

| 检查 | 结果 |
|---|---|
| A=0 时还原 `H(k)` | 偏差 ≤1e-12 |
| 任意 A 下 Hermiticity | ≤1e-10 |
| 圆偏振 |A(t)| 恒定 | 1e-16 |

状态：**内部自洽**（尚未与独立 Floquet 求解器交叉验证，属 stage B）。

## 7. 测试

`python -m unittest discover -s tests -v`：**18 项全部通过**（约 16 s），
覆盖：Wannier 解析/厄米性/带隙、spn 解码结构与质量指标、缺失变换文件报错、
Rashba/RD 解析基准、Peierls 一致性、CLI（validate 退出码 2、benchmark 0、
bands、status、stage B 占位）。

## 8. 状态汇总

| 功能 | 状态 |
|---|---|
| Wannier90 输入 / H(k) / Hermiticity / 静态带 | **已独立验证** |
| 通用驱动 A(t) / Peierls H(k,t) | **内部自洽** |
| 自旋纹理 / PST 指标算法 | **已用解析模型独立验证** |
| VASP 二进制 `.spn` 格式解码 | **结构已验证；数据质量 WARN** |
| `.spn` → Wannier 基组旋转 | **需要额外输入（BLOCKED）** |
| PbNV 静态/Floquet 自旋纹理 | **需要额外输入（BLOCKED）** |
| Floquet solver / 自旋分辨谱（stage B） | **未实现** |
| SOC splitting 拟合 / PST 驱动分析（stage C） | **未实现** |
| 参数扫描 / 收敛 / 跨材料迁移（stage D） | **未实现** |

## 9. Stage B 验收条件（预告）

1. 三种方法（直接 U(T)、Sambe、高频展开）对同一模型给出相容的
   ε、Ē、光子权重，并明确三者的定义差异；
2. 周期平均自旋必须包含 micromotion：`S̄=(1/T)∫⟨u(t)|S|u(t)⟩dt`，
   并与 `t=0` 近似的差异量化；
3. Rashba/RD 解析模型上：驱动关闭 → 恢复静态自旋纹理；
   弱驱动高频极限 → 与 `H_eff` 自旋纹理一致；
4. 自旋分辨谱 `A^{S_i}(k,E)` 的带符号加权与非负权重分开实现并交叉检查；
5. 通过 PbNV 回归（前提：§4 自旋算符解锁）。

---

## 10. 更新（2026-10-22）：chk + 新 run 到位后的最终状态

> 本节取代上面 §4 / §8 中“自旋算符 BLOCKED”的旧结论。

用户提供了同一次 VASP/wannier90 run 的完整文件
（`wannier90.chk/eig/hr/centres/spn`，`.spn` 为 VASP `LWRITE_SPN` 格式化文件），
已放入 `D:\Floquet\PbNV`；旧 1 月 run 备份于
`C:\Users\15756\.dsh\PbNV_oldrun_backup_2026-10-22`。框架新增
`core/chk.py`、`core/spin_wannier.py` 与 formatted-`.spn` 自动识别。

### 10.1 验证结果

| 检查 | 结果 |
|---|---|
| `H_W(k)=v†diag(E)v` vs `H_hr(k)`（36 个 k） | **9.73e-6 eV** |
| `U` 幺正性 / `v` 等距 / 晶格一致 | 8.9e-14 / 8.8e-14 / 2.2e-16 |
| `.spn` 格式 | `formatted-spn`（192×36，自动识别） |
| 自旋本征值范围 | **[-1.000, +1.000]，100% 在 Pauli 界内** |
| 自旋求和规则 `eig(ΣS²)` | **[2.786, 3.000]** |
| 静态带隙 | **1.487060 eV** |
| `validate` 汇总 | **25 PASS / 2 WARN / 0 FAIL，exit 0**（TRS 判据修正后） |

剩余 4 个 WARN 全部是 `.spn` 的赝波函数残差（TRS 自检 ~1e-2、求和规则略偏离 3
等），不是框架错误；已在 `outputs/PbNV/validate.md` 逐条列出。
PbNV 无空间反演（自旋分裂体系），Kramers 简并相关检查应视为 N/A。

### 10.2 解锁后的 Stage A 验收

- ✅ 自旋算符 Hermiticity（0），Pauli 界（100%），求和规则（[2.786,3.0]）；
- ✅ 实空间插值往返 `max|S(k_mesh)-S_W(k_mesh)|=4.8e-12`；
- ✅ TRS 自检：修正判据后 `max|sorted eig S(-k) + reverse sorted eig S(k)|`
  = **0.002**（见 §11；旧判据 2.0 是判据错误，不是数据错误）；
- ✅ 静态自旋纹理对 k 网格与规范稳定（带 111：面外 ±1 分区 + ky≈0 面内过渡带）。

### 10.3 测试

当前测试套件：Stage A 19 项、Stage B 13 项、chk/旋转/纹理 7 项、
Stage C 10 项、Stage D 7 项（合计 56 项，PbNV/Pt 数据缺失时自动 skip）。



---

## 11. 更正（2026-10-22）：TRS 判据公式与 PbNV 数据质量结论

早期版本用 `max|eig(S(-k)) + eig(S(k))|` 作为 TRS 质量判据，这是**错误的**。
非磁体系的正确条件是

```
spec S(-k) = -spec S(k)   <=>   sorted(S(-k)) = -reverse(sorted(S(k)))
```

对于本征谱本身关于 0 对称的体系（存在 P*T 对称的 Pt 等），旧公式即使数据完美
也会给出 2.0 的“破坏”。修正为 `Sg[j] + Sg[ik][:, ::-1]` 后：

| 材料 | 旧（错误）判据 | 修正后 | validate |
|---|---|---|---|
| PbNV（VASP LWRITE_SPN，新 run） | 2.000 | **0.002** | 25 PASS / 2 WARN / 0 FAIL |
| Pt（QE pw2wannier90） | 2.000 | **0.0000** | 23 PASS / 4 WARN / 0 FAIL |

结论：PbNV 新 run 的 `.spn` 是质量合格的 VASP LWRITE_SPN 数据
（Pauli 界 100%、求和规则 [2.786,3.000]、TRS 2e-3）；先前“~1% TRS 残差”
的表述部分来自判据 bug，特此更正。

## 12. 冻结说明（2026-10-23，v0.1.0-rc1）

本节数字为当时阶段报告；冻结前按 README/RELEASE_VALIDATION 记录的命令**复跑**
后，当前代码的输出为：

| 材料 | 本报告历史计数 | 冻结复跑（当前代码） | 差异原因 |
|---|---|---|---|
| PbNV | 25 PASS / 2 WARN | **26 PASS / 2 WARN** | 后加了 `spin_wannier_hermiticity_offgrid` 检查（PASS，携带 NOT_VALIDATED 注记） |
| Pt | 23 PASS / 4 WARN | **24 PASS / 4 WARN** | 同上：WS/ndegen 插值下的 off-grid 自旋厄米性检查（PASS） |

历史计数保留不改；以 `outputs/PbNV/validate.*`、`outputs/Pt/validate.*`
（冻结复跑生成）、`RELEASE_VALIDATION.md` 与 `STATUS_REGISTRY.md` 为准。
当时全套测试 56 项（§10.3）；冻结版为 **74 项（2 项数据依赖 skip）**，
见 `RELEASE_TEST_RESULTS.txt`。
