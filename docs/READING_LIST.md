# floquet_soc 使用者文献清单

定位：你是**使用者**——用这个代码的输出（准能量带、自旋纹理、有效 Rashba 系数、
dressed hr.dat）做凝聚态文章计算，不是要读代码。所以清单按"读它是为了看懂哪个
输出"来组织，而不是按代码模块。

## 第 0 层：先读自己仓库的（半天）

- `README.md` + `docs/CONFIG_REFERENCE.md` —— 每个输入键什么含义
- `docs/TROUBLESHOOTING.md` —— 报错时往哪查
- `tests/test_commercial_p0.py` / `test_a4_continuity.py` 的 docstring —— 代码里
  锁定的物理约定（准能量折叠区间、H_p(k) Fourier 约定）比任何论文都权威

---

## 第 1 层：Floquet 物理入门（按顺序读，约一周）

1. **Eckardt, Rev. Mod. Phys. 89, 011004 (2017)**
   "Colloquium: Atomic quantum gases in periodically driven optical lattices"
   → 最易读的入门。看懂三件事就够：准能量是模 ħω 的量（第一 Floquet 区），
   高频驱动下 U(T) 的有效哈密顿量，边带结构。对应你所有输出的基本语言。

2. **Rudner & Lindner, Nat. Rev. Phys. 2, 229 (2020)**
   "Band structure engineering and non-equilibrium dynamics in Floquet
   topological insulators"
   → Floquet 能带、微运动（micromotion）、准能带边界的处理方式——正是
   `floquet_bands.py` 画谱函数而不直接画准能量的原因（区边缘排序问题）。

3. **Oka & Kitamura, Annu. Rev. Condens. Matter Phys. 10, 387 (2019)**
   "Floquet engineering of quantum materials"
   → 最贴近"做文章"的一篇：光驱动 2D/拓扑材料、圆偏振光打开能隙、
   **prethermalization（加热问题）**。审稿人问"你的准稳态能维持多久"时，
   引文就在这里。

4. **Bukov, D'Alessio & Polkovnikov, Adv. Phys. 64, 139 (2015)**（最厚，放最后）
   → `solvers/effective.py` 的 van Vleck 高频展开和有效性判据的理论出处。
   只需要精读高频展开那一章。

---

## 第 2 层：代码里每个方法的原始文献（写方法部分时引）

| 代码 | 输出 | 文献 |
|------|------|------|
| `evolution.py` | U(T) 的本征相位 = 准能量 | Shirley, Phys. Rev. **138**, B979 (1965)（Floquet 理论开山） |
| `sambe.py` | Sambe 矩阵对角化交叉验证 | Sambe, Phys. Rev. A **7**, 2203 (1973)（扩展 Hilbert 空间） |
| `peierls.py` | Peierls 替代 + 贝塞尔展开 | Goldman & Dalibard, PRX **4**, 031027 (2014)（Floquet-Bloch 带 / 人工规范场）；Oka & Aoki, PRB **79**, 081406(R) (2009)（圆偏振光把石墨烯变成 Chern 绝缘体——圆偏振驱动的"标准模型"） |
| `effective.py` | H_eff = H0 + Σ[H_n,H_-n]/(nħω) | Bukov 2015（见上）；van Vleck 展开原名检索 "van Vleck contact transformation" |
| `spectral.py` | Floquet 谱函数 A(k,E) | Wang et al., Science **342**, 453 (2013)（tr-ARPES 观测 Floquet 边带——你的谱函数图对应的实验量） |
| `spin_operator.py` / spn 解析 | Wannier 基自旋算符 | Ryoo, Park & Souza, PRB **99**, 235113 (2019)（spn 文件方法学，你 VASP `LWRITE_SPN` 的出处） |
| wannier 插值约定 | hr.dat / wsvec | Marzari et al., RMP **84**, 391 (2012)；Pizzi et al., J. Phys.: Condens. Matter **32**, 165902 (2020)（wannier90 社区版，use_ws_distance 约定） |
| `benchmarks/rashba*.py` | 解析基准 | Rashba, Sov. Phys. Solid State **2**, 1109 (1960)；Dresselhaus, Phys. Rev. **100**, 580 (1955) |

---

## 第 3 层：拿输出写文章时的对标文献（2D + SOC + 光）

- **Basov et al.**, "Light-induced emergent phenomena in 2D materials and
  topological materials"（综述, arXiv:2302.00443）→ 全景图：Floquet 边带、
  能隙开合、光-自旋/谷耦合、非线性响应，一张图定位你的计算在哪个象限。
- **McIver et al., Nat. Phys. 16, 38 (2020)** "Light-induced anomalous Hall
  effect in graphene" → 圆偏振光 Floquet 输运的实验锚点；你的理论计算可以
  直接对标这类实验。
- **"Dynamic control of quantum phases in two-dimensional materials via
  Floquet engineering"**（2026 综述, ORNL）→ 最新进展：TMD 自旋-谷物理、
  磁性 2D 材料的光控拓扑（MnBi₂Te₄ 圆偏振相变等）——找"前人做到哪、
  你的材料新在哪"时的地图。
- **谷/自旋方向**：你的 K/K′ 谷不对称输出（`floquet_points` 的 valley
  asymmetry）对标的物理是 TMD 谷电子学 + 圆二色性选择定则；检索
  "valley-selective circular dichroism TMD"（Xiao et al., PRX 2, 041001
  (2012) 是选择定则的原始文献）。

---

## 输出 → 文献速查（写作时用）

| 你手里的输出 | 物理量 | 必引 |
|---|---|---|
| `scan.npz` 的 Ebar/S_bar vs γ | 光修饰能级与自旋 | Eckardt 2017; Goldman-Dalibard 2014 |
| `spectrum.npz` 的 A(k,E) | Floquet 谱函数 / 边带 | Rudner & Lindner 2020; Wang 2013 |
| `rashba_fit` 的 alpha_eff | 光致 Rashba 系数 | Rashba 1960 + 第 3 层任一篇应用 |
| `dressed/wannier90_hr_p*.dat` | 光修饰紧束缚模型 H_p(R) | Oka & Aoki 2009（高频极限即 Haldane 项）；Oka & Kitamura 2019（prethermalization——审稿必问） |
| 准能量带 + Berry 曲率 | Floquet 拓扑 | Kitagawa et al., PRB **84**, 235108 (2011)（驱动系统的拓扑不变量定义） |

## 阅读策略建议

- 第 1 层的 1、2 篇精读，3、4 篇通读——两周内可以建立完整的物理图像；
- 第 2 层**不用通读**，写文章方法节时按表引，只读与你用到的命令对应的段落；
- 第 3 层按你材料体系挑 2-3 篇对标结构（图怎么排、物理量怎么定义）；
- 审稿人最常打的两张牌：**加热/prethermalization**（Oka & Kitamura 2019 回答）
  和**与静态 DFT 结果的一致性**（你仓库里的 dft-spin-check 和 STAGE 报告就是
  为这个准备的）。
