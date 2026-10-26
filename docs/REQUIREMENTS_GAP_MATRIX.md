# 需求–代码差距矩阵（floquet_soc v0.1.7-rc1 工作区）

> 判定口径：**PASS** = 有独立 benchmark/数值容差；**WARN** = 有实现与内部自洽检查，但缺独立验证或输入；**GAP** = 未实现；**N/A** = 明确不在当前阶段。
> 证据以仓库代码/测试/`floquet.py` 命令为准，不以“图看起来合理”为依据。

| 规格条目 | 当前状态 | 证据 | 缺口 / 下一步 |
|---|---|---|---|
| 2 输入解析：win/hr/centres | **PASS** | `core/wannier.py`；`tests/test_stage_a.py`；`validate` | — |
| 2 输入解析：spn 基组识别 | **PASS** | `core/spin_operator.py::read_wannier_spn` 读 VASP 记录；`num_bands != num_wann` → `vasp-bands` | `.u.mat/.u_dis.mat` 路径未支持（只支持 `.chk`+`.eig`） |
| 2 自旋旋转 `S_W=v†S v` | **PASS（三材料验证）** | PbNV 156、Pt 18、MoSSe 34；H_W−H_hr ≤2.3e−5 eV；离网格厄米 ≤1.1e−15 | `.u.mat` 已显式拒绝（提示用 chk+eig）；磁性子空间未测 |
| 2 自旋矩阵检查 | **PASS** | `validate_spin_matrices`：厄米、实对角、特征值界、ΣS² | 自旋算符来源/哈希已由 `provenance` 记录 |
| 2 与 DFT 自旋交叉验证 | **WARN** | `MOSSE_DFT_SPIN_COMPARE.md`：方向 RMS 0.0082、Sz corr 0.9991；能量偏移 −3.5 meV；DFT \|S\|≤0.78 | 需确认 DFT 输出归一化；K′ 不在 DFT 网格 |
| 3 方法 A 直接演化 | **PASS** | `solvers/evolution.py`；`tests/test_stage_b.py` | — |
| 3 方法 B Sambe | **PASS** | `solvers/sambe.py`；与 A 对照测试 | 高维截断自动选择可加强 |
| 3 方法 C 高频有效 | **PASS（有适用条件）** | `solvers/effective.py`；收敛 benchmark | 低频/强共振自动拒绝策略需显式化 |
| 3 ε_F / Ē_F / E_eff 区分 | **PASS（代码/文档）** | `spectral.py`、`floquet_spin.py` 返回分字段 | 报告层需统一避免混称 |
| 4 驱动接口 | **PASS（Peierls）** | `drives/peierls.py`：ω、E0/γ、linear/circular/elliptical、φ、e1/e2、sign | `drives/custom.py`、onsite/hopping 调制、周期 Zeeman、H_drive(t) 未实现 |
| 4 参数批量扫描 | **PASS（部分）** | `scan`、`converge`、`floquet_points` | 扫描结果统一 schema/图集需完善 |
| 5.1 SOC-related splitting | **WARN** | `soc`、`soc_splitting.py`、Rashba benchmark | 自动区分 SOC/轨道杂化/光子复制未实现；报告只能用“spin-related spectral splitting” |
| 5.2 静态/Floquet 自旋纹理 | **PASS** | `spin_texture`、`spin_floquet`；Rashba benchmark；micromotion 平均 | 多带近简并子空间表示已有 `reliable` 标记，需在所有输出统一 |
| 5.2 spin-resolved A^{S_i} | **PASS（实现）** | `spin_spectral_density`；带符号/非负分离 | 缺独立解析 benchmark（建议加 2-band 磁模型） |
| 5.3 PST 指标/旋转 | **PASS** | `pst.py`、Rashba/RD benchmark；主轴本征矢 | 缺真实材料 PST 对照（无可靠 DFT 自旋纹理覆盖） |
| 5.4 有效 SOC 项/对称性约束 | **WARN/GAP** | 高频 `H_eff` 已实现；Rashba 拟合已实现 | 对称性允许项枚举、Dresselhaus/高阶拟合模板、磁交换分离未实现 |
| 5.5 SOC-on/off 与 λ 分解 | **GAP** | 无接口 | 需在基组/规范匹配后实现；否则拒绝相减 |
| 6 Benchmark 1（无 SOC 自旋简并） | **PASS** | `benchmarks/spin_degenerate.py` 3/0 | — |
| 6 Benchmark 2（Rashba） | **PASS** | `benchmarks/rashba.py`；ΔE=2αk、自旋锁定 | — |
| 6 Benchmark 3（Rashba+RD） | **PASS** | `benchmarks/rashba_dresselhaus.py`；PST/PSH | — |
| 6 Benchmark 4（Kane–Mele） | **GAP** | 无 | 新增内禀 SOC/光场能隙测试 |
| 6 Benchmark 5（磁性 SOC） | **PASS** | `benchmarks/magnetic_soc.py` 4/0（Γ 交换 2M、有限 k SOC） | 多带磁性 Wannier 仍缺 |
| 6 Benchmark 6（PbNV 真实 Wannier） | **PASS/WARN** | `benchmarks/pbnv_spin.py`；PbNV 回归 | 156 轨道全量开销大，CI 用小 kmesh |
| 6 Benchmark 7（跨材料） | **PASS（Pt）** | `benchmarks/pt.py`；18 Wannier、chk+eig | 还需一个不同轨道数/点群材料做“无需改核心”验证 |
| 7 数值检查 1–15 | **PASS（多数）** | `validate`、stage A–D tests、benchmark reports | 光子扇区归一化/折叠规范不变性需单列回归；近简并统一用子空间 |
| 8 架构目录 | **PASS（近似）** | `core/drives/solvers/observables/workflows/benchmarks/tests` | `drives/custom.py` 缺失；`observables/soc_fitting.py` 未独立 |
| 8 CLI 命令 | **PASS** | `validate bands spectrum spin_static spin_floquet soc pst scan converge benchmark`（+ `spin_texture dft-spin-check floquet_points provenance`） | 命令统一 `--config` schema 文档化 |
| 8 元数据/可复现 | **PASS（新增）** | `floquet.py provenance`：版本、config/输入 SHA256、numerics、spin source | 需把 provenance 注入所有 JSON 输出（目前已注入 `floquet_points`、`dft-spin-check`） |
| 9 阶段 A 验收 | **PASS（带 WARN）** | 见 `PHASE_A_ACCEPTANCE.md` | DFT 幅度归一化、`.u.mat` 兼容 |
| 9 阶段 B | **PASS（解析+Rashba）/ WARN（PbNV,Pt）** | `PHASE_B_ACCEPTANCE.md`、`PHASE_B_SUBSPACE_AND_BENCHMARKS.md`；子空间自旋已实现 | PbNV/Pt 子空间量可算，但单带 Ising 标签不可用；磁性 Wannier 缺 |
| 9 阶段 C | **WARN** | SOC splitting/PST/Rashba 拟合 | 对称性约束有效项、SOC-on/off 缺 |
| 9 阶段 D | **WARN** | 多 benchmark、CLI、scan/converge | 新材料教程、统一 PASS/FAIL 总报告、跨材料第二例缺 |

## 下一阶段优先级（按规格顺序）

1. **Phase A 收尾**：`.u.mat` 读取/拒绝策略；provenance 注入全部输出；自旋基组元数据写入 `validate` 报告。
2. **Benchmark 1 + 5**：无 SOC 自旋简并、Rashba ferromagnet；建立“交换 vs SOC”判据。
3. **SOC-on/off 接口**：先做基组/规范匹配检查，不匹配直接拒绝矩阵相减。
4. **有效 SOC 项**：从高频 `H_eff` 出发，按点群枚举允许项；Rashba/Dresselhaus 拟合做成可配置模型，不套万能公式。
5. **跨材料第二例**：选一个不同轨道数/点群 Wannier 数据，验证只改 config 即可。
