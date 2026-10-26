# 目录说明（DIRECTORY）

本文件描述仓库的目录分布与“新文件该放哪里”。**根目录只放仓库级元数据**，所有报告与历史文档都在 `docs/` 下。

## 顶层

| 路径 | 作用 |
|---|---|
| `floquet.py` | 源码树 CLI 入口（等价于 `python -m floquet_soc`；`pip install` 后为 `floquet-soc`） |
| `pyproject.toml` | 打包元数据、依赖、console 入口、pytest / ruff 配置 |
| `LICENSE` | MIT |
| `VERSION` | 版本号（与 `floquet_soc/version.py` 保持同步） |
| `requirements.txt` | 最小运行依赖（等价于 `pip install .`） |
| `README.md` | 快速开始与状态摘要 |
| `CHANGELOG.md` | 版本变更合并摘要（详见 `docs/release/`） |
| `DIRECTORY.md` | 本文件 |
| `configs/` | 材料配置（JSON）：`pbnv` / `pt` / `rashba` / `rashba_dresselhaus` |
| `benchmarks_data/` | 随仓库发布的 Benchmark 输入数据（`Pt/`：来源与授权见其 `README.md`） |
| `floquet_soc/` | Python 包（见下） |
| `tests/` | 全部测试 |
| `docs/` | 全部文档（见下） |
| `outputs/` | 运行产物（`.gitignore` 忽略，不入库） |
| `.github/` | CI（GitHub Actions） |
| `.workbuddy/` | **本地 AI 记忆/笔记**（`.gitignore` 忽略，`MANIFEST.in` 排除，绝不进入仓库或发行包） |

## `floquet_soc/` 包

| 模块 | 职责 |
|---|---|
| `cli.py` | argparse、20 个子命令、绘图（待拆分） |
| `__main__.py` | `python -m floquet_soc` 入口 |
| `version.py` | 版本号单一来源 |
| `plot_style.py` | 统一图样式（字号 / 颜色 / 尺寸 / `save_figure`） |
| `core/` | `wannier`（win/hr/chk/wsvec 解析、`WannierModel`、H(k)、厄米性断言）、`chk`、`spin_operator`（spn 解析、`spin_at_k`）、`spin_wannier`（chk 旋转的 S_W）、`symmetry`（TRS 配对）、`validation`（PASS/WARN 报告）、`cache`（模型缓存）、`provenance`（运行溯源）、`model_factory`（`load_config` / `build_model` / `validate_cfg`） |
| `drives/peierls.py` | 通用 A(t)、Peierls H(k,t)、贝塞尔展开与截断守卫 |
| `solvers/` | `evolution`（U(T) 演化）、`sambe`（Sambe 交叉验证）、`effective`（van Vleck 展开） |
| `observables/` | `spin_texture` / `spin_texture_map`（静态自旋纹理）、`spectral`、`floquet_spin`、`floquet_bands`、`floquet_points`、`soc_splitting`、`soc_scan`、`pst`、`rashba_fit`、`kramers`、`spin_plot`（绘图） |
| `benchmarks/` | `rashba`、`rashba_dresselhaus`、`floquet_rashba`、`magnetic_soc`、`spin_degenerate`、`pbnv_spin`、`pt` |
| `workflows/` | `scan`、`converge`、`parallel`（k 点并行） |
| `validation/dft_spin.py` | DFT 自旋纹理与模型对比 |

## `docs/`

| 路径 | 内容 |
|---|---|
| `docs/release/` | 发布说明与测试结果存档（`RELEASE_NOTES_v0.1.x`、`RELEASE_TEST_RESULTS*`、`RELEASE_VALIDATION.md`） |
| `docs/reports/` | 专题报告：`STAGE_A–D_REPORT.md`、`PT_*_REPORT.md`、`PBNV_PST_REPORT.md`、`PHOTON_TRUNCATION_REPORT.md` |
| `docs/STATUS_REGISTRY.md` | 功能状态登记表 |
| `docs/NEW_MATERIAL_GUIDE.md` | 接新材料 / 新 benchmark 的步骤 |
| `docs/CONFIG_REFERENCE.md` | 配置文件参考（每个键的含义与取值范围） |
| `docs/TROUBLESHOOTING.md` | 错误排查指南（异常 / 警告 / 退出码） |
| `docs/REQUIREMENTS_GAP_MATRIX.md` | 需求覆盖矩阵 |
| `docs/PHASE_*_ACCEPTANCE.md` | 各阶段验收记录 |

## `tests/`

- 阶段回归：`test_stage_a/b/c/d.py`
- 数值/物理：`test_ws_grid.py`、`test_wsvec.py`、`test_pt_physics.py`、`test_spin_wannier.py`、`test_spin_texture_map.py`、`test_target_bands.py`、`test_result_correctness.py`、`test_hermiticity_guard.py`
- 工程：`test_cache.py`、`test_cache_key.py`、`test_commercial_p0.py`、`test_floquet_points_outputs.py`、`test_wsvec_flags.py`、`test_provenance.py`、`test_workflows_new.py`、`test_release_regressions.py`

运行：`pytest -q`（CI 入口），或 `python -m unittest discover -s tests -v`。
**测试产物只写系统临时目录**（`tempfile.mkdtemp`，测试结束自动删除）；少数既有参考数据在
`outputs/_tests*`。**源码树（`tests/`、`floquet_soc/`、根目录）内不生成任何文件。**

## `outputs/`（运行产物，git 忽略）

默认输出根：`$FLOQUET_SOC_OUT`，否则当前工作目录下的 `./runs`。

```
<out>/<material>/static/     # validate、bands、spin_static、spin_texture、soc、pst、rashba…
<out>/<material>/floquet/    # floquet_points、spectrum、spin_floquet、scan、converge
<out>/<material>/benchmark/  # benchmark 报告与图
```

> **计划中的变更（下一批）**：改为“一个计算一个子目录” —— `<out>/<material>/<command>/`，
> 并在每个目录写入 `run.json`（argv、生成模块、版本、config 快照哈希、输入哈希），
> 使每个输出目录自描述、可复现。当前 `static/` 子目录被多个命令共用，属于已知问题。

## 新文件放哪里

| 新增内容 | 位置 |
|---|---|
| 新材料配置 | `configs/<name>.json`（键名见 `core/model_factory.py::_KNOWN_TOP_LEVEL_KEYS`，未知键会被拒绝） |
| 新 benchmark 代码 | `floquet_soc/benchmarks/<name>.py`，并在 `benchmarks/__init__.py` 注册 |
| 新观测量 / 求解器 | `observables/`、`solvers/` |
| 新测试 | `tests/test_<topic>.py` |
| 阶段/专题报告 | `docs/reports/` |
| 发布说明 | `docs/release/`（不要放根目录） |
| Benchmark 输入数据 | `benchmarks_data/<name>/`，并附来源与授权说明 |

## 一次计算的数据流

```
configs/<material>.json
  → core.model_factory.load_config（路径解析 + validate_cfg）
  → build_model → WannierModel(win/hr/chk/eig/wsvec) 或解析模型
  → core.cache 缓存（$FLOQUET_SOC_CACHE，schema 版本进键）
  → solvers / observables 计算
  → outputs/<material>/<subdir>/ 下的 npz + png + report.json/md
```
