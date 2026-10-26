# Changelog

所有版本的重大变更按倒序汇总。完整细节与回归数据见
`docs/release/RELEASE_NOTES_v*.md` 与 `docs/release/RELEASE_TEST_RESULTS*`。

## v0.1.7-rc1（当前工作树）

- 静态自旋纹理图与 DFT 自旋纹理交叉对比升级为通用框架命令
  （`spin_texture`、`dft-spin-check`）。
- 工程硬化：结果级 bug 修复（`floquet_points` overlap 行/列笔误、缓存键
  缺 `spin.source`/`wsvec` 开关、Peierls 贝塞尔截断守卫）；`pyproject.toml` +
  `LICENSE` + `python -m floquet_soc` 入口 + console script；`core/model_factory`
  解耦；缓存 schema 版本 / 原子写 / LRU；并行层改 `-m` 入口并清理临时目录；
  测试去桌面绝对路径；CI（ubuntu + windows 矩阵）；配置校验（顶层 + 子键 +
  数值范围）；目录重组与 `DIRECTORY.md`。
- 详见 `docs/release/RELEASE_NOTES_v0.1.7-rc1.md`。

## v0.1.6-rc1

- Wannier90 `use_ws_distance`（`wannier90_wsvec.dat`）应用于 H(k)、Wannier 基
  自旋算符与 Peierls 驱动；修正 off-grid 插值（MoSSe Rashba α 64.8 → 44.3 meV·Å）。

## v0.1.5-rc1

- 修复奇 `mp_grid` 的 WS 超胞；自旋算符与 H(k) 使用同一套 hr.dat R/ndegen；
  新增 `wsvec` 解析器与回归；修复解析模型缓存键碰撞。

## v0.1.4-rc1

- 连续带曲线（去掉 NaN 掩膜）；bands/soc 图统一字号。

## v0.1.3-rc1

- SOC 分裂口径收紧为 gap-edge 两对（VBM-1/VBM、CBM/CBM+1）；双面板 PRL 图。

## v0.1.2-rc1

- 图样式更新（PRL 风格）。

## v0.1.1-rc1

- cleanup-only：测试产物与发布输出分离。
