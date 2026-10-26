# RELEASE NOTES — v0.1.2-rc1（figure style update）

日期：2026-10-23　｜　基线：v0.1.1-rc1（commit `e819463`，cleanup release）
上一基线 v0.1.0-rc1（commit `a9b8e9a`）均保留且不可改写。

## 1. 这个版本改了什么

**只改图形呈现 + 重新生成对应图片；物理代码、数值结果、状态登记均未变。**

1. **新增 `floquet_soc/plot_style.py`**：参考
   `C:\Users\15756\Desktop\phonkit\src\phonkit\band\plot.py` 的 PRL 风格
   （Arial/stixsans、轴线宽 1.2、字号 14/11/12、Wong 2011 色盲友好配色、
   高对称点竖虚线 + Γ 标签、600 dpi、单栏 8 cm 方形图）。
2. **`cli._plot_bands` 改为 PRL 风格**：k 轴由 k index 改为累积笛卡尔 k 距离（1/Å）；
   半导体能量零点取 VBM（`num_electrons` 存在时），y 轴标注 `Energy − E_VBM`；
   输出 `outputs/<material>/bands.png`（600 dpi）。
3. **`cli._plot_soc` 改为三面板 PRL 风格并修正尺度混用**：
   (a) SOC 量级的相邻带分裂（`max ΔE ≤ 0.5 eV`，meV）；
   (b) 最大相邻对 = 直接带隙（eV，单独面板，标注 min gap）；
   (c) SOC 窗口内能带（带边对高亮）。
   **带隙不再与 SOC 分裂画在同一 ΔE 轴上**（此前 111-112 的 ~2.3 eV 会把
   ≤0.1 eV 的 SOC 曲线压成 0，是已知读图陷阱）。
4. **`cli._plot_spin_pairs` 改为 PRL 风格**：k 距离横轴、配对曲线用 Wong 扩展配色、
   图例带框；仅画 spin-paired 配对。
5. **重新生成并覆盖以下图片**（命令与冻结版相同，数据文件逐字节不变）：
   - `outputs/PbNV/bands.png`（`floquet.py bands --config configs/pbnv.json`）
   - `outputs/PbNV/soc_splitting.png`、`outputs/PbNV/soc_spin_splitting.png`
     （`floquet.py soc --config configs/pbnv.json`）
   - PbNV 的 `bands.npz`、`soc_splitting.npz` 与 Pt 的 `soc_splitting.npz`
     **git diff 为空（byte-identical）**，物理数据未变。

## 2. 什么没有变

- `floquet_soc/{core,drives,solvers,observables,benchmarks,workflows}/` **零改动**；
- 全部 npz 数据、`STATUS_REGISTRY.md`、WARN / NOT_APPLICABLE / NOT_VALIDATED /
  UNCONVERGED / BLOCKED / NEEDS_INPUT 状态**原样保留**；
- 自旋纹理图（bwr + 箭头）、`scan/converge/spectrum` 图、benchmark 图本轮**未改**，
  保持原风格（后续如需统一可再发版）；
- `outputs/_tests*` 不进入发布（v0.1.1-rc1 起的 `.gitignore` 规则）。

## 3. 测试

在 v0.1.2-rc1 树的干净副本上运行：
`FLOQUET_SOC_JOBS=1 python -m unittest discover -s tests -v`
→ 结果记录在 `RELEASE/_v0.1.2_test_run.txt`（预期与 v0.1.1-rc1 相同：
Ran 74 tests, OK (skipped=2)）。

## 4. 与历史版本的关系

| 版本 | 内容 | tag |
|---|---|---|
| v0.1.0-rc1 | 初始冻结基线（物理验证） | `v0.1.0-rc1`（保留） |
| v0.1.1-rc1 | 清理测试产物（cleanup-only） | `v0.1.1-rc1`（保留） |
| v0.1.2-rc1 | 本版：能带/SOC 图改为 phonkit PRL 风格 | `v0.1.2-rc1` |

引用规则：物理结论三版相同；**正式图引用 v0.1.2-rc1 快照**；
历史复核用对应 tag。
