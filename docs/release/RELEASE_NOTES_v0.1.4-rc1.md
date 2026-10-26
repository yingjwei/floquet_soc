# RELEASE NOTES — v0.1.4-rc1（纯绘图修复）

日期：2026-10-23　｜　基线：v0.1.3-rc1（commit `f462abd`，gap-edge SOC splitting）
历史 tag `v0.1.0-rc1` ～ `v0.1.3-rc1` 均保留且不可改写。

## 1. 这个版本改了什么

**只改绘图；数据、SOC 劈裂定义、物理代码全部与 v0.1.3-rc1 相同。**

1. **能带连续绘制（修复"断开"）**：此前 `bands` / `soc` 图为了只显示 SOC 窗口，
   把窗口外能量设成 `NaN`，导致能带在 k 路径中途消失/再出现，看起来像断线。
   现在**整条 k 路径画连续曲线**，仅用 y 轴范围裁剪可见窗口；窗口外的能带自然被
   坐标轴边缘裁掉，不再在路径内部断裂。
2. **字号统一**：此前标题 9–11 pt、刻度 11、轴标签 14、图例 8 混用，视觉上大小不一。
   现在统一为一套：轴标签 **11**、标题/suptitle **10**、刻度 **10**、高对称点 **10**、
   图例 **9**（`floquet_soc/plot_style.py::FONTS`，所有绘图函数引用同一组常量）。
3. **重新生成并覆盖**：
   - `outputs/PbNV/bands.png`（连续能带、统一字号）；
   - `outputs/PbNV/soc_splitting.png`（(a) 连续带边能带；(b) 两对劈裂）；
   - `outputs/PbNV/soc_spin_splitting.png`（两对劈裂，统一字号）。
   **npz 数据文件逐字节不变**（`git diff` 为空），SOC 数值仍是
   110-111 max 29.98 / median 0.51 meV；112-113 max 107.66 / median 0.58 meV。

## 2. 什么没有变

- `floquet_soc/{core,drives,solvers,observables,benchmarks,workflows}/` 零改动；
- SOC 劈裂口径（只报 VBM-1/VBM 与 CBM/CBM+1）与数值不变；
- 全部 npz、`STATUS_REGISTRY.md`、WARN/NOT_APPLICABLE 等状态原样保留；
- 其它输出（Pt、benchmarks、scan/converge/spin_texture 图）不变。

## 3. 测试

本次改动只涉及绘图函数（`cli._plot_bands/_plot_soc/_plot_spin_pairs` 与
`plot_style.FONTS`），数据路径未变。因此本版用**定向测试**验证受影响代码路径：

```bat
python -m unittest discover -s tests -p "test_stage_c.py"      :: soc CLI（含 PbNV）
python -m unittest discover -s tests -p "test_stage_a.py"      :: bands CLI / validate
python -m unittest discover -s tests -p "test_release_regressions.py"
```

结果记录在 `RELEASE/_v0.1.4_tests.txt`。完整 74 项套件最后一次全绿为 v0.1.3-rc1
（commit `f462abd`，426 s）；自那以后仅绘图函数与版本字符串变化，未触碰
数据/求解/分析路径。

## 4. 与历史版本关系

| 版本 | 内容 | tag |
|---|---|---|
| v0.1.0-rc1 | 初始物理冻结基线 | `v0.1.0-rc1` |
| v0.1.1-rc1 | 清理 outputs 测试产物 | `v0.1.1-rc1` |
| v0.1.2-rc1 | 能带/SOC 图改 PRL 风格 | `v0.1.2-rc1` |
| v0.1.3-rc1 | SOC 劈裂只看 gap-edge 两对 | `v0.1.3-rc1` |
| v0.1.4-rc1 | 本版：能带连续 + 字号统一（纯绘图） | `v0.1.4-rc1` |
