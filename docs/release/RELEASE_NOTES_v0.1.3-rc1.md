# RELEASE NOTES — v0.1.3-rc1（SOC 劈裂口径：只看 gap-edge 两对）

日期：2026-10-23　｜　基线：v0.1.2-rc1（commit `7682d0c`，figure-style release）
历史 tag `v0.1.0-rc1`、`v0.1.1-rc1`、`v0.1.2-rc1` 均保留且不可改写。

## 1. 这个版本改了什么

**SOC 劈裂分析口径 + 两张 SOC 图，按要求收敛到带边两对；Flouqet/能带计算核心不变。**

1. **新增 gap-edge 配对定义**（`observables/soc_splitting.py`）：
   - `gap_edge_pairs(num_electrons, num_wann)` → 0-based 对：
     **价带 `(VBM-1, VBM) = (ne-2, ne-1)`**、**导带 `(CBM, CBM+1) = (ne, ne+1)`**；
   - `edge_pair_splitting(energies, pairs, spins)` → 逐 k 直接能差 + 自旋反平行度
     （不做全窗口贪心配对）。
2. **`floquet.py soc` 改为只报这两对**：
   - 不再对 SOC 窗口内所有相邻带贪心配对（旧行为会带出 114-115 等与带边无关的对）；
   - 打印：两对的 max / median / spin-antiparallel median；
   - npz 新增 `edge_pairs`、`edge_pair_split(nk×2)`、`edge_pair_rows`；
     `spin_pair_rows` 现在只含这两对（旧格式 5 列保留）；
   - PT 体系（`pt_nonmagnetic`）分支不变（仅 Kramers 诊断，无 `spin_pair_rows`）。
3. **两张图重做**（PRL 风格，仅两对）：
   - `soc_splitting.png`：2 面板 —— (a) 带边能带（价带对蓝、导带对红，其余灰，
     能量零点 VBM，标注 min direct gap）；(b) 两对劈裂 ΔE(k)（meV）；
   - `soc_spin_splitting.png`：单面板 —— 两对劈裂曲线（含 max/median 图例），
     k 轴为累积 k 距离，Γ/X/S/Y 竖虚线。
4. **重新生成 PbNV 输出**：`outputs/PbNV/soc_splitting.{png,npz}`、
   `outputs/PbNV/soc_spin_splitting.png`；`outputs/PbNV/README.md` 更新为新口径。
5. **PbNV 新数值**（81 k，0-based 带号）：

   | 配对 | 含义 | max | median |
   |---|---|---|---|
   | 110-111 | VBM-1 / VBM（价带边） | **29.98 meV** | 0.51 meV |
   | 112-113 | CBM / CBM+1（导带边） | **107.66 meV** | 0.58 meV |

   对比旧贪心配对口径：112-113 为 103.28 meV、另报 114-115 71.98 meV。
   **旧值属 v0.1.0-rc1 历史口径**，在 `STAGE_C_REPORT.md` 中保留并加注；
   新值以 `edge_pair_rows`/`spin_pair_rows` 为准。

## 2. 什么没有变

- `floquet_soc/{core,drives,solvers,benchmarks,workflows}/` 零改动（`observables` 只新增
  两个纯分析函数，不改变已有 `spin_paired_splitting` 的行为/返回值）；
- Floquet 求解、Wannier 解析、自旋算符、PST、能带计算均未动；
- 其它输出（Pt、benchmarks、scan/converge/spin_texture 图）内容与 rc3 相同；
- 全部 WARN / NOT_APPLICABLE / NOT_VALIDATED / UNCONVERGED / BLOCKED / NEEDS_INPUT
  状态保持登记不变。

## 3. 复现

```bat
python floquet.py soc --config configs/pbnv.json --outdir outputs/PbNV
:: 打印两对 gap-edge 劈裂；生成 soc_splitting.png / soc_spin_splitting.png
python -m unittest discover -s tests -v      :: 结果见 RELEASE/_v0.1.3_test_run.txt
```

## 4. 与历史版本关系

| 版本 | 内容 | tag |
|---|---|---|
| v0.1.0-rc1 | 初始物理冻结基线 | `v0.1.0-rc1`（保留） |
| v0.1.1-rc1 | 清理 outputs 测试产物 | `v0.1.1-rc1`（保留） |
| v0.1.2-rc1 | 能带/SOC 图改 PRL 风格 | `v0.1.2-rc1`（保留） |
| v0.1.3-rc1 | 本版：SOC 劈裂只看 VBM-1/VBM 与 CBM/CBM+1 | `v0.1.3-rc1` |
