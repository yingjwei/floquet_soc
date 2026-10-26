# RELEASE NOTES — v0.1.1-rc1（cleanup-only release）

日期：2026-10-23　｜　基线：v0.1.0-rc1（commit `a9b8e9a…`，tag `v0.1.0-rc1`，保留且不可改写）

## 1. 这个版本改了什么

**只做发布卫生清理，物理代码与数值结果与 v0.1.0-rc1 完全相同。**

1. **从 `outputs/` 移除全部测试运行产物**：`outputs/_tests/`、`_tests_b/`、`_tests_c/`、
   `_tests_d/`、`_tests_spin/`。原因：这些是 `unittest` 运行产物，部分文件与材料主结果
   **同名**（例如 `_tests_c/pbnv/soc_spin_splitting.png` 用 `--nk 5` 粗 k 路径，
   峰值 ~88 meV；主结果 `outputs/PbNV/soc_spin_splitting.png` 用 `--nk 20`，
   峰值 103.28 meV），会被误当主结果引用。
2. **`.gitignore` 增加 `outputs/_tests*/`**：测试产物不再进入版本库与发布快照；
   `unittest` 仍会照常生成它们（在本地工作目录，便于调试），但不再污染发布输出。
3. **版本号**：`VERSION` / `floquet_soc/version.py` / `tests/test_release_regressions.py`
   由 `0.1.0-rc1` → `0.1.1-rc1`；`python floquet.py --version` 输出相应变化。
4. `outputs/_cache/`（可再生成模型缓存，约 159–323 MB）不进入快照（rc1 起已被
   `.gitignore` 忽略；本次同时清理了工作目录中的缓存副本）。

## 2. 什么没有变

- `floquet_soc/{core,drives,solvers,observables,benchmarks,workflows}/` **零改动**；
- 全部材料主结果（`outputs/PbNV/`、`outputs/Pt/`、`outputs/Rashba/`、
  `outputs/benchmarks/`、`outputs/*_verify/`）与状态登记（`STATUS_REGISTRY.md`）
  内容不变；WARN / NOT_APPLICABLE / NOT_VALIDATED / UNCONVERGED / BLOCKED /
  NEEDS_INPUT 全部保留；
- 测试内容除版本断言外不变：`python -m unittest discover -s tests -v` →
  **Ran 74 tests, OK (skipped=2)**，0 failure / 0 error（在 v0.1.1-rc1 上复跑）。

## 3. 与 rc1 的关系

| 项 | v0.1.0-rc1 | v0.1.1-rc1 |
|---|---|---|
| git tag | `v0.1.0-rc1`（保留） | `v0.1.1-rc1`（本版） |
| 物理代码 | 基线 | 完全相同 |
| `outputs/_tests*` | 随快照包含（4.6 MB） | **已移除**（可由测试再生成） |
| 版本字符串 | 0.1.0-rc1 | 0.1.1-rc1 |
| 快照位置 | `RELEASE/floquet_soc-v0.1.0-rc1/` | `RELEASE/floquet_soc-v0.1.1-rc1/` |

引用旧结果时用 rc1 tag；做新一轮材料研究时用 rc1 或 rc2 的代码均可（代码相同），
但**引用 outputs 结果时以 rc2 快照为准**（不含测试产物）。

## 4. 复现

```bat
git -C <repo> archive --format=tar v0.1.1-rc1 -o snapshot.tar
python floquet.py --version      :: floquet_soc 0.1.1-rc1
python -m unittest discover -s tests -v   :: Ran 74, OK (skipped=2)
```
