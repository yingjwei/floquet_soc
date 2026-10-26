# 错误排查指南（TROUBLESHOOTING）

异常、警告与退出码的含义和应对。所有解析/数值失败都遵循 fail-fast 原则：
宁可抛错，不静默给出错数。

## 退出码

| 码 | 含义 |
|---|---|
| `0` | 成功 |
| `2` | 输入/配置缺失或不可用（BLOCKED），打印了原因与缺失文件 |
| `3` | 命令未实现（历史遗留） |

## 异常（`floquet_soc.core.wannier`）

- **`WannierInputError`**：win/hr/chk/eig/wsvec 解析失败或格式矛盾
  （截断、维度不匹配、未知元素、自相矛盾的 wsvec 头）。→ 检查文件是否完整、
  `files` 路径是否正确；hr.dat 每块必须恰好 `nw*nw` 个矩阵元。
- **`HermiticityError`**：`H(k)` 非厄米超过相对容差 `1e-10`。→ 检查 `_hr.dat`
  的 R/-R 配对，或 `wannier90_wsvec.dat` 是否只给了单边块集。
  **不要**用 `0.5*(H+H†)` 掩盖——那会给出看似正常但错误的本征值。

## 异常（`floquet_soc.core.spin_operator`）

- **`SpinInputError`**：自旋源未知/不匹配（如 `pauli` 用于非 2 带模型）。
- **`MissingTransformError`**：需要 chk/eig 旋转但文件缺失。→ 补 `spin.chk`/`spin.eig`，
  或确认 `spin.source` 正确。

## 异常（`floquet_soc.solvers.effective`）

- **`HighFrequencyNotValid`**：van Vleck 展开的 1/ω 修正不满足高频条件
  （与能级间距同量级）。→ 提高 `drive.hw`，或用 `check_validity=False` 显式忽略
  （仅在你确定高频极限成立时）。

## 警告

- **`TruncationWarning`**（`drives/peierls`）：Peierls 贝塞尔级数在 `pmax` 处
  截断，丢弃的 `|J_p(g)|` 超过 1e-3。→ 提高 `numerics.pmax`/`m_sp`，或降低
  `drive.gamma`。同进程内同量级只告警一次。
- **`validate` 的 `WARN` 项**：数据质量问题（自旋求和规则 <3、TRS 残差、
  离网格插值未验证等）。`WARN` 不阻塞，但报告数值应一并引用。

## 配置校验失败（`core/model_factory.validate_cfg`）

- `unknown config key(s): ...` → 顶层键拼写错误，或用了已废弃键。
- `unknown key(s) in config['<block>']: ...` → 子键拼写错误（如 `sorce` vs `source`）。
- `drive.gamma must be >= 0` / `numerics.m_sp must be >= 2` / `jobs must be >= 1` → 数值越界。

若某键确实合法，把它加入 `_KNOWN_TOP_LEVEL_KEYS` / `_KNOWN_SUBKEYS` 并同步本页与
`docs/CONFIG_REFERENCE.md`。

## 常见症状

| 症状 | 可能原因 |
|---|---|
| `Wannier model needs one of: win+hr, chk+hr, or chk+eig` | `root` 相对路径解析到了错误目录（复制配置后最常见）；`files` 文件名拼写 |
| `missing files: ...`（BLOCKED） | `spin` 变换文件路径错；`spin.source` 与文件不符 |
| 离网格 H(k)/自旋非厄米 | wsvec 头未识别 → 静默退回折叠约定；显式设 `files.wsvec_use_ws_distance` |
| 并行 worker 失败但看不到原因 | 现在 `run_chunks` 会把 `log_i.txt` 尾部打进报错信息 |
| 图字体回退（Liberation/DejaVu） | Linux 无 Arial；只影响外观不影响数值 |

## 报告与溯源

每个命令尽量在输出目录写 `provenance.json` / `validate.json`（版本、配置哈希、
输入文件 sha256、数值设置）。报障时附上这两个文件 + 完整命令 + 报错输出，
通常可直接定位。
