# 配置文件参考（CONFIG_REFERENCE）

材料配置是一个 JSON 文件，路径作为 `--config` 传入。顶层键与子键**由代码校验**
（未知键直接报错，见 `floquet_soc/core/model_factory.py` 的 `_KNOWN_TOP_LEVEL_KEYS`
与 `_KNOWN_SUBKEYS`）。本页是那份 allowlist 的可读版，也是唯一需要维护的键表来源。

路径规则：`root` 与 `files`/`spin` 里的路径先做 `${VAR}` 与 `~` 展开；相对路径
相对**配置文件所在目录**（所以 `configs/pt.json` 的 `"../benchmarks_data/Pt"` 指向
仓库内 `benchmarks_data/Pt`）。复制配置到别处时请用绝对路径或 `${VAR}`。

## 顶层键

| 键 | 类型 | 必需 | 含义 |
|---|---|---|---|
| `name` | str | 是 | 材料名（输出目录与缓存键的一部分） |
| `kind` | str | 是 | `"wannier90"` 或 `"analytic"` |
| `root` | str | 否 | 输入文件基目录（相对配置目录或绝对/`${VAR}`） |
| `files` | dict | Wannier 时是 | 输入文件路径，见下 |
| `spin` | dict | 否 | 自旋算符来源，见下 |
| `drive` | dict | 否 | 驱动参数，见下 |
| `numerics` | dict | 否 | 数值参数，见下 |
| `num_electrons` | int | 否 | 占据带数（VBM/CBM 判定的依据） |
| `num_wann` | int | 否 | Wannier 轨道数（仅无 win 时兜底） |
| `kpath` | list | 否 | `[[label, [kx,ky,kz]], ...]` 高对称路径 |
| `nk_seg` | int | 否 | 每段采样点数（默认 20） |
| `band_window` | [lo, hi] | 否 | 带图能量窗口 (eV) |
| `kmesh` / `kmesh_pst` / `kmesh_spin_floquet` | dict | 否 | k 网格，见下 |
| `kplane` | str | 否 | 平面切法（如 `"xy"`） |
| `scan` / `converge` / `rashba_fit` / `soc` / `pst` / `spectrum` / `system` | dict | 否 | 各命令专用块，见下 |
| `jobs` | int | 否 | 默认并行 worker 数（≥1） |
| `floquet_match_tol` | float | 否 | Floquet 态匹配容差 (eV) |
| `spin_inplane_ratio` | float | 否 | 面内/面外自旋画点阈值 |
| `spin_plot_band` | int | 否 | 自旋纹理默认带指标 |

## `files`

`win` `hr` `chk` `eig` `spn` `centres` `amn` `mmn` `wsvec` `wout`
`bands_ref` `bands_ref_kpt` —— 均为文件路径（相对 `root`）。
`wsvec_use_ws_distance`: bool，显式覆盖 `wsvec` 文件的 `use_ws_distance` 语义。

模型构建优先级（见 `build_model`）：`win+hr` → `chk+hr` → `chk+eig`。

## `spin`

| 键 | 含义 |
|---|---|
| `source` | `"spn"`（Wannier90 .spn + chk/eig 旋转）或 `"pauli"`（2 带解析） |
| `file` / `chk` / `eig` / `u_mat` / `u_dis_mat` | 变换所需文件路径 |
| `notes` | 自描述注释 |

## `drive`

| 键 | 含义 | 默认 |
|---|---|---|
| `hw` | ħω (eV) | 0.4 |
| `gamma` | 无量纲驱动强度 e·A0·a_ref/ħ，≥0 | 0.0 |
| `polarization` | `circular` / `linear` | circular |
| `sign` | ±1 旋向 | +1 |
| `a_ref` | 参考晶格常数（Å），覆盖模型推断 | 模型 lattice 首个矢量 |
| `ratio` / `phi` / `e1` / `e2` / `ellipticity` / `amplitude` | 椭圆/线性偏振参数 | — |

## `numerics`

`n_step`（时间步数，≥1，默认 150）、`m_sp`（Floquet 谐波阶，≥2，默认 16）、
`pmax`（Peierls 贝塞尔截断阶，默认 Wannier 3 / 解析 8）。

## 网格块 `kmesh` / `kmesh_pst` / `kmesh_spin_floquet`

`n`（每方向点数）、`kmax`（|k| 上界）、或 `dk` + `kind`（自动 2π/Å 间距）。

## 命令专用块

- `scan`：`hw` `gamma`（列表）、`kpt`（[[kx,ky,kz], …]）、`n_step`、`m_sp`、`window`。
- `converge`：`kpt`、`n_step`/`m_sp`（列表）、`window`、`tol_ebar`、`tol_spin`、`tol_wsum`。
- `rashba_fit`：`anchor` `direction`（高对称点标签）、`nk`、`kmax`、`kmin`、`step`、`r2_thr`、`alpha_ratio_thr`、`split_r2_thr`。
- `soc`：`emin` `emax`（能窗）、`max_gap`。
- `pst`：`bands`（带指标列表）、`min_polarization`、`min_splitting`。
- `spectrum`：`nE`、`emin`、`emax`、`sigma`。
- `system`：`type`、`inversion`（bool）、`time_reversal`（bool）。
- `model`（仅 `kind="analytic"`）：`type`（`rashba`/`rashba_dresselhaus`）、`hbar2_over_2m`、`alpha`、`beta`。

## 最小示例

```json
{
  "name": "PbNV",
  "kind": "wannier90",
  "root": "../../PbNV",
  "files": {"win": "wannier90.win", "hr": "wannier90_hr.dat",
            "chk": "wannier90.chk", "eig": "wannier90.eig",
            "spn": "wannier90.spn"},
  "spin": {"source": "spn", "file": "wannier90.spn",
           "chk": "wannier90.chk", "eig": "wannier90.eig"},
  "num_electrons": 112,
  "drive": {"hw": 0.4, "gamma": 0.1, "polarization": "circular"}
}
```
