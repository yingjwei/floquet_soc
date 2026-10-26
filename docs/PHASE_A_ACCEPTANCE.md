# 阶段 A 验收报告（Phase A acceptance）

对象：MoSSe 15×15×1（num_wann=34, num_bands=64）与 PbNV 回归；核心 v0.1.7-rc1 工作区。

## 1. 输入文件事实

| 文件 | 事实 |
|---|---|
| `wannier90.win` | num_wann=34, num_bands=64, spinors=.true., mp_grid 15×15×1；`dis_win=[-63.36,10.21]`, `dis_froz=[-6.07,1.93]` |
| `wannier90_hr.dat` | 225 个 R（折叠），ndegen=1；`use_ws_distance=.true.`，`write_ndegen_applied` 未设置 |
| `wannier90_wsvec.dat` | 260100 blocks；11626 个非零超胞平移；展开后 261 R |
| `wannier90.spn` | VASP `vaspspn` 格式，**64×64 in WAVECAR band basis**，不是 Wannier 基 |
| `wannier90.chk` + `.eig` | 提供 `v(k)` 与本征值，用于把 spn 旋转到 34 维 Wannier 基 |
| `.centres.xyz` | 34 个 Wannier 中心，用于 Peierls 位移 |

## 2. 自旋算符构造与检查

1. 读取 spn：`num_bands=64 != num_wann=34` → 判定 `basis="vasp-bands"`，**不会**把 64×64 直接当 Wannier 自旋矩阵。
2. 旋转：`S_W(k)=v†(k) S_band(k) v(k)`（`core/spin_wannier.py`），使用 chk 的 `v`。
3. 一致性：自旋算符与 H(k) 使用同一套 Wannier90 R/ndegen + wsvec 展开（v0.1.6 起）。
4. 数值检查（MoSSe）：
   - `S_i` 厄米性：机器精度；
   - 离网格厄米性：**1.11e−15**；
   - on-grid 回代误差：2.4e−11；
   - H_W 与 `hr.dat` 对比：9e−6 eV（PbNV 回归）。
5. 自旋来源/哈希：`floquet.py provenance` 记录 spn/chk/eig/config 的 SHA256、版本、numerics。

## 3. 静态自旋纹理与 DFT 对照

- K 谷：Wannier `S=(0.0008,0.0014,−0.9999)`；DFT band 26 `(0,0,−0.78)` → **方向一致**，幅度为 DFT 归一化差异。
- 全网格方向：归一化 RMS **0.0082**（劈裂 >50 meV）；Sz 相关 **0.9991**（>5 meV）。
- 能量：DFT−Wannier 常数偏移 **−3.50 meV**，去偏移残差 RMS **0.73 meV**。
- 近 Γ 简并：93 个 DFT 点 `Sz=0`，单带自旋规范依赖，已从方向统计排除。
- **WARN**：DFT 文件 `|S|≤0.78`（约 0.73×Pauli），需确认 DFT 输出归一化；不据此缩放 Wannier。

## 4. 阶段 A 判定

| 项 | 判定 |
|---|---|
| 输入格式/尺寸/基组识别 | **PASS** |
| spn 基组判定 + chk 旋转 | **PASS** |
| 自旋算符厄米/维度/网格一致性 | **PASS** |
| 静态自旋期望/纹理 | **PASS** |
| 与 DFT 自旋交叉验证 | **PASS（方向）/ WARN（幅度）** |
| `.u.mat`/`.u_dis.mat` 兼容 | **GAP**（当前只支持 chk+eig；缺文件时明确报错） |
| 磁性子空间自旋算符 | **GAP**（未验证） |

**结论：阶段 A 通过（带 WARN）**；可以进入阶段 B/C 的解析模型验证，但真实磁性材料与 `.u.mat` 路径尚未验收。

## 5. 复现

```powershell
Set-Location D:\Floquet\floquet_soc
python floquet.py validate --config D:\Floquet\materials\MoSSe\configs\mosse.json
python floquet.py provenance --config D:\Floquet\materials\MoSSe\configs\mosse.json
python floquet.py dft-spin-check --config ... --dft D:\Floquet\MoSSe\dft\SPINTEXTURE_2D_GAMMA_CENTER.dat
```
