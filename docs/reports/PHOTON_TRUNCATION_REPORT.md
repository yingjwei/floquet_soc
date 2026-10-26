# Photon truncation report (m_sp convergence)

The previous convergence report used only global dWsum / matched deviations.
dEbar and dS are *mathematically independent* of the photon window by
construction, so they cannot demonstrate m_sp convergence. The report now
uses per-state photon-weight deficits  D_i = 1 - sum_{m=-m_sp}^{m_sp} w_i(m)
inside the target energy window.

Pt convergence run (k=(0.1,0.05,0), hbar*omega=1 eV, gamma=0.02, 240 steps,
window [-30,30] eV):

| m_sp | max D_i | median D_i | #states D_i>1e-3 | dWsum |
|---|---|---|---|---|
| 16 | **1.00** | 5.6e-9 | **2** | 1.34e-8 |
| 32 | 3.9e-14 | 1.2e-14 | 0 | 0 |
| 48 | 3.9e-14 | 1.2e-14 | 0 | 0 |

Interpretation: at m_sp=16 two states in the target window lose ~100%% of
their photon weight (the window is too narrow for those states), while the
median state is already converged. m_sp=32 gives full photon-sector weights
(accepted m_sp >= 32). Production setting [n_step=160, m_sp=32] is consistent
with the photon-deficit criterion.

Implication for SOC / spin-texture / PST analyses: any state with
D_i > 1e-3 must not be used for spin spectral weights or PST-style sums; the
new output (`w_deficit_max_m`, `w_deficit_median_m`, `n_deficit_gt_1e3`,
`accepted_m_sp` in `converge.npz`) makes this check explicit. The
spin-texture/PST runs themselves use converged observables (S_bar, S_t0 are
window-independent), so they are not affected by the m_sp=16 deficit; only
photon-resolved spectral quantities need m_sp >= 32.
