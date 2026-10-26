# -*- coding: utf-8 -*-
"""Floquet solvers (stage B).

Implemented and validated against each other / analytic limits:
  * evolution.py  — direct time evolution U(T), quasienergies, Ebar, wspec;
  * sambe.py      — extended-space Sambe Hamiltonian (independent cross-check);
  * effective.py  — high-frequency van Vleck expansion (validity estimated).
"""
from .evolution import (floquet_modes, quasienergies, mean_energies,
                        fold_quasienergy)
from .sambe import (fourier_components, sambe_matrix, sambe_modes,
                    sambe_vs_evolution)
from .effective import (effective_hamiltonian, effective_bands,
                        effective_split_error, HighFrequencyNotValid)

__all__ = ["floquet_modes", "quasienergies", "mean_energies",
           "fold_quasienergy",
           "fourier_components", "sambe_matrix", "sambe_modes",
           "sambe_vs_evolution", "effective_hamiltonian", "effective_bands",
           "effective_split_error", "HighFrequencyNotValid"]
