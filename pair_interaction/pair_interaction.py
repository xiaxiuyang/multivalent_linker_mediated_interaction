"""Analytical EE‑mediated pair‑interaction calculator.

This module implements the statistical‑mechanical formalism developed in

    X. Xia *et al.*,
    *Designed self‑assembly of programmable colloidal atom‑electron equivalents*,
    arXiv:2410.23784 (2024).

The central class :class:`PairInteraction` evaluates the effective free energy
and potential between two programmable atom equivalents (PAEs) connected by
multivalent electron‑equivalent (EE) linkers of arbitrary valency.
"""

from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np
from scipy.optimize import fsolve

from .constants import RHO_0
from . import utils


class PairInteraction:
    """Compute EE‑mediated thermodynamics for a single PAE–PAE pair.

    All energies are expressed in units of *k*\_B*T*.

    Parameters
    ----------
    n : Sequence[int]
        Number of grafted linkers on colloids **A** and **B** (``[n_A, n_B]``).
    lambda_ : Sequence[int]
        Linker valency on **A** and **B** (``[\\lambda_A, \\lambda_B]``).
        A valency ``λ=1`` corresponds to monovalent (one‑to‑one) binding, while
        larger values allow a single linker to bridge multiple sites.
    dG : Sequence[float]
        Standard hybridisation free energies for linker attachment to **A**
        and **B** (positive for unfavourable binding).
    mu : float
        Chemical potential of free linkers in bulk solution.
    radius : float
        Core radius of each PAE (nm).
    rc : float
        Contour length of a grafted ssDNA (nm). Defines the shell thickness.
    d : float
        Centre‑to‑centre distance between the two colloids (nm).

    Notes
    -----
    The implementation follows the notation of the Supplementary Information.
    The derivation assumes a grand‑canonical reservoir of linkers and an ideal
    gas approximation inside the accessible volume.
    """

    # --------------------------------------------------------------------- #
    # Public API                                                            #
    # --------------------------------------------------------------------- #

    def __init__(
        self,
        n: Sequence[int],
        lambda_: Sequence[int],
        dG: Sequence[float],
        mu: float,
        radius: float,
        rc: float,
        d: float,
    ):
        # Store geometric volumes that depend only on the pair separation *d*
        self.Va_prime: float = utils.shell_volume(radius, rc)
        self.Va: float = self._calc_Va(radius, rc, d)
        self.Vb: float = self._calc_Vb(radius, rc, d)

        self.n = np.asarray(n, dtype=float)  # grafted linker counts on A and B
        self.lambda_ = np.asarray(lambda_, dtype=int)  # valencies

        # Chi factors encode the Boltzmann weight of forming a bond in the
        # *current* (Va) or reference (Va_prime) geometry.
        dG = np.asarray(dG, dtype=float)
        self.chi = np.exp(-(dG + np.log(RHO_0 * self.Va)))
        self.chi_prime = np.exp(-(dG + np.log(RHO_0 * self.Va_prime)))

        # Grand‑canonical fugacities for linkers bound to A, B, or unbound.
        exp_mu = np.exp(mu)
        self.m0_A = np.array([exp_mu * self.Va, exp_mu * self.Va])
        self.m0_B = exp_mu * self.Vb
        self.m0_A_prime = np.array([exp_mu * self.Va_prime, exp_mu * self.Va_prime])

        # Solve the self‑consistent equations for mean linker numbers.
        self._solve_barn_prime()
        self._solve_barn()

        # Free energy components
        self.M = np.array(
            [
                self._sum_m(self.barn[0], self.lambda_[0], self.chi[0], self.m0_A[0]),
                self._sum_m(self.barn[1], self.lambda_[1], self.chi[1], self.m0_A[1]),
            ]
        )
        self.M_prime = np.array(
            [
                self._sum_m(
                    self.barn_prime[0],
                    self.lambda_[0],
                    self.chi_prime[0],
                    self.m0_A_prime[0],
                ),
                self._sum_m(
                    self.barn_prime[1],
                    self.lambda_[1],
                    self.chi_prime[1],
                    self.m0_A_prime[1],
                ),
            ]
        )
        self.Q = np.array(
            [
                self._sum_q(
                    self.barn,
                    self.lambda_,
                    self.chi,
                    self.m0_A,
                    self.m0_B,
                )
            ]
        )
        # Hard‑core repulsion between grafted tethers (see Eq. S16)
        self.F_rep = self._calc_F_rep(self.n[0], self.n[1], self.Va, self.Va_prime)

        # Absolute free energies
        self.F = (
            self.n[0] * (np.log(self.barn[0] / self.n[0]) + 1)
            + self.n[1] * (np.log(self.barn[1] / self.n[1]) + 1)
            - self.barn[0]
            - self.barn[1]
            - self.M[0]
            - self.M[1]
            - self.Q
        )
        self.F_ref = (
            self.n[0] * (np.log(self.barn_prime[0] / self.n[0]) + 1)
            + self.n[1] * (np.log(self.barn_prime[1] / self.n[1]) + 1)
            - self.barn_prime[0]
            - self.barn_prime[1]
            - self.M_prime[0]
            - self.M_prime[1]
        )

        # Effective potential of mean force (reduced units)
        self.U_eff = self.F - self.F_ref + self.F_rep

    # ------------------------------------------------------------------ #
    # Convenience properties                                             #
    # ------------------------------------------------------------------ #

    @property
    def effective_potential(self) -> float:
        """Return *U*\_eff in *k*\_B*T* at the specified separation *d*."""
        return float(self.U_eff)

    # ------------------------------------------------------------------ #
    # Internal helpers                                                   #
    # ------------------------------------------------------------------ #

    # --------- combinatorial sums ------------------------------------ #

    @staticmethod
    def _sum_m(barn: float, lambda_: int, chi: float, m0: float) -> float:
        """Eq. (S10): ∑ₘ m C(λ,m) (barn χ)^m."""
        return m0 * ((barn * chi + 1) ** lambda_ - 1)

    @staticmethod
    def _sum_km(barn: float, lambda_: int, chi: float, m0: float) -> float:
        """Derivative of :func:`_sum_m` with respect to *barn*."""
        return m0 * lambda_ * barn * chi * (barn * chi + 1) ** (lambda_ - 1)

    @staticmethod
    def _sum_q(
        barn: Sequence[float],
        lambda_: Sequence[int],
        chi: Sequence[float],
        m0_A: Sequence[float],
        m0_B: float,
    ) -> float:
        """Eq. (S14): cross‑linker contribution."""
        return (
            m0_B
            / (m0_A[0] * m0_A[1])
            * PairInteraction._sum_m(barn[0], lambda_[0], chi[0], m0_A[0])
            * PairInteraction._sum_m(barn[1], lambda_[1], chi[1], m0_A[1])
        )

    @staticmethod
    def _sum_kq1(
        barn: Sequence[float],
        lambda_: Sequence[int],
        chi: Sequence[float],
        m0_A: Sequence[float],
        m0_B: float,
    ) -> float:
        """∂Q/∂barn₁ (see Supplementary)."""
        return (
            m0_B
            / (m0_A[0] * m0_A[1])
            * PairInteraction._sum_km(barn[0], lambda_[0], chi[0], m0_A[0])
            * PairInteraction._sum_m(barn[1], lambda_[1], chi[1], m0_A[1])
        )

    @staticmethod
    def _sum_kq2(
        barn: Sequence[float],
        lambda_: Sequence[int],
        chi: Sequence[float],
        m0_A: Sequence[float],
        m0_B: float,
    ) -> float:
        """∂Q/∂barn₂ (see Supplementary)."""
        return (
            m0_B
            / (m0_A[0] * m0_A[1])
            * PairInteraction._sum_m(barn[0], lambda_[0], chi[0], m0_A[0])
            * PairInteraction._sum_km(barn[1], lambda_[1], chi[1], m0_A[1])
        )

    # --------- self‑consistent solutions ----------------------------- #

    def _solve_barn(self) -> None:
        """Solve Eq. (S12) for mean number of *A*‑attached linkers."""

        def _sc(barn):
            return [
                barn[0]
                + self._sum_km(barn[0], self.lambda_[0], self.chi[0], self.m0_A[0])
                + self._sum_kq1(barn, self.lambda_, self.chi, self.m0_A, self.m0_B)
                - self.n[0],
                barn[1]
                + self._sum_km(barn[1], self.lambda_[1], self.chi[1], self.m0_A[1])
                + self._sum_kq2(barn, self.lambda_, self.chi, self.m0_A, self.m0_B)
                - self.n[1],
            ]

        # Use half of grafted linkers as an initial guess
        self.barn = fsolve(_sc, self.n / 2.0)

    def _solve_barn_prime(self) -> None:
        """Solve the reference (non‑overlapping) state."""

        def _sc(barn):
            return [
                barn[0]
                + self._sum_km(
                    barn[0], self.lambda_[0], self.chi_prime[0], self.m0_A_prime[0]
                )
                - self.n[0],
                barn[1]
                + self._sum_km(
                    barn[1], self.lambda_[1], self.chi_prime[1], self.m0_A_prime[1]
                )
                - self.n[1],
            ]

        self.barn_prime = fsolve(_sc, self.n / 2.0)

    # --------- geometry --------------------------------------------- #

    @staticmethod
    def _calc_Va_prime(radius: float, rc: float) -> float:
        """Volume accessible to a linker grafted on an *isolated* colloid."""
        return utils.shell_volume(radius, rc)

    def _calc_Va(self, radius: float, rc: float, d: float) -> float:
        """Accessible volume of a tether on **A** in presence of **B**."""
        r1, r2 = radius, radius + rc
        overlap = utils.overlap_sphere(r1, r2, d)
        return self._calc_Va_prime(radius, rc) - overlap

    def _calc_Vb(self, radius: float, rc: float, d: float) -> float:
        """Volume available to an EE linker bridging **A** and **B**."""
        shell_vol = utils.overlap_sphere(radius + rc, radius + rc, d)
        excluded = utils.overlap_sphere(radius + rc, radius, d)
        return shell_vol - 2.0 * excluded

    # --------- steric repulsion ------------------------------------- #

    @staticmethod
    def _calc_F_rep(nl: int, nr: int, Va: float, Va_prime: float) -> float:
        """Repulsive free‑energy of grafted polymers (Eq. S16)."""
        return -(nl + nr) * np.log(Va / Va_prime)