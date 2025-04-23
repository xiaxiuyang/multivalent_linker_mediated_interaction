"""Physical and numerical constants used throughout the ``pae_ee_interaction`` package."""

# Avogadro‑like number that converts molar concentration (mol L⁻¹) to number
# density in nm⁻³ (1 M ≈ 0.6022 nm⁻³).
RHO_0 = 0.6022

# Boltzmann constant in reduced units: all free energies are expressed in k_BT,
# therefore k_B ≡ 1. This symbol is kept for clarity and future extensions where
# an explicit energy scale might be desirable.
K_B = 1.0