"""Minimal radiative-transfer calculation without FITS I/O."""

import numpy as np

from gtlmapping import compute_extinction, convert_surface_density

observed = np.array([[60.0]])
foreground = np.array([[20.0]])
background = np.array([[100.0]])

tau, sigma, saturated, invalid_background, bright = compute_extinction(
    observed,
    background,
    foreground,
    kappa_cm2_g=7.5,
)
sigma_msun_pc2 = convert_surface_density(sigma)

print("tau:", tau[0, 0])
print("Sigma [g/cm2]:", sigma[0, 0])
print("Sigma [Msun/pc2]:", sigma_msun_pc2[0, 0])
print(
    "saturated:",
    saturated[0, 0],
    "invalid background:",
    invalid_background[0, 0],
    "bright:",
    bright[0, 0],
)
