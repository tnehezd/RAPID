# T Tauri stars: indicative parameters

The values below are order-of-magnitude starting points for low-mass,
pre-main-sequence stars a few million years old. They are not strict limits
or values from a single stellar evolutionary track, and they do not replace
measurements of a specific object.

## Typical ranges

| Quantity | Indicative range | Notes |
|---|---|---|
| Stellar mass, M_star | approximately 0.2-2 M_sun | For low-mass T Tauri stars. |
| Effective temperature, T_eff | approximately 3000-5500 K | Depends on mass and age. |
| Stellar radius, R_star | approximately 1-4 R_sun | Changes as the star contracts. |
| Photospheric luminosity, L_phot | approximately 0.1-10 L_sun | Excludes accretion luminosity. |
| Accretion rate, Mdot_star | often 1e-9-1e-7 M_sun/yr | For accreting T Tauri stars outside outbursts; not a universal range. |

Values in different rows cannot be chosen independently. In particular,
R_star, T_eff and L_phot must satisfy:

```text
L_phot = 4*pi*R_star^2*sigma_SB*T_eff^4
```

The main-sequence mass-luminosity relation cannot be applied directly to
T Tauri stars: age and accretion state also matter. A mass of 1 M_sun does
not automatically imply a luminosity of 1 L_sun.

## An idealized, internally consistent example

| Quantity | Value | Role |
|---|---|---|
| M_star | 1 M_sun | Prescribed stellar mass. |
| R_star | 2 R_sun | Prescribed stellar radius. |
| T_eff | 4300 K | Prescribed effective temperature. |
| L_phot | approximately 1.23 L_sun | Calculated from the radius and T_eff. |
| Mdot_star | 1e-8 M_sun/yr | Prescribed accretion rate. |
| R_in | 5 R_star | Idealized starting radius of the accretion flow, not the inner boundary of the simulation grid. |
| L_acc | approximately 0.13 L_sun | From the potential-energy estimate below. |
| L_phot + L_acc | approximately 1.36 L_sun | Simple sum; the two components can have different spectra and absorption efficiencies. |

```text
L_acc ~= G*M_star*Mdot_star/R_star * (1 - R_star/R_in)
```

This is an illustrative parameter set, not a point taken from a stellar
evolution model at a specified age. It assumes that the accretion potential
energy given by the expression above is radiated away.

## Suggested RAPID irradiation parameter study

Keep the stellar mass at 1 M_sun and `irradiation_factor` at 0.05.
Vary the luminosity using
`thermodynamics_parameters.irradiation_luminosity_lsun`:

| L_star [L_sun] | YAML value | T_irr / T_irr(L_sun) | Change in irradiation temperature |
|---|---|---|---|
| 0.5 | `0.5` | 0.841 | -15.9% |
| 1 | `1.0` | 1 | 0% |
| 2 | `2.0` | 1.189 | +18.9% |

These ratios apply at the same radius and irradiation factor:

```text
T_irr^4 = irradiation_factor * L_star / (4*pi*sigma_SB*r^2)
T_irr proportional to L_star^(1/4)
```

These are not necessarily the changes in the total midplane temperature:
viscous heating and background radiation also contribute to the equilibrium.
The current model does not automatically calculate L_star from stellar mass
or accretion rate. The prescribed luminosity is a bolometric irradiation
parameter, not the L_X used for X-ray photoevaporation.

For details and limitations of the irradiation approximation, see
[Thermodynamics](thermodynamics.md#prescribed-stellar-irradiation).
