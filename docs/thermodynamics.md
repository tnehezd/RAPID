# Thermodynamics

## Configuration and running

The optional YAML section `thermodynamics_parameters` selects the model.
Omitting it, or its `model` key, preserves existing locally isothermal behavior.
`lociso` and `locally_isotherm` are aliases: the prescribed aspect ratio and
flaring determine the sound speed and pressure. `polytropic` is not implemented
and is rejected, as are unknown model names and settings.

For radiative equilibrium:

```yaml
thermodynamics_parameters:
  model: radiative
  mean_molecular_weight: 2.34
  background_temperature_k: 10.0
  minimum_temperature_k: 1.0
  maximum_temperature_k: 1500.0
  irradiation_luminosity_lsun: 0.0   # default: off; use e.g. 1.0 for irradiation
  irradiation_factor: 0.05          # effective geometry/absorption factor, [0, 1]
```

Numerical parameters are optional, with the defaults above. They may remain
when switching to `lociso`; the wrapper reports that they are inactive.
Run the [example configuration](../configs/config_thermal.yaml) from the project
root:

```bash
make all
python run_simulation.py -c configs/config_thermal.yaml
```

The C CLI accepts `-thermal_model`, `-thermal_mu`, `-thermal_background`,
`-thermal_tmin`, and `-thermal_tmax`. Existing internal benchmark modes require
locally isothermal thermodynamics. Linear and logarithmic grids work, but ghost
radii must be positive; the example uses a logarithmic grid at its 0.1 AU inner
edge.

The startup terminal summary includes a `THERMODYNAMICS` panel next to the
boundary-condition summary. It reports the active model, prescribed geometry
for locally isothermal runs, or heating/cooling, opacity, molecular weight,
background temperature and solver bracket for radiative runs.
It also identifies instantaneous equilibrium and whether prescribed stellar
irradiation is active, including luminosity and irradiation factor when active.
A separate `SELF-GRAVITY` panel reports whether SG is enabled, its axisymmetric
rotation model and, for radiative runs with SG, the thermal coupling.

## Equations and integration

The code-unit EOS and vertical closure are in [eos.c](../src/physics/eos.c),
with the public interface in [eos.h](../include/physics/eos.h).
`calculateEOS(Sigma, r, disk_params)` explicitly dispatches between locally
isothermal and radiative closures: the lociso branch computes prescribed H(r)
and cs = Omega*H, while the radiative branch reads equilibrium H and cs from
the thermal cache. It then evaluates Gaussian midplane density and returns
P_mid = rho_mid*cs^2.
It does not evolve internal energy. The existing
gas_physics.h interface includes eos.h for compatibility.

Parabolic surface-density boundary extrapolation uses the adjacent physical
density (zero gradient) when a ghost value would undershoot `density_floor`.
A warning is printed on the
first activation at each boundary. This applies during initialization and
evolution, in both thermodynamic modes, without changing physical cells or
the extrapolation of pressure, pressure gradient or velocity. Absorbing density
boundaries still prescribe zero. This limiter changes the boundary flux when
active; boundary sensitivity must still be checked for scientific results.
Run `make test-boundary` for regression tests covering negative extrapolated
ghosts, both grid types and boundaries, thermal closure and unchanged other
fields/boundary types.

The independent [CGS solver](../src/physics/thermal_model.c) is connected to gas
evolution through [disk_thermal.c](../src/physics/disk_thermal.c). Each trial
temperature determines the isothermal sound speed, scale height, density,
viscosity, opacity and optical depth:

```text
cs^2 = k_B*T/(mu*m_H)
H = cs/Omega
rho = Sigma/(sqrt(2*pi)*H)
nu = alpha*cs*H
q = |d ln Omega / d ln r|
Qplus = nu*Sigma*(q*Omega)^2
tau = kappa*Sigma/2
tau_eff = 3*tau/8 + sqrt(3)/4 + 1/(4*tau)
T_irr^4 = irradiation_factor * L_star/(4*pi*sigma_SB*r^2)
T_rad^4 = T_bg^4 + T_irr^4
Qminus = 2*sigma_SB*(T^4 - T_rad^4)/tau_eff
Qminus = Qplus
```

### Prescribed stellar irradiation

Set `irradiation_luminosity_lsun` to the stellar luminosity in solar units
(`L_sun = 3.828e33 erg/s`); zero disables irradiation and preserves the
previous model. `irradiation_factor` defaults to 0.05, must be in [0, 1],
and combines the prescribed geometry, absorption and reprocessing efficiency.
Setting it to zero also disables irradiation. Both settings are inactive in
locally isothermal mode. Direct CLI equivalents are `-irradiation_luminosity`
and `-irradiation_factor`.

This is an **effective midplane radiation-temperature prescription**, not a
calculation of where stellar photons are absorbed. Its normalization is defined
by the equation above: no additional factor of two or albedo is applied.
With constant luminosity/factor, T_irr scales as r^(-1/2). The background and
irradiation temperatures add in fourth powers, not linearly. Viscous heating
raises T above this radiation floor; with zero viscosity T equals T_rad.
The thermal/SG iteration uses this same temperature, so irradiation affects
pressure, scale height, viscosity, gas evolution and dust dynamics.

`Qplus` remains **viscous heating only** and `Qminus` is cooling **net of both
background and irradiation**. Thus the existing nine-column snapshots and
energy-balance plots retain Qplus = Qminus; they do not separately plot the
intercepted stellar luminosity. Thermal snapshot headers and run metadata
record the irradiation parameters.

There is no self-shadowing, evolving grazing angle, separate stellar-wavelength
opacity or vertical radiative transfer. In particular, this prescription does
not reduce stellar absorption when the disk becomes transparent; use caution
for depleted/optically thin disks. It is a controlled first approximation,
not a photon-energy-conserving irradiation solver.

Fluxes include both disk faces. The bridge uses RAPID's AU, solar-mass and
`year/(2*pi)` units for conversion to CGS and back. The standalone example uses
its own explicitly stated physical constants. The low-level API requires CGS
inputs, temperature in K and dimensionless `q` (1.5 for a Keplerian disk).

Equilibrium is recomputed at initialization and after each gas-density update,
before refreshing pressure and its gradient. Viscosity is frozen within each
existing `nu*Sigma` update, then recomputed. Radiative mode also handles zero
local viscosity without division by zero. Timesteps, dust drift/growth and
pressure-trap widths use the thermally determined sound speed and scale height.
Dead-zone widths retain their prescribed values.

With `simulation_parameters.enable_selfgravity: true`, the local frequency and
shear include the radial SG field. Without SG, heating reduces to
`(9/4)*nu*Sigma*Omega_K^2`. The stellar shear is differentiated analytically.
SG softening depends on H, so SG and thermal states are iterated with
under-relaxation, relative H tolerance 1e-8 and at most 64 iterations.

The temperature solve uses bounded bisection, relative temperature and
heating/cooling tolerances of 1e-8, and at most 128 iterations. An already
balanced bracket endpoint is also accepted. The lower bound is raised to
`T_rad` (or `T_bg` without irradiation). The upper bound must accommodate
the irradiated temperature, including ghost radii. A missing sign change or non-convergence is an explicit error, not a
fallback to locally isothermal behavior. Multiple equilibria are possible:
callers must isolate the desired branch with the bounds; the solver neither
enumerates roots nor determines thermal stability.

## Opacity and limitations

The log-space opacity implementation uses smoothed Lin-Papaloizou power laws,
cross-checked against the
[FARGO-CPT Lin implementation](https://github.com/rometsch/fargocpt/blob/master/src/opacity.cpp).
It is not Bell-Lin (1994), and does not reuse the historical code's converted
prefactors. It retains the high-temperature electron-scattering branch and
the reference implementation's low-density exception at `rho <= 1e-10 g/cm^3`.
The blends soften transitions; branch switches are not guaranteed exactly
continuous.

This is instantaneous local equilibrium, not an energy-evolution or full
radiative-transfer solver. There is no self-consistent stellar irradiation geometry, radial
thermal transport or temperature-dependent molecular weight. The bridge uses
the same opacity for thick/thin terms, rather than separate Rosseland/Planck
means. Its numerical domain (`rho > 0`, `0 < T < 1e9 K`) does not imply physical
validity throughout; constant molecular weight cannot model dissociation or
ionization.

Using radial Omega in `H=cs/Omega` is an effective closure, **not** a solution of
vertical SG. Gas transport retains RAPID's existing Keplerian viscous operator;
generalized heating does not make it a full non-Keplerian angular-momentum
solver. Absorbing ghosts use the positive gas floor only for their thermal
closure, without changing the mass boundary.

## Tests and output

From the project root:

```bash
make test-thermal
python -m unittest discover -s tests -p 'test_thermal_config.py'
# Requires make all; exercises evolution, SG, grids, dust, photoevaporation, I/O.
python -m unittest discover -s tests -p 'test_radiative_simulation.py'
```

Generate the synthetic 1--100 AU standalone example (not YAML-driven):

```bash
./bin/test_thermal_model --profile > temperature_profile.dat
python plot_temperature.py temperature_profile.dat
```

It uses 1 Msun, `Sigma=100*(r/AU)^-1 g/cm^2`, `alpha=1e-3`, `mu=2.34`,
`T_bg=10 K` and Keplerian shear.

Radiative runs write ASCII files
`<run>/LOGS/temperature_profile_00000000.dat` initially, at scheduled
snapshots and at the final state, for both ASCII and HDF5 runs. The suffix uses
the same eight-digit, zero-padded integer output-time label as gas snapshots
(for example, `temperature_profile_00010000.dat`). The actual simulation time
in years is retained in the header. Existing gas
output schemas are unchanged. The run-information log records the thermal
model and parameters. The nine columns are:

```text
r_AU T_K H_over_r rho_g_cm3 nu_cm2_s kappa_cm2_g tau Qplus Qminus
```

They can be plotted directly:

```bash
python plot_temperature.py output_thermal/LOGS/temperature_profile_00000000.dat
MPLBACKEND=Agg python plot_temperature.py temperature_profile.dat -o temperature_profile.png
```

For line plots of T(r) at all saved times and the percentage change from the
earliest profile, with a shared time colorbar:

```bash
python plot_temperature_evolution.py output_thermal_sg
python plot_temperature_evolution.py output_thermal_sg -o temperature_evolution.png
```

This uses thermal snapshots only. The upper panel shows absolute temperature;
the lower panel shows 100*(T/T(initial)-1) [%], without clipping at +/-3%.
Dashed black curves identify the initial profile and thick curves the final one.

Combined time-radius heatmaps from matching gas and thermal ASCII snapshots:

```bash
python plot_thermal_heatmap.py output_thermal
MPLBACKEND=Agg python plot_thermal_heatmap.py output_thermal -o thermal_heatmaps.png
```

Nine panels in a 3-by-3 layout show Sigma/Sigma(initial), symlog dP/dr,
100*(T/T(initial)-1) [%], H/r, opacity, optical depth, heating, cooling and
100*(cs/cs(initial)-1) [%]. Sound speed is recovered from the midplane pressure and
density as sqrt(P/rho), with gas pressure converted from code units to CGS.
The relative equilibrium residual is omitted.
Heating and cooling share a color normalization. Positive fields use logarithmic
colors; fluxes containing zeros use linear colors. Temperature and sound-speed changes use a
linear diverging color scale centered on zero (negative for cooling, positive
for increases), fixed at -3% to +3%; values outside this range saturate at the
endpoint colors without changing the data. `--linthresh` sets the
pressure-gradient linear threshold (default 1e-12). White dP/dr=0 contours and
black positive-to-negative crossing markers follow the gas heatmap convention.
At least two paired snapshots are required; times and radial grids are checked.
Initial density, temperature and sound speed mean the earliest snapshot in that run,
normalized separately at each radius. Other thermal panels remain absolute.
