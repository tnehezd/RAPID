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

## Equations and integration

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
Qminus = 2*sigma_SB*(T^4 - T_bg^4)/tau_eff
Qminus = Qplus
```

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
`T_bg`. A missing sign change or non-convergence is an explicit error, not a
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
radiative-transfer solver. There is no stellar irradiation geometry, radial
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
