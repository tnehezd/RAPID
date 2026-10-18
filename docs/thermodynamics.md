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
For time-dependent internal energy, select `model: radiative_energy` and
optionally set `gamma: 1.4` (default; finite and strictly greater than 1).
The existing `radiative` mode retains instantaneous equilibrium.
Run the [example configuration](../configs/config_thermal.yaml) from the project
root:

```bash
make all
python run_simulation.py -c configs/config_thermal.yaml
```

The C CLI accepts `-thermal_model`, `-thermal_mu`, `-thermal_gamma`, `-thermal_background`,
`-thermal_tmin`, and `-thermal_tmax`. Existing internal benchmark modes require
locally isothermal thermodynamics. Linear and logarithmic grids work, but ghost
radii must be positive; the example uses a logarithmic grid at its 0.1 AU inner
edge.

The startup terminal summary includes a `THERMODYNAMICS` panel next to the
boundary-condition summary. It reports the active model, prescribed geometry
for locally isothermal runs, or heating/cooling, opacity, molecular weight,
background temperature and solver bracket for radiative runs.
It also identifies instantaneous equilibrium versus energy evolution (including
gamma and the transport scheme), and whether prescribed stellar
irradiation is active, including luminosity and irradiation factor when active.
A separate `SELF-GRAVITY` panel reports whether SG is enabled, its axisymmetric
rotation model and, for radiative runs with SG, the thermal coupling.

## Equations and integration

The code-unit EOS and vertical closure are in [eos.c](../src/physics/eos.c),
with the public interface in [eos.h](../include/physics/eos.h).
An internal `calculateEOSState` helper selects between locally
isothermal and radiative closures: the lociso branch computes prescribed H(r)
and cs = Omega*H, while the radiative branch reads current H and cs from
the thermal cache. All public closure functions use this shared selection.
`calculateEOS(Sigma, r, disk_params)` evaluates Gaussian midplane density and
returns P_mid = rho_mid*cs^2. Height-only requests do not read or compute cs:
SG initialization uses an initial H guess before the first thermal solve.
The EOS itself does not advance internal energy. The existing
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
equilibrium energy-balance plots retain Qplus = Qminus; in energy mode their
difference is a physical source term, not a solver residual. They do not separately plot the
intercepted stellar luminosity. Thermal snapshot headers and run metadata
record the irradiation parameters.

In prescribed mode there is no self-shadowing, evolving grazing angle, separate stellar-wavelength
opacity or vertical radiative transfer. In particular, this prescription does
not reduce stellar absorption when the disk becomes transparent; use caution
for depleted/optically thin disks. It is a controlled first approximation,
not a photon-energy-conserving irradiation solver.

### Geometry-dependent irradiation

Both `radiative` and `radiative_energy` support the opt-in
`irradiation_model: geometry`. The default remains `prescribed`, so existing
configurations retain their normalization and behavior. For example:

```yaml
thermodynamics_parameters:
  model: radiative_energy       # also supported: radiative
  irradiation_luminosity_lsun: 1.0
  irradiation_model: geometry
  irradiation_surface_height: 4.0
  irradiation_stellar_radius_rsun: 2.0
```

The new positive parameters default to 4 and 2, respectively; CLI equivalents
are `-irradiation_model`, `-irradiation_surface_height` and
`-irradiation_stellar_radius`. The stellar radius must be smaller than every
physical grid radius. Luminosity zero still disables the stellar contribution.
The fixed `irradiation_factor` is **not used** in geometry mode.

The absorbing surface is prescribed as `z_s = chi*H`. Its radial slope is a
backward difference between adjacent physical cells. At the inner physical
cell a conical continuation (`dz_s/dr = z_s/r`) closes the unresolved interior;
there is no inner rim or disk material inside the computational domain.
For a surface at `(r,z_s)` with slope `s`, define
`g = (r*s-z_s)/(sqrt(r^2+z_s^2)*sqrt(1+s^2))`.
The star is a uniformly bright circular disk centered at the origin, with
small angular size. A stellar strip at height `y` has incidence
`max(0, g + y/(sqrt(r^2+z_s^2)*sqrt(1+s^2)))`. Rays from each strip to the
surface are blocked if they pass below **any** interior surface node.
The integral over visible, positively illuminating stellar strips gives the
intercepted fraction. Half is assigned to downward reprocessing, defining
`f_geom = 0.5*<positive visible incidence>`. The existing closure then uses
`T_irr^4 = f_geom*L_star/(4*pi*sigma_SB*r^2)`. Partial stellar occultation and
the finite stellar size avoid a binary point-source shadow transition.

Equilibrium cells are solved from inside out: each temperature trial
recalculates its own H, surface slope, occultation and irradiation, while the
interior solved surfaces remain fixed. This is a causal, one-sided spatial
closure, not a global vertical-structure solution; multiple equilibrium
branches still require care. With SG, this radial sweep participates in the
existing damped SG/thermal iteration. In energy mode initialization uses this
equilibrium, then every substep recalculates geometry at the evolved T.
Ghost temperatures and irradiation follow their adjacent physical cells.
The local relaxation diagnostic holds the resulting irradiation fixed, so it
does **not** diagnose geometry-driven/nonlocal thermal instabilities.

Geometric runs write `irradiation_profile_<time>.dat` alongside each thermal
snapshot, including for HDF5 output:

```text
r_AU surface_AU effective_factor visible_stellar_fraction Tirr_K
```

Visibility describes occultation alone; an unocculted part of the star can
still lie below the local tangent and contribute no direct heating. Headers
and run metadata identify the model and surface/stellar-radius parameters.
The radial sweep incrementally retains the upper convex hull of interior
surface points. Points below an upper chord cannot set the stellar horizon
for any exterior cell, so they can be discarded without a neighbor-distance
cutoff. Each trial surface height queries this hull with binary search.
Hull construction is amortized linear per sweep, and each shadow query is
O(log N), rather than scanning all interior cells. The hull is rebuilt for
every thermal/SG sweep from the newly evaluated interior heights; geometry
is still refreshed at every energy substep. This optimization preserves
distant occulting peaks and does not change the physical shadow model.

This is an opaque prescribed-surface approximation: no scattering, limb
darkening, stellar-wavelength opacity, diffuse heating in shadows, inner-rim
transfer or optical-depth-dependent surface height is solved. Depleted cells
still act as opaque screens, so do not use this approximation to predict
optically thin clearing or cavity illumination. It retains the effective
midplane-temperature closure, rather than separately transporting and
conserving stellar photon energy. Thin disks and `R_star/r << 1` are assumed;
check spatial/time convergence, particularly near shadow boundaries.

**Caveat: thermal optical depth is not stellar-light optical depth.** The
existing opacity and effective optical depth already control thermal emission
in both optically thick and optically thin limits. They do not determine how
much incident stellar light a cell transmits. Geometry mode assumes an opaque
surface at `z_s = chi*H`, without attenuation `exp(-tau_star)` along stellar
rays. Consequently, a depleted inner annulus can overestimate shadowing and
underestimate irradiation of the outer disk. Interpret this model as a
first approximation for disks opaque to stellar light, not as a model of
transparent gaps, cavities or late-stage clearing.

Stellar-light transmission is deferred. Implementing it requires an explicitly
chosen stellar-spectrum-weighted absorption opacity, a consistent density
column along each ray and accounting for absorbed versus transmitted light.
The current local-temperature thermal opacity is not automatically an
appropriate stellar-light opacity; no additional stellar opacity parameter
or optical-depth-dependent irradiation correction is implemented.

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
The softened ring kernel, SG cache update and cache interpolation are implemented
in [self_gravity.c](../src/physics/self_gravity.c), with their public interface in
[self_gravity.h](../include/physics/self_gravity.h). The gas module retains the
orbital-frequency calculation and adds the cached SG contribution when enabled.

The temperature solve uses bounded bisection, relative temperature and
heating/cooling tolerances of 1e-8, and at most 128 iterations. The energy
residual additionally permits a floating-point roundoff allowance of
`64 * DBL_EPSILON * (emitted flux + absorbed flux)`, evaluated with the
effective optical depth. Near irradiation equilibrium the net cooling can
be many orders smaller than either flux; a tolerance relative only to the
small viscous heating can otherwise be unattainable even at adjacent
representable temperatures. This allowance does not relax the temperature
bracket criterion or suppress genuine non-convergence. An already
balanced lower endpoint is also accepted. The lower bound is raised to
`T_rad` (or `T_bg` without irradiation). The upper bound must accommodate
the irradiated temperature, including ghost radii.

Opacity transitions make the thermal residual non-monotone: both bounds can
have net heating while an interior cold equilibrium exists. The solver searches
256 logarithmically spaced intervals from cold to hot for the first sampled
heating-to-cooling crossing, then bisects that local bracket. This selects the
cold restoring branch rather than following the previous temperature; it applies
to prescribed and geometric irradiation, including energy-mode initialization.
It is not a complete root enumeration: pairs of crossings within a sampling
interval can be missed, and a local restoring crossing does not establish
stability of the coupled irradiation/SG disk. Narrower explicit bounds can
isolate a desired restoring branch. A missing sampled crossing or
non-convergence remains an explicit error, never a locally isothermal fallback.

## Opacity and limitations

The log-space opacity implementation uses smoothed Lin-Papaloizou power laws,
cross-checked against the
[FARGO-CPT Lin implementation](https://github.com/rometsch/fargocpt/blob/master/src/opacity.cpp).
It is not Bell-Lin (1994), and does not reuse the historical code's converted
prefactors. It retains the high-temperature electron-scattering branch and
the reference implementation's low-density exception at `rho <= 1e-10 g/cm^3`.
The blends soften transitions; branch switches are not guaranteed exactly
continuous.

The `radiative` mode is instantaneous local equilibrium; `radiative_energy`
advances internal energy as described below. Neither is a full radiative-transfer
solver. Irradiation can use a prescribed height surface, but no self-consistent
vertical absorption structure is solved. There is no radial
radiative diffusion or temperature-dependent molecular weight. The bridge uses
the same opacity for thick/thin terms, rather than separate Rosseland/Planck
means. Its numerical domain (`rho > 0`, `0 < T < 1e9 K`) does not imply physical
validity throughout; constant molecular weight cannot model dissociation or
ionization.

Using radial Omega in `H=cs/Omega` is an effective closure, **not** a solution of
vertical SG. Gas transport retains RAPID's existing Keplerian viscous operator;
generalized heating does not make it a full non-Keplerian angular-momentum
solver. Absorbing ghosts use the positive gas floor only for their thermal
closure, without changing the mass boundary.

## Time-dependent energy mode

```yaml
thermodynamics_parameters:
  model: radiative_energy
  gamma: 1.4
  mean_molecular_weight: 2.34
  background_temperature_k: 10.0
  minimum_temperature_k: 1.0
  maximum_temperature_k: 1500.0
  irradiation_luminosity_lsun: 1.0
  irradiation_factor: 0.05
```

With gas evolution enabled, [energy_equation.c](../src/physics/energy_equation.c)
advances vertically integrated internal energy E = Sigma*u:

```text
u = cv*T,  cv = k_B / [mu*m_H*(gamma - 1)]
P_integrated = (gamma - 1)*E = Sigma*cs^2
dE/dt + (1/r)*d(r*v_r*E)/dr
    = -P_integrated*(1/r)*d(r*v_r)/dr + Qplus - Qminus - u*Sigma_dot_wind
```

P_integrated is distinct from the midplane pressure returned by `calculateEOS`.
The hydrostatic sound speed remains isothermal, cs^2 = k_B*T/(mu*m_H);
it is not multiplied by gamma. Vertical hydrostatic balance remains instantaneous.
E is represented by the stored Sigma and T; neither T nor E is reset to thermal
equilibrium after a step. Initialization alone uses the equilibrium solution.
Thereafter the thermal module evaluates opacity, cooling, viscosity and H at
the evolved temperature, including the coupled SG iteration when enabled.

Only the new mode uses a conservative finite-volume mass update. At face j
between adjacent radial grid points:

```text
F_mass = r_face*Sigma*v_r
       = -3*sqrt(r_face)*d(nu*Sigma*sqrt(r))/dr
F_energy = F_mass*u_upwind
cell_area_per_radian = (r_right_face^2 - r_left_face^2)/2
```

The identical face mass flux transports gas and energy on linear or logarithmic
grids. Pressure work uses the divergence of the face velocity (mass flux divided
by arithmetic-mean face Sigma); the reported cell velocity averages the adjacent
mass fluxes. This remains a Keplerian viscous transport law even with SG, not a
general non-Keplerian angular-momentum solver. Because the old modes retain their
original discretization, comparisons with equilibrium mode require resolution
checks as well as thermal-timescale checks.

The scheme is first-order upwind in energy transport and explicit in time.
Internal substeps limit diffusive transport, outgoing mass, photoevaporation
and thermal/compressional changes. The thermal bound uses viscous heating,
gross radiative emission and absorption separately, not their nearly cancelling
net balance (including near the radiation floor). Excessive
subcycling, invalid states and temperatures outside the configured bounds
terminate with an explicit error; temperatures are not silently clipped.
The existing radiative floor is a net heating term, not an imposed minimum T
in evolved states: expansion can transiently cool below the radiation temperature.

Thermal boundaries use zero-gradient specific internal energy (temperature).
The selected density boundary still controls mass; reflecting boundaries have
exactly zero face flux. Absorbing ghosts carry no mass. Wind sinks remove the
transported local specific energy. If the density floor adds mass, it adds
energy at the transported temperature and emits a warning; floor activation is
not conservative and must be checked in depleted-disk studies.
This is an internal-energy model, not a conservative total-energy solver;
it does not evolve kinetic, orbital or gravitational energy.

Density extrapolation alone does **not** prohibit external mass supply:
zero-gradient and parabolic ghosts can produce inward boundary fluxes.
For an isolated disk in `radiative_energy` mode, select:

```yaml
boundary_conditions:
  inner_boundary_condition: outflow
  outer_boundary_condition: outflow
```

`outflow` uses zero-gradient ghosts but limits the actual boundary face flux
to outward-only (negative at the inner edge, positive at the outer edge).
Energy advection, pressure work and reported velocities use the same limited
mass flux. This prevents boundary mass supply without forcing mass to leave.
It is supported only by `radiative_energy`; other modes reject it explicitly.
Existing boundary choices keep their original behavior. The energy example
configuration selects outflow at both ends; restart the simulation for a changed
boundary configuration to take effect.

Thermal snapshots retain the existing nine columns, so existing temperature
plotters continue to work. Headers and run info identify the energy mode and
gamma. E in CGS can be reconstructed from Sigma and T using cv above.
Qplus and Qminus need not match outside equilibrium.

### Energy diagnostics

Energy-mode runs also write `energy_budget_<time>.dat` and
`thermal_timescales_<time>.dat` in LOGS, with the same cadence, time labels and
actual-time headers as thermal snapshots (also for HDF5 runs).
The original nine-column temperature output is unchanged.

The budget contains one row:

```text
time_yr E_initial E_current heating cooling pressure_work boundary wind_loss floor_added residual relative_residual
```

All energies are in erg, integrated over exact finite-volume annuli and over
all internal substeps since initialization. The discrete balance is
`E_current - E_initial = heating - cooling + pressure_work + boundary
- wind_loss + floor_added + residual`. Boundary and pressure work are signed
energy gains; net radiative cooling can be negative. Wind loss removes local
specific internal energy; floor energy is accounted separately, even if a wind
attempts to remove more gas than the allowed floor. The relative residual divides
by max(abs(E_initial), abs(E_current)). This tests the implemented internal-energy
update, not total (orbital/gravitational) energy conservation.

The timescale file has one row per physical cell:

```text
r_AU E_area t_cool_yr cool_valid t_relax_yr restoring_sign dnet_dT
```

`E_area = Sigma*cv*T` is in erg/cm^2. `t_cool = E_area/Qminus` is reported only
for positive net cooling; otherwise `cool_valid=0` and the stored time is zero
(undefined, not instantaneous cooling). It excludes heating and transport.
`dnet_dT` is a centered temperature derivative of Qminus-Qplus with fractional
temperature offset 1e-4, holding Sigma, Omega, alpha, shear and irradiation fixed
but recalculating hydrostatic H, density, viscosity and opacity.
`t_relax = Sigma*cv/abs(dnet_dT)` is the local linear thermal response estimate:
sign +1 means restoring, -1 means growing/unstable, and 0 means neutral with
an undefined time stored as zero. Away from equilibrium it is a local Jacobian
estimate, not a measured relaxation time or a global stability verdict; it
excludes advection, pressure work and nonlocal SG feedback. Opacity branch
transitions can affect this finite-difference estimate.

## Tests and output

From the project root:

```bash
make test-thermal
make test-energy
python -m unittest discover -s tests -p 'test_energy_plots.py'
python -m unittest discover -s tests -p 'test_thermal_config.py'
# Requires make all; exercises evolution, SG, grids, dust, photoevaporation, I/O.
python -m unittest discover -s tests -p 'test_radiative_simulation.py'
```

The energy tests independently reconstruct face transport, check mass and
internal-energy budgets including pressure work and radiative sources, exercise
reflecting/open boundaries on both grids, verify finite cooling and gamma
dependence, and check first-order timestep convergence. Integration tests cover
irradiation with SG enabled/disabled and reconstruct the evolved EOS pressure.

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

For energy-mode profiles the last panel shows the net thermal source
Qplus - Qminus rather than an equilibrium residual. Negative net cooling
(radiative heating) is supported and uses a signed symlog flux axis.

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
colors; nonnegative fluxes containing zeros use linear colors. When net cooling
is negative, both flux heatmaps share a diverging symlog scale with a linear
threshold of 0.001 times the maximum absolute flux. Temperature and sound-speed changes use a
linear diverging color scale centered on zero (negative for cooling, positive
for increases), fixed at -3% to +3%; values outside this range saturate at the
endpoint colors without changing the data. `--linthresh` sets the
pressure-gradient linear threshold (default 1e-12). White dP/dr=0 contours and
black positive-to-negative crossing markers follow the gas heatmap convention.
At least two paired snapshots are required; times and radial grids are checked.
Initial density, temperature and sound speed mean the earliest snapshot in that run,
normalized separately at each radius. Other thermal panels remain absolute.
