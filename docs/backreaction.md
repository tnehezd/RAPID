# Stationary dust backreaction (experimental)

Enable with `simulation_parameters.enable_backreaction: true` in YAML or
`-backreaction 1` on the CLI. The default is off and retains the existing
passive-dust solver. See [the example](../configs/config_backreaction.yaml).

The first version supports evolved **lociso gas and one drifting particle
population**, with optional grain growth, on linear or logarithmic grids.
Self-gravity, photoevaporation, a second population, legacy benchmark modes,
and non-CIC smoothing are rejected explicitly. Gas boundaries must be
`zero_gradient` or `reflecting`. Reflecting closes the drag-correction face flux; it does
not reflect Lagrangian grains, which retain the existing absorbing radial
domain treatment. Zero-gradient gas boundaries may supply mass.

## Local mutual drag

This is a steady Keplerian drag balance, not time-dependent two-fluid momentum
evolution. Each grain size contributes separately to the local moments:

```text
epsilon_k = deposited dust surface density / gas surface density
A = sum_k epsilon_k / (1 + St_k^2)
B = sum_k epsilon_k * St_k / (1 + St_k^2)
D = (1 + A)^2 + B^2
P = cs * H * (dP_g/dr) / P_g = -2 * eta * v_K
u_g = ((1 + A) * u_nu - B * P) / D
v_g = (B * u_nu + (1 + A) * P) / (2 * D)
u_d,k = (u_g + 2 * St_k * v_g) / (1 + St_k^2)
```

`u` is radial velocity (positive outward), `v` is the azimuthal offset from
Keplerian rotation, and `u_nu` is the existing dust-free viscous velocity.
The pressure-force convention matches the existing drift model. Stokes
numbers use the local Epstein surface-density prescription.

CIC deposits each representative particle's mass into adjacent physical
gas cells. Each contribution's Stokes number is evaluated using the receiving
cell's gas density; there is no mean-size approximation. Annular cell areas
use midpoint faces, including the two boundary faces. Escaped particles
do not contribute. The drag moments do not discard mass using the diagnostic
dust-density floor. Only the growth/diagnostic dust density uses that floor.
Particle queries interpolate A, B, the viscous velocity, and pressure gradient
before solving the local closure, preserving the passive drift limit.

The loading uses **surface densities**, an unstratified/vertically averaged
approximation. It is not the midplane density enhancement of a settled dust
layer. See [geometry](geometry.md); settling and vertical averaging of the
drag coefficients need a separate model.

## Gas transport and time integration

The original viscous FTCS update is retained. Backreaction adds only the
transport associated with the difference between coupled and dust-free gas
velocities, on the same beginning-of-step state:

```text
delta_u = u_g - u_nu
F_BR,face = 2*pi*r_face * delta_u_face * Sigma_donor
Sigma_new = Sigma_FTCS - dt*(F_BR,right - F_BR,left)/annular_area
```

The correction velocity at a face is the arithmetic average of adjacent
`u_g - u_nu` values. The donor is selected by the correction's direction,
not by the total gas velocity. Internal **correction** fluxes cancel exactly
in the discrete mass budget; reflecting correction fluxes are zero.
The original FTCS boundary treatment and its mass-budget limitations remain.
The complete update is not claimed to conserve annular gas mass exactly:
only the added drag transport has that discrete guarantee. Invalid combined
densities cause an explicit error instead of silent floor mass injection.

The particle RK4, viscous FTCS and drag correction use the frozen beginning-of-step
drag state. Growth, if enabled, is split **after transport**, then the moments
and velocities are refreshed. Gas-outflow and dust-crossing safety factors
are 0.2; growth limits the step to 0.1 of the local growth time. The legacy
viscous timestep bound remains active. The explicit gas diffusion bound is
also capped at `0.2 * min(dr_left, dr_right)^2 / (2 * diffusivity)`, with
`diffusivity = 3*nu*(1 + abs((1+A)/D - 1)) + cs*H*B/D`.
This conservative bound does not assume cancellation between the FTCS
viscosity and its drag correction. It includes pressure feedback even
when viscosity is zero; a small instantaneous flux alone is not a diffusion
stability criterion. Temporal smoothing cannot exceed the instantaneous bounds.

The local closure satisfies equal-and-opposite drag balance and
`u_g + sum(epsilon_k * u_d,k) = u_nu`. The complete disk scheme is nevertheless
operator-split, not a globally momentum-conserving two-fluid integrator:
azimuthal momentum is not an evolved conserved field, grains are Lagrangian
while gas uses FTCS plus an upwind drag correction, and deposition/interpolation introduce
spatial errors. RK4 particle positions do not make the coupled scheme RK4.
No drag heating, streaming instability, or transient stopping-time dynamics
is modeled. Resolve radial pressure structure and check grid/time convergence.

At zero loading, the correction is exactly zero and the original FTCS gas
density update is recovered at identical timesteps, including its existing
density-ghost treatment. Refreshing the drag cache does not change density
ghosts. Adaptive timestep histories may differ, so independent runs need not
be bitwise identical at finite time resolution. Dust splitting and the
feedback growth cap also remain distinct from the legacy passive-dust path.
Initial feedback pressure is reconstructed from EOS rather than the rounded
input pressure file.

Earlier experimental versions replaced FTCS with total-gas upwind fluxes.
Their saved runs are not a same-scheme backreaction control and should not be
mixed with the FTCS-based results. Run metadata now explicitly names
`Original viscous FTCS + conservative drag-only flux correction`.

Run metadata identifies the enabled closure and gas transport. Existing
ASCII and HDF5 snapshots carry the resulting gas/dust states; the initial
`drift_timecale.dat` uses the coupled velocities. Dedicated HDF5 metadata
attributes and full momentum-budget outputs are not implemented here.

## Validation

`make test-backreaction` checks the single-size analytic solution, mixed-size
force balance, dust-free velocity limit, conservative mass deposition and
drag transport with open/closed boundaries on both grids, with fewer/equal/more
particles than gas cells. It also checks exclusion of escaped grains, the
deposited size moments, exact multistep zero-loading agreement with the
original FTCS update at equal timesteps, the correction mass budget relative
to that update, and gas-only spatial convergence to the viscous continuum equation.

`tests/test_backreaction_simulation.py` checks default-off compatibility,
gas response and reduced dust drift under loading, growth and OpenMP
determinism, unsupported-input errors, and the YAML/HDF5 workflow.
