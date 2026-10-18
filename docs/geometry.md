# Disk geometry and flaring

## Locally isothermal disks

In `lociso` (alias `locally_isotherm`), geometry is prescribed:

```text
H/r = h0 * (r / AU)^f
```

`aspect_ratio_at_1au` sets `h0`; `flaring_index` sets `f`. These affect the
sound speed, pressure, viscosity, dust drift, and deadzone transition widths.
With `f = 0`, the prescribed aspect ratio is constant.

## Radiative disks

In both `radiative` and `radiative_energy`, scale height follows the thermal
state through `H = cs / Omega`, with `cs` calculated from temperature.
The resulting `H/r` need not be constant or a single power law.

**Prescribed flaring is disabled in these modes.** A YAML `flaring_index`
or CLI `-flind_init` value is ignored and normalized to zero before profile
generation, including deadzone-width calculation. A nonzero supplied value
produces a warning. The run summary labels flaring as inactive.
CLI option order does not affect this behavior.

`aspect_ratio_at_1au` remains an initialization and deadzone-width parameter:
the initial SG softening guess uses `H = h0 * r`, and deadzone transition
widths use the unflared reference height at their edges. It does **not**
prescribe the radiative disk's evolved aspect ratio. Deadzone widths are
not continuously rescaled with the changing thermal height.

For SG runs, rotation and thermal height are iterated together using softened
radial self-gravity. This is not a self-consistent vertical SG solution.

## Prescribed versus geometric irradiation

- `irradiation_model: prescribed` sets an effective irradiation temperature
  from stellar luminosity and a prescribed factor. It does not adjust the
  illumination to the resulting flaring or shadows.
- `irradiation_model: geometry` uses the surface `z_s = chi * H`, its radial
  slope, the finite stellar disk and occultation by inner disk surfaces.
  Temperature, height and illumination are coupled. Flaring is an outcome,
  not an additional imposed exponent.

The geometric surface is assumed opaque even in depleted rings; thermal
optical depth does not determine stellar-light transmission. See
[thermodynamics](thermodynamics.md) for the solver and irradiation limitations.

## Comparing models

Varying flaring is a physical parameter study in lociso mode. It is not an
independent radiative parameter study: radiative comparisons should inspect
the resulting height, temperature, illumination and shadow profiles instead.
The equilibrium model solves temperature instantaneously; the energy model
evolves it in time. Neither implies an imposed power-law height profile.
