# TODO

- [ ] Restart a simulation from a selected snapshot, restoring the saved time
  and all gas/dust state needed to continue consistently.
- [ ] Add central-star irradiation to the thermal model, with consistent
  absorption geometry and midplane heating.
  - [ ] Stellar irradiation is connected to mass? 
  - [ ] Connect irradiation angle to disk geometry and stellar position.
  - [ ] Self-consistent absorption/vertical geometry, including shadowing and
    reduced absorption in optically thin regions.
- [ ] Add vertical dust settling, initially using an equilibrium dust scale height.
- [ ] Extend the disk model to 1+1D with a local vertical structure.
- [ ] Make alpha respond to temperature and ionization.
- [ ] Include stellar evolution and couple changing stellar properties to the disk.
  - [ ] Measure inward accretion from the same discrete inner-boundary mass flux
    used by the gas evolution; the continuum diagnostic is
    `Mdot_inner = -2*pi*r*Sigma*v_r` (positive inward).
  - [ ] Model delayed stellar accretion through an unresolved inner-disk reservoir.
    Gas crossing a grid boundary near 1 AU still carries orbital angular
    momentum; it does not immediately reach the star by free fall.
    - Start with `dM_res/dt = Mdot_inner - M_res/t_drain`,
      `Mdot_star = (1-f_wind)*M_res/t_drain` and
      `Mdot_wind = f_wind*M_res/t_drain`, with `0 <= f_wind <= 1`.
      Define the reservoir's initial mass and treatment of any outward boundary
      flow explicitly; track reservoir, stellar and expelled mass consistently.
    - Estimate the unresolved transport timescale from
      `t_nu ~ r^2/nu = 1/(alpha*h^2*Omega)`, where `h = H/r`.
      At 1 AU around 1 M_sun, alpha=0.003 and h=0.03 give roughly 59,000 yr.
      This is an order-of-magnitude viscous timescale, not an exact travel time.
      Reservoir drainage gives a distributed response, not a fixed time shift.
    - Reserve a near-free-fall delay for the final magnetospheric flow to the
      stellar surface. Radial free fall from rest at 1 AU takes about 65 days
      for 1 M_sun, but neglects the orbital angular momentum of disk gas.
    - Track angular momentum as well as mass, including transport across the
      unresolved region and any wind losses. Do not count photoevaporated gas,
      outer-boundary losses or planetary accretion as stellar accretion.
    - Couple `dM_star/dt = Mdot_star` back to the gravitational model and derive
      `L_acc ~ G*M_star*Mdot_star/R_star * (1-R_star/R_mag)` from the delayed
      stellar accretion rate, not directly from the flux at 1 AU. Here R_mag is
      the magnetospheric truncation radius, distinct from the grid boundary.
      Assess how the accretion luminosity contributes to irradiation.
    - Later replace the reservoir with a resolved inner-disk model if needed.
- [ ] Add planet-induced gaps and evolving planet mass.
- [ ] Add pebble accretion.
- [ ] Add a two-fluid gas/dust model.
- [ ] Add dust backreaction on the gas through mutual drag, with consistent
  momentum exchange between gas and dust.
- [ ] Extend collisional dust evolution with a Smoluchowski solver and/or
  Monte Carlo collisions, including coagulation and fragmentation.
- [ ] Explore a finite-width inner accretion sink with a stellar-distance-dependent
  removal timescale instead of imposing an abrupt zero-density edge. Account for
  removed mass and angular momentum and check sensitivity to sink parameters.
- [ ] Explore an outer sink modulated by binary orbital phase as a phenomenological
  boundary model. Specify circumstellar versus circumbinary geometry and distinguish
  mass removal from companion-driven tidal torques.
- [ ] Explore episodic accretion and FUor/EXor-like outbursts: start with prescribed
  accretion-luminosity events coupled to irradiation, then assess self-consistent
  disk-instability models. Treat burst amplitudes, durations and recurrence as
  model-dependent rather than universal stellar properties.
- [ ] Add synthetic spectral energy distributions (SEDs), initially as disk
  post-processing using dust temperatures and wavelength-dependent dust opacities,
  including stellar emission, viewing inclination and distance.
- [ ] Track species-dependent ice/snow lines (e.g. H2O, CO2 and CO), initially as
  temperature-based diagnostics, then using pressure-dependent sublimation and
  condensation. Couple ice fractions to dust fragmentation properties and opacity;
  assess vapor/ice transport and finite condensation times during outbursts.
- [ ] Update the documentation to match the implemented physics, configuration,
  outputs and limitations.
- [ ] Update rapidsim to support the new simulation capabilities.
- [ ] Particle data into one consolidated structure? 
- [ ] Move self-gravity session to a different module.