# TODO

- [ ] Restart a simulation from a selected snapshot, restoring the saved time
  and all gas/dust state needed to continue consistently.
- [ ] Add central-star irradiation to the thermal model, with consistent
  absorption geometry and midplane heating.
- [ ] Add a time-dependent internal-energy equation with gas-consistent energy
  transport, compressional heating/expansion cooling, viscous heating and
  radiative cooling. Retain the instantaneous-equilibrium mode for comparison.
- [ ] Add vertical dust settling, initially using an equilibrium dust scale height.
- [ ] Extend the disk model to 1+1D with a local vertical structure.
- [ ] Make alpha respond to temperature and ionization.
- [ ] Include stellar evolution and couple changing stellar properties to the disk.
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
