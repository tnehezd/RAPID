#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include "simulation_types.h"

/*
 * ============================================================
 *  BOUNDARY CONDITIONS MODULE — FULL DESCRIPTION
 * ============================================================
 *
 *  PURPOSE:
 *    This module implements boundary conditions for the
 *    EULERIAN GAS FIELDS of the 1D viscous disk model:
 *
 *        - gas surface density Σ(r)
 *        - gas pressure P(r)
 *        - gas pressure gradient dP/dr
 *        - gas radial velocity v_gas(r)
 *
 *    These quantities live on the Euler grid and require
 *    ghost-cell boundary conditions for numerical stability.
 *
 *
 *  IMPORTANT:
 *    DUST IS NOT TREATED AS AN EULERIAN FIELD.
 *
 *    Dust is represented by LAGRANGIAN PARTICLES.
 *    Therefore:
 *
 *        *** NO boundary conditions are applied to dust here. ***
 *
 *    Dust “boundary behaviour” (infall, escape, reflection)
 *    is handled entirely inside the particle integrator.
 *
 *
 * ============================================================
 *  FIELD-SPECIFIC BC LOGIC
 * ============================================================
 *
 *  Some BC types are only physically meaningful for Σ(r).
 *  For example:
 *
 *    - Absorbing BC (3)      → only valid for Σ (dust/gas sink)
 *    - Fixed-flux BC (2)     → only valid for Σ (viscous flux)
 *
 *  Pressure and pressure-gradient cannot use these BCs.
 *
 *  Therefore the dispatcher automatically corrects BC types:
 *
 *    FIELD TYPE:
 *       0 = Σ(r)
 *       1 = P(r)
 *       2 = dP/dr
 *
 *    CORRECTION RULES:
 *       absorbing (3):
 *           Σ → absorbing
 *           P → zero-gradient
 *           dP/dr → zero-gradient
 *
 *       fixed_flux (2):
 *           Σ → fixed_flux
 *           P → parabolic
 *           dP/dr → parabolic
 *
 *       all other BCs:
 *           valid for all fields
 *
 *
 * ============================================================
 *  INNER BOUNDARY CONDITION TYPES
 * ============================================================
 *
 *   0 = Zero-gradient inner
 *       Neumann BC. Ghost cell = first interior cell.
 *       Simple, stable, non-physical.
 *
 *   1 = Parabolic inner
 *       Smooth extrapolation using 3 interior points.
 *       Non-reflecting, good for gas outflow-like inner edges.
 *
 *   2 = Fixed-flux inner (Lynden–Bell)
 *       Enforces Σ r^(1/2) = const.
 *       Physical viscous disk inner boundary.
 *
 *   3 = Absorbing inner
 *       Ghost cell = 0. Dust/gas sink at inner edge.
 *       Physical for dust; gas fallback handled by correction.
 *
 *   4 = Reflecting inner
 *       Ghost cell mirrors interior. Non-physical; test only.
 *
 *   5 = Linear extrapolation inner
 *       v0 = 2*v1 - v2. Simple, stable continuation.
 *
 *   6 = Log-grid extrapolation inner
 *       Extrapolation respecting logarithmic spacing.
 *
 *
 * ============================================================
 *  OUTER BOUNDARY CONDITION TYPES
 * ============================================================
 *
 *   0 = Zero-gradient outer
 *       Neumann BC. Simple outflow.
 *
 *   1 = Parabolic outer
 *       Smooth, non-reflecting outflow.
 *       Excellent for viscous diffusion stability.
 *
 *   2 = Fixed-flux outer
 *       Viscous disk outer-edge flux condition.
 *
 *   3 = Absorbing outer
 *       Ghost cell = 0. Dust sink at outer boundary.
 *
 *   4 = Reflecting outer
 *       Ghost cell mirrors interior. Non-physical.
 *
 *   5 = Linear extrapolation outer
 *       v(N+1) = 2*v(N) - v(N-1). Very stable outflow.
 *
 *   6 = Log-grid extrapolation outer
 *       Extrapolation respecting logarithmic spacing.
 *
 *
 * ============================================================
 *  NOTE ON RADIAL GRID:
 * ============================================================
 *
 *   The radial grid r[i] is GEOMETRY, not a physical field.
 *   Boundary conditions DO NOT apply to the radial grid.
 *
 *   Ghost-cell coordinates are set by createRadialGrid()
 *   and must not be modified by BC routines.
 *
 * ============================================================
 */



/* Dispatcher */
void applyBoundaryConditions(double *v,
                             const DiskParameters *dp,
                             const SimulationOptions *opt);

/* Parabolic helper */
void parabolicExtrapolationToGhostCells(double *input_vector,
                                        int i1, int i2, int i3,
                                        double *a, double *b, double *c,
                                        double grid_spacing,
                                        const DiskParameters *dp);

/* INNER BCs */
void applyZeroGradientInner(double *v, const DiskParameters *dp);
void applyParabolicInner(double *v, const DiskParameters *dp);
void applyFixedFluxInner(double *v, const DiskParameters *dp);
void applyAbsorbingInner(double *v, const DiskParameters *dp);
void applyReflectingInner(double *v, const DiskParameters *dp);
void applyLinearExtrapolationInner(double *v, const DiskParameters *dp);
void applyLogGridExtrapolationInner(double *v, const DiskParameters *dp);

/* OUTER BCs */
void applyZeroGradientOuter(double *v, const DiskParameters *dp);
void applyParabolicOuter(double *v, const DiskParameters *dp);
void applyFixedFluxOuter(double *v, const DiskParameters *dp);
void applyAbsorbingOuter(double *v, const DiskParameters *dp);
void applyReflectingOuter(double *v, const DiskParameters *dp);
void applyLinearExtrapolationOuter(double *v, const DiskParameters *dp);
void applyLogGridExtrapolationOuter(double *v, const DiskParameters *dp);

#endif
