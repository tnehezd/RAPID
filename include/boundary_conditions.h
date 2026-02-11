/**
 * @file boundary_conditions.h
 * @brief Boundary condition handlers and ghost-cell extrapolation.
 */

#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include "simulation_types.h" 

/**
 * @brief Parabolic ghost-cell extrapolation for smooth outflow boundaries.
 * * Fits a second-order polynomial \f$y(r) = a \cdot r^2 + b\cdot r + c\f$ through three interior points
 * to provide a smooth numerical continuation beyond the computational domain.
 *
 * @param[in]  input_vector             The physical field to be extrapolated.
 * @param[in]  reference_index1         First interior grid index for the fit.
 * @param[in]  reference_index2         Second interior grid index for the fit.
 * @param[in]  reference_index3         Third interior grid index for the fit.
 * @param[out] out_coefficient_quadratic Pointer to store the quadratic coefficient (a).
 * @param[out] out_coefficient_linear    Pointer to store the linear coefficient (b).
 * @param[out] out_coefficient_constant  Pointer to store the constant coefficient (c).
 * @param[in]  grid_spacing             Radial grid resolution (dr).
 * @param[in]  disk_params              Pointer to global disk parameters.
 * @return void
 */
void parabolicExtrapolationToGhostCells(double *input_vector, int reference_index1, int reference_index2, int reference_index3, 
	                                    double *out_coefficient_quadratic, double *out_coefficient_linear, double *out_coefficient_constant, 
	                                    double grid_spacing, const DiskParameters *disk_params);

/**
 * @brief Applies physical and numerical boundary conditions to a given field.
 * * Updates ghost cells at the inner and outer boundaries to ensure stability.
 *
 * @param[in,out] input_vector          The array to be updated with boundary values.
 * @param[in]     disk_params           Pointer to global disk parameters.
 * @return void
 */
void applyBoundaryConditions(double *input_vector,const DiskParameters *disk_params);

#endif // BOUNDARY_CONDITIONS_H
