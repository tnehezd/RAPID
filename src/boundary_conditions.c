#include "utils.h"
#include "config.h"
#include <math.h>
#include <stdlib.h>
#include "particle_data.h"
#include "simulation_types.h" 
#include "dust_physics.h"
#include "gas_physics.h"
#include "boundary_conditions.h"


void parabolicExtrapolationToGhostCells(double *input_vector, int reference_index1, int reference_index2, int reference_index3, double *out_coefficient_quadratic, double *out_coefficient_linear, double *out_coefficient_constant, double grid_spacing, const DiskParameters *disk_params) {


/**
 * @brief Parabolic ghost-cell extrapolation for smooth outflow boundaries.
 *
 * This routine computes a second-order polynomial (parabolic) extrapolation
 * based on three interior grid points and evaluates the resulting polynomial
 * at the ghost-cell location. The extrapolated value is then used to fill
 * the ghost cell.
 *
 * This is NOT a physical boundary condition by itself. It does not impose
 * Dirichlet or Neumann constraints. Instead, it provides a smooth numerical
 * continuation of the interior solution beyond the computational domain.
 *
 * In practice, this acts as a non-reflecting, smooth outflow boundary:
 * the field is allowed to "leave" the domain by following its interior trend,
 * without enforcing a fixed value (Dirichlet) or a strict zero gradient
 * (Neumann). This reduces spurious reflections and maintains second-order
 * accuracy at the boundary.
 *
 * Use this method whenever a numerically stable, trend-following outflow
 * closure is desired. For true Dirichlet or Neumann boundaries, use the
 * corresponding explicit ghost-cell formulas instead.
 */


	double x_coordinate_1, x_coordinate_2, x_coordinate_3;
	double value_1, value_2, value_3;	
	double local_coefficient_quadratic, local_coefficient_linear, local_coefficient_constant;

	x_coordinate_1 = disk_params->r_min + (reference_index1-1) * grid_spacing;
	x_coordinate_2 = disk_params->r_min + (reference_index2-1) * grid_spacing;
	x_coordinate_3 = disk_params->r_min + (reference_index3-1) * grid_spacing;
 
	value_1 = input_vector[reference_index1];
	value_2 = input_vector[reference_index2];
	value_3 = input_vector[reference_index3];

	local_coefficient_quadratic = ((value_1 - value_3) / (x_coordinate_1 - x_coordinate_3) - (value_1 - value_2) / (x_coordinate_1 - x_coordinate_2)) / (x_coordinate_3 - x_coordinate_2);
	local_coefficient_linear = (value_1 - value_2) / (x_coordinate_1 - x_coordinate_2) - local_coefficient_quadratic * (x_coordinate_1 + x_coordinate_2);
	local_coefficient_constant = value_1 - local_coefficient_quadratic * x_coordinate_1 * x_coordinate_1 - local_coefficient_linear * x_coordinate_1;

	*out_coefficient_quadratic = local_coefficient_quadratic;
	*out_coefficient_linear = local_coefficient_linear;
	*out_coefficient_constant = local_coefficient_constant;

}


void applyBoundaryConditions(double *input_vector, const DiskParameters *disk_params) {


//  OPEN BOUNDARY: Parabolic extrapolation is applied to both velocity and all other physical quantities.
	

//	parabolicExtrapolationToGhostCells(veinput_vectorc, 1, 2, 3, &a, &b, &c, disk_params->delta_r,disk_params);
//	input_vector[0] =  a * (disk_params->r_min - disk_params->delta_r) * (disk_params->r_min - disk_params->delta_r) + b * (disk_params->r_min - disk_params->delta_r) + c;
	input_vector[0] = input_vector[1];
//	parabolicExtrapolationToGhostCells(input_vector, disk_params->grid_number - 2, disk_params->grid_number - 1, disk_params->grid_number, &a, &b, &c, disk_params->delta_r,disk_params);
//	input_vector[disk_params->grid_number+1] = a * (disk_params->r_max + disk_params->delta_r) * (disk_params->r_max + disk_params->delta_r) + b * (disk_params->r_max + disk_params->delta_r) + c;
	input_vector[disk_params->grid_number+1] = input_vector[disk_params->grid_number];
}
