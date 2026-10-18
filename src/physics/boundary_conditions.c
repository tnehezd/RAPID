// FIELD TYPES FOR BC CORRECTION
// 0 = Sigma, radial gas vel.
// 1 = Pressure
// 2 = Pressure Gradient
// (set by caller before applyBoundaryConditions)


#include "boundary_conditions.h"
#include "logger.h"
#include <math.h>
#include <stdlib.h>


/**
 * @brief Correct BC type for fields where some BCs are invalid.
 *
 * absorbing (3)  → only valid for Sigma
 * fixed_flux (2) → only valid for Sigma
 *
 * For Pressure and dP/dr:
 *   absorbing → fallback to zero-gradient (0)
 *   fixed_flux → fallback to parabolic (1)
 */
int correctBCForField(int bc_type, int field_type)
{
    // ABSORBING BC
    if (bc_type == 3) {
        if (field_type == 0) return 3;   // Sigma OK
        return 0;                        // Pressure, dP/dr → zero-gradient
    }

    // FIXED-FLUX BC
    if (bc_type == 2) {
        if (field_type == 0) return 2;   // Sigma OK
        return 1;                        // Pressure, dP/dr → parabolic
    }

    // All other BCs valid for all fields
    return bc_type;
}


/************************************************************
 * PARABOLIC EXTRAPOLATION (kept exactly as requested)
 ************************************************************/
void parabolicExtrapolationToGhostCells(double *input_vector,
                                        int i1, int i2, int i3,
                                        double *a, double *b, double *c,
                                        double dr,
                                        const DiskParameters *dp)
{
    (void)dr;
    double x1 = dp->radial_grid[i1];
    double x2 = dp->radial_grid[i2];
    double x3 = dp->radial_grid[i3];

    double v1 = input_vector[i1];
    double v2 = input_vector[i2];
    double v3 = input_vector[i3];

    double aq = ((v1 - v3) / (x1 - x3) - (v1 - v2) / (x1 - x2)) / (x3 - x2);
    double bl = (v1 - v2) / (x1 - x2) - aq * (x1 + x2);
    double cn = v1 - aq * x1 * x1 - bl * x1;

    *a = aq;
    *b = bl;
    *c = cn;
}

/************************************************************
 * DISPATCHER
 ************************************************************/
void applyBoundaryConditions(double *v,
                             const DiskParameters *dp,
                             const SimulationOptions *opt)
{
    // Determine corrected BC for this field
    int bc_inner = correctBCForField(opt->inner_boundary_condition_type,
                                     opt->current_bc_target);

    int bc_outer = correctBCForField(opt->outer_boundary_condition_type,
                                     opt->current_bc_target);

    // INNER BC
    switch(bc_inner) {
        case 0: applyZeroGradientInner(v, dp); break;
        case 1: applyParabolicInner(v, dp); break;
        case 2: applyFixedFluxInner(v, dp); break;
        case 3: applyAbsorbingInner(v, dp); break;
        case 4: applyReflectingInner(v, dp); break;
//        case 5: applyFreeOutflowInner(v, dp); break;
        default: applyZeroGradientInner(v, dp); break;
    }

    // OUTER BC
    switch(bc_outer) {
        case 0: applyZeroGradientOuter(v, dp); break;
        case 1: applyParabolicOuter(v, dp); break;
        case 2: applyFixedFluxOuter(v, dp); break;
        case 3: applyAbsorbingOuter(v, dp); break;
        case 4: applyReflectingOuter(v, dp); break;
        default: applyZeroGradientOuter(v, dp); break;
    }
}


/************************************************************
 * INNER BC IMPLEMENTATIONS
 ************************************************************/
void applySurfaceDensityBoundaryConditions(double *v,
                                          const DiskParameters *dp,
                                          const SimulationOptions *opt)
{
    SimulationOptions density_options = *opt;
    density_options.current_bc_target = 0;
    applyBoundaryConditions(v, dp, &density_options);

    static int warned[2] = {0, 0};
    const int indices[2] = {0, dp->grid_number + 1};
    const int types[2] = {opt->inner_boundary_condition_type,
                          opt->outer_boundary_condition_type};
    for (int side = 0; side < 2; ++side) {
        if (types[side] != 1) continue;
        int i = indices[side];
        if (!isfinite(v[i]) || !isfinite(dp->density_floor) || dp->density_floor <= 0.0) {
            LOG_ERROR("Invalid parabolic density boundary or density floor at ghost cell %d.", i);
            exit(EXIT_FAILURE);
        }
        if (v[i] < dp->density_floor) {
            int adjacent = side == 0 ? 1 : dp->grid_number;
            if (!isfinite(v[adjacent]) || v[adjacent] < dp->density_floor) {
                LOG_ERROR("Invalid physical density %.9g next to ghost cell %d.",
                          v[adjacent], i);
                exit(EXIT_FAILURE);
            }
            if (!warned[side]) {
                LOG_WARN("%s parabolic density ghost %.9g undershoots floor %.9g at r=%.9g AU; "
                         "using adjacent density %.9g (zero gradient). "
                         "subsequent activations on this side will not repeat this warning.",
                         side == 0 ? "Inner" : "Outer", v[i], dp->density_floor,
                         dp->radial_grid[i], v[adjacent]);
                warned[side] = 1;
            }
            v[i] = v[adjacent];
        }
    }
}

void applyZeroGradientInner(double *v, const DiskParameters *dp)
{
	(void)dp; // Unused parameter
    v[0] = v[1];
}

void applyParabolicInner(double *v, const DiskParameters *dp)
{
    double a,b,c;
    parabolicExtrapolationToGhostCells(v, 1, 2, 3, &a, &b, &c, dp->delta_r, dp);
    double xg = dp->radial_grid[0];
    v[0] = a*xg*xg + b*xg + c;
}

void applyFixedFluxInner(double *v, const DiskParameters *dp)
{
    double r0 = dp->radial_grid[0];
    double r1 = dp->radial_grid[1];
    v[0] = v[1] * sqrt(r1 / r0);
}

void applyAbsorbingInner(double *v, const DiskParameters *dp)
{
    (void)dp; // Unused parameter
    v[0] = 0.0;
}

void applyReflectingInner(double *v, const DiskParameters *dp)
{
    (void)dp; // Unused parameter
    v[0] = v[2];
}

/*void applyFreeOutflowInner(double *v, const DiskParameters *dp)
{
    (void)dp;
    if (v[1] < v[2]) {
        v[0] = v[1];
    } else {
        v[0] = 2.0 * v[1] - v[2];
    }
}*/


/************************************************************
 * OUTER BC IMPLEMENTATIONS
 ************************************************************/
void applyZeroGradientOuter(double *v, const DiskParameters *dp)
{
    int N = dp->grid_number;
    v[N+1] = v[N];
}

void applyParabolicOuter(double *v, const DiskParameters *dp)
{
    int N = dp->grid_number;
    double a,b,c;
    parabolicExtrapolationToGhostCells(v, N-2, N-1, N, &a, &b, &c, dp->delta_r, dp);
    double xg = dp->radial_grid[N + 1];
    v[N+1] = a*xg*xg + b*xg + c;
}

void applyFixedFluxOuter(double *v, const DiskParameters *dp)
{
    int N = dp->grid_number;
    v[N+1] = v[N];  // simple outflow
}

void applyAbsorbingOuter(double *v, const DiskParameters *dp)
{
    int N = dp->grid_number;
    v[N+1] = 0.0;
}

void applyReflectingOuter(double *v, const DiskParameters *dp)
{
    int N = dp->grid_number;
    v[N+1] = v[N-1];
}
