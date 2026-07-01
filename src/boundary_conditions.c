// FIELD TYPES FOR BC CORRECTION
// 0 = Sigma, radial gas vel.
// 1 = Pressure
// 2 = Pressure Gradient
// (set by caller before applyBoundaryConditions)


#include "boundary_conditions.h"
#include <math.h>


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
    double x1 = dp->r_min + (i1 - 1) * dr;
    double x2 = dp->r_min + (i2 - 1) * dr;
    double x3 = dp->r_min + (i3 - 1) * dr;

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
        case 5: applyLinearExtrapolationInner(v, dp); break;
        case 6: applyLogGridExtrapolationInner(v, dp); break;
        default: applyZeroGradientInner(v, dp); break;
    }

    // OUTER BC
    switch(bc_outer) {
        case 0: applyZeroGradientOuter(v, dp); break;
        case 1: applyParabolicOuter(v, dp); break;
        case 2: applyFixedFluxOuter(v, dp); break;
        case 3: applyAbsorbingOuter(v, dp); break;
        case 4: applyReflectingOuter(v, dp); break;
        case 5: applyLinearExtrapolationOuter(v, dp); break;
        case 6: applyLogGridExtrapolationOuter(v, dp); break;
        default: applyZeroGradientOuter(v, dp); break;
    }
}


/************************************************************
 * INNER BC IMPLEMENTATIONS
 ************************************************************/
void applyZeroGradientInner(double *v, const DiskParameters *dp)
{
	(void)dp; // Unused parameter
    v[0] = v[1];
}

void applyParabolicInner(double *v, const DiskParameters *dp)
{
    double a,b,c;
    parabolicExtrapolationToGhostCells(v, 1, 2, 3, &a, &b, &c, dp->delta_r, dp);
    double xg = dp->r_min - dp->delta_r;
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

void applyLinearExtrapolationInner(double *v, const DiskParameters *dp)
{
    (void)dp; // Unused parameter
    v[0] = 2*v[1] - v[2];
}

void applyLogGridExtrapolationInner(double *v, const DiskParameters *dp)
{
    double dr_left = dp->radial_grid[1] - dp->radial_grid[0];
    double dr_right = dp->radial_grid[2] - dp->radial_grid[1];
    v[0] = v[1] + dr_left/dr_right * (v[1] - v[2]);
}

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
    double xg = dp->r_max + dp->delta_r;
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

void applyLinearExtrapolationOuter(double *v, const DiskParameters *dp)
{
    int N = dp->grid_number;
    v[N+1] = 2*v[N] - v[N-1];
}

void applyLogGridExtrapolationOuter(double *v, const DiskParameters *dp)
{
    int N = dp->grid_number;
    double dr_left  = dp->radial_grid[N]   - dp->radial_grid[N-1];
    double dr_right = dp->radial_grid[N-1] - dp->radial_grid[N-2];
    v[N+1] = v[N] + dr_left/dr_right * (v[N] - v[N-1]);
}
