#ifndef PHOTOEVAPORATION_COEFFS_H
#define PHOTOEVAPORATION_COEFFS_H

#define WIND_CONV   1.0
#define WIND_CRIT   1e-20
#define L_X_DEFAULT 1e30   /* erg/s */

/* Owen et al. 2012 coefficients */
static const double OWEN_A1 = 0.15138;
static const double OWEN_B1 = -1.2182;
static const double OWEN_C1 = 3.4046;
static const double OWEN_D1 = -3.5717;
static const double OWEN_E1 = -0.32762;
static const double OWEN_F1 = 3.6064;
static const double OWEN_G1 = -2.4918;

static const double OWEN_A2 = -0.438226;
static const double OWEN_B2 = -0.10658387;
static const double OWEN_C2 = 0.5699464;
static const double OWEN_D2 = 0.010732277;
static const double OWEN_E2 = -0.131809597;
static const double OWEN_F2 = -1.32285709;

/* Picogna et al. 2019 coefficients */
static const double PICO_A  = -0.5885;
static const double PICO_B  = 4.313;
static const double PICO_C  = -12.1214;
static const double PICO_D  = 16.3587;
static const double PICO_E  = -11.4721;
static const double PICO_F  = 5.7248;
static const double PICO_G  = -2.8562;

/* Picogna hole-case */
static const double PICO_AA = 0.11843;
static const double PICO_BB = 0.99695;
static const double PICO_CC = 0.48835;

/* Picogna Lx-fit */
static const double PICO_A_L = -2.7326;
static const double PICO_B_L = 3.3307;
static const double PICO_C_L = -0.0029868;
static const double PICO_D_L = -7.258;

#endif
