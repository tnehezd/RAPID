// src/init_tool_module.h

#ifndef INIT_TOOL_MODULE_H
#define INIT_TOOL_MODULE_H

// Struktúra az init_tool paramétereihez
typedef struct init_tool_options {
    int n;
    double ri;
    double ro;
    double sigma0;
    double sigma0cgs;
    double index;
    double rdze_i;
    double rdze_o;
    double drdze_i;
    double drdze_o;
    double alphaParam;
    double amod;
    double h;
    double flind;
    double m0;
    double md;
    double eps;
    long double ratio;
    long double mic;
    long double onesize;
    // int run_init; // EZT A TAGOT TÖRÖLD! Már nincs rá szükség.
} init_tool_options_t;

// Függvény deklarációk
void create_default_init_tool_options(init_tool_options_t *opt);

int run_init_tool(init_tool_options_t *init_opts);

#endif // INIT_TOOL_MODULE_H