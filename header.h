#ifndef HEADER_H
#define HEADER_H

#include "inmost.h"

using namespace INMOST;

class TwwoPhase
{
private:
    // Physical parameters
    double K0      = 1e-18;        // Initial intrinsic permeability, m^2
    double mul     = 1.e-3;        // Liquid dynamic viscoity, Pa*s
    double mug     = 1.e-3;        // Gas    dynamic viscoity, Pa*s
    double rhol    = 7.0e2;        // Liquid density, kg/m^3
    double rhog    = 7.0e2;        // Gas    density, kg/m^3
    double rhos    = 2.0e3;        // Solid density, kg/m^3
    double g       = 9.81;         // m/s^2
    double c_f     = 1./22e9;      // parameter from van Noort and Yarushina, 1/Pa
    double c_phi   = 9e-3*1e-6;    // parameter from van Noort and Yarushina, 1/Pa
    double Pt      = 43e6;         // Confining pressure, Pa
    double P0      = 1e5;          // Atmospheric pressure, Pa
    double gamma   = 0.028*1e-6;   // Exponent factor for permeability function
    double vg_a    = 1.0;          // van Genuchten pore parameter
    double vg_n    = 5;            // van Genuchten pore parameter
    double vg_m    = 1. - 1./vg_n;
    double phi0    = 0.16;         // Initial porosity

    // Numerical
    Mesh *mesh;

public:
    TwwoPhase(TwwoPhase)
};

#endif // HEADER_H
