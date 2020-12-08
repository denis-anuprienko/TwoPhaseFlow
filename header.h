#ifndef HEADER_H
#define HEADER_H

#include "inmost.h"

using namespace INMOST;

class TwoPhaseFlow
{
private:
    // Physical parameters
    double K0;     // Initial intrinsic permeability, m^2
    double mul;    // Liquid dynamic viscoity, Pa*s
    double mug;    // Gas    dynamic viscoity, Pa*s
    double rhol;   // Liquid density, kg/m^3
    double rhog;   // Gas    density, kg/m^3
    double rhos;   // Solid density, kg/m^3
    double g;      // m/s^2
    double c_f;    // parameter from van Noort and Yarushina, 1/Pa
    double c_phi;  // parameter from van Noort and Yarushina, 1/Pa
    double Pt;     // Confining pressure, Pa
    double P0;     // Atmospheric pressure, Pa
    double gamma;  // Exponent factor for permeability function
    double vg_a;   // van Genuchten pore parameter
    double vg_n;   // van Genuchten pore parameter
    double vg_m;
    double phi0;   // Initial porosity

    // Numerical
    Mesh *mesh;

public:
    TwoPhaseFlow();
    ~TwoPhaseFlow();
    void readParams(std::string path);
};

#endif // HEADER_H
