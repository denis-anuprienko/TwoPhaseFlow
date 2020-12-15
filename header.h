#ifndef HEADER_H
#define HEADER_H

#include "inmost.h"

using namespace INMOST;

enum{
    T_ASSEMBLE = 0,
    T_SOLVE,
    T_PRECOND,
    T_LINITER,
    T_IO,
    T_INIT,
    T_UPDATE
};

// Primary variable types
enum{
    PV_PRES = 1,
    PV_SAT
};

// BC type
enum{
    BC_DIR = 1,
    BC_NEUM
};

// BC category: liquid or gas
enum{
    BCAT_L = 1,
    BCAT_G
};

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
    double Sl0;    // Initial liquid saturation
    double Sl0_c;  // Initial liquid saturation in center (for 2phase_center)
    double Pg0;    // Initial gas pressure
    double phi0;   // Initial porosity

    double dt, T;
    std::string save_dir;     // results folder
    std::string problem_name; // can be "2phase_center" for a p. with centr. sat. zone
    int maxit;                // Newton iteration limit
    double rtol, atol;        // Newton tolerances
    std::string solver_type;  // linear solver type
    double w;                 // relaxation parameter
    double inflowFluxL;       // liquid flux at top boundary ! DIVIDED BY INFLOW AREA
    double outflowPresL;      // liquid pressure at bottom boundary
    int saveIntensity;        // save VTK every ... step
    bool loadMesh;            // if we load a ready mesh
    std::string meshName;     // mesh file location
    int Nx, Ny, Nz;           // dimensions for Cartesian mesh to generate

    // Mesh
    Mesh *mesh;

    // Tags
    Tag Sl;      // Liquid saturation
    Tag Sl_old;  // Liquid saturation at the previous time step
    Tag Pl;      // Liquid pressure
    Tag Pg;      // Gas pressure
    Tag Pc;      // Capillary pressure
    Tag Pf;      // Phase-averaged fluid pressure
    Tag Pf_old;  // Phase-averaged fluid pressure at the previous time step
    Tag Phi;     // Porosity
    Tag Phi_old; // Porosity at the previous time step
    Tag X;       // Cellwise variable X (either Sl or Pl)
    Tag PV;      // Primary variable type indicator
    Tag TCoeff;  // TPFA coefficient for faces
    Tag Grav;    // TPFA approximated gravity term

    Tag Sltmp;   // Temporary tag for Sl
    Tag Xtmp;    // Temporary tag for Xl
    Tag Pgtmp;   // Temporary tag for Pg
    Tag Phitmp;  // Temporary tag for Phi

    Tag BCtype;  // BC type: Neumann, Dirichlet
    Tag BCval;   // BC value: pressure or flux

    // Autodiff things
    Automatizator *aut;
    dynamic_variable varX;   // variable X (Sl or Pl)
    dynamic_variable varPg;  // variable Pg
    dynamic_variable varPhi; // variable Phi
    Residual R;

    // Solver
    Solver *S;

    // Auxiliary
    double times[7];
    double ttt; // global timer
    int iterLinear;
    int iterNewton;
    double mass;

public:
    TwoPhaseFlow();
    ~TwoPhaseFlow();
    void setDefaultParams();
    void readParams(std::string path);
    void readMesh(std::string path);
    void cleanMesh();
    void initTags();
    void computeTPFAcoeff();
    variable get_Sl(variable Pcc);
    variable get_Pc(variable S);
    void get_Kr(variable S, variable &Krl, variable &Krg);
    void assembleResidual();
    void copyTagReal(Tag Dest, Tag Src, ElementType mask);
    void setMesh();
    void createMesh();
    void setInitialConditions();
    void setBoundaryConditions();
    void setPrimaryVariables();
    void initAutodiff();
    void countMass();
    void makeTimeStep();
    void runSimulation();
};

#endif // HEADER_H
