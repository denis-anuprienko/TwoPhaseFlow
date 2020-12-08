#include <cstdio>
#include "header.h"

TwoPhaseFlow::TwoPhaseFlow()
{
    ttt = Timer();
    setDefaultParams();
    std::fill_n(times, 7, 0.0);
}

TwoPhaseFlow::~TwoPhaseFlow()
{
    delete varX;
    delete varPg;
    delete varPhi;
    delete aut;

    printf("+=========================\n");
    printf("| T_assemble = %lf\n", times[T_ASSEMBLE]);
    printf("| T_solve    = %lf\n", times[T_SOLVE]);
    printf("| T_precond  = %lf\n", times[T_PRECOND]);
    printf("| T_liniter  = %lf\n", times[T_LINITER]);
    printf("| T_IO       = %lf\n", times[T_IO]);
    printf("| T_update   = %lf\n", times[T_UPDATE]);
    printf("| T_init     = %lf\n", times[T_INIT]);
    printf("+-------------------------\n");
    printf("| T_total    = %lf\n", Timer() - ttt);
    printf("+=========================\n");
}

void TwoPhaseFlow::setDefaultParams()
{
     K0      = 1e-18;
     mul     = 1.e-3;
     mug     = 1.e-3;
     rhol    = 7.0e2;
     rhog    = 7.0e2;
     g       = 9.81;
     c_f     = 1./22e9;
     c_phi   = 9e-3*1e-6;
     Pt      = 43e6;
     P0      = 1e5;
     gamma   = 0.028*1e-6;
     vg_a    = 1.0;
     vg_n    = 5;
     vg_m    = 1. - 1./vg_n;
     phi0    = 0.16;
}

void TwoPhaseFlow::readParams(std::string path)
{
    double t = Timer();
    std::ifstream file(path);
    std::string line;
    while (getline(file, line))
    {
        std::istringstream iss(line);
        std::string firstword;
        iss >> firstword;
        if(firstword[0] == '#'){
            //std::cout << "Found a comment line" << std::endl;
            continue;
        }

        if(firstword == "K0")
            iss >> K0;
        if(firstword == "c_phi")
            iss >> c_phi;
        if(firstword == "vg_n"){
            iss >> vg_n;
            vg_m = 1.-1./vg_n;
        }
        if(firstword == "vg_a")
            iss >> vg_a;
    }
    times[T_IO] += Timer() - t;
}

void TwoPhaseFlow::readMesh(std::string path)
{
    double t = Timer();
    mesh.Load(path);
    times[T_IO] += Timer()-t;

    t = Timer();
    Mesh::GeomParam param;
    param[ORIENTATION]  = FACE;
    param[MEASURE]      = FACE | CELL;
    param[BARYCENTER]   = FACE | CELL;
    param[NORMAL]       = FACE;
    mesh.PrepareGeometricData(param);
    mesh.AssignGlobalID(CELL|FACE);
    times[T_INIT] += Timer() - t;
}

void TwoPhaseFlow::cleanMesh()
{
    // Tags that are likely to be on mesh
    double t = Timer();
    std::vector<std::string> tagNames;
    tagNames.push_back("Liquid_Saturation");
    tagNames.push_back("Liquid_Saturation_Old");
    tagNames.push_back("Gas_Saturation");
    tagNames.push_back("Gas_Pressure");
    tagNames.push_back("Liquid_Pressure");
    tagNames.push_back("Capillary_Pressure");
    tagNames.push_back("Porosity");
    tagNames.push_back("Porosity_Old");
    tagNames.push_back("Fluid_Pressure");
    tagNames.push_back("Fluid_Pressure_Old");
    tagNames.push_back("Water_Head");
    tagNames.push_back("Water_Saturation");
    tagNames.push_back("Water_Head_Prev");
    tagNames.push_back("Moisture_Content");
    tagNames.push_back("_BLK_0_Offset");
    tagNames.push_back("_BLK_1_Offset");
    tagNames.push_back("_BLK_2_Offset");
    tagNames.push_back("Primary_Variable");
    tagNames.push_back("Primary_Variable_Type");

    for(unsigned i = 0; i < tagNames.size(); i++){
        if(mesh.HaveTag(tagNames[i]))
            mesh.DeleteTag(mesh.GetTag(tagNames[i]));
    }
    times[T_INIT] += Timer() - t;

    //mesh.Save("out.vtk");
}

void TwoPhaseFlow::initTags()
{
    double t = Timer();
    Sl       = mesh.CreateTag("Liquid_Saturation",     DATA_REAL,    CELL, false, 1);
    Pl       = mesh.CreateTag("Liquid_Pressure",       DATA_REAL,    CELL, false, 1);
    Pc       = mesh.CreateTag("Capillary_Pressure",    DATA_REAL,    CELL, false, 1);
    Pg       = mesh.CreateTag("Gas_Pressure",          DATA_REAL,    CELL, false, 1);
    Pf_old   = mesh.CreateTag("Fluid_Pressure_Old",    DATA_REAL,    CELL, false, 1);
    X        = mesh.CreateTag("Primary_Variable",      DATA_REAL,    CELL, false, 1);
    Phi      = mesh.CreateTag("Porosity",              DATA_REAL,    CELL, false, 1);
    Phi_old  = mesh.CreateTag("Porosity_Old",          DATA_REAL,    CELL, false, 1);
    PV       = mesh.CreateTag("Primary_Variable_Type", DATA_INTEGER, CELL, false, 1);
    Sl_old   = mesh.CreateTag("Liquid_Saturation_Old", DATA_REAL,    CELL, false, 1);
    TCoeff   = mesh.CreateTag("TPFA_Coefficient",      DATA_REAL,    FACE, false, 1);
    Sltmp    = mesh.CreateTag("Sl_tmp",                DATA_REAL,    CELL, false, 1);
    Xtmp     = mesh.CreateTag("X_tmp",                 DATA_REAL,    CELL, false, 1);
    Pgtmp    = mesh.CreateTag("Pg_tmp",                DATA_REAL,    CELL, false, 1);
    Phitmp   = mesh.CreateTag("Phi_tmp",               DATA_REAL,    CELL, false, 1);

    // Some tags don't need to be printed
    X.SetPrint(false);
    TCoeff.SetPrint(false);
    Sl_old.SetPrint(false);
    Phi_old.SetPrint(false);
    Sltmp.SetPrint(false);
    Xtmp.SetPrint(false);
    Pgtmp.SetPrint(false);
    Phitmp.SetPrint(false);

    times[T_INIT] += Timer() - t;
}

void TwoPhaseFlow::computeTPFAcoeff()
{
    double t = Timer();
    for(auto iface = mesh.BeginFace(); iface != mesh.EndFace(); iface++){
        if(iface->GetStatus() != Element::Ghost){
            Face face = iface->getAsFace();

            if(face.Boundary()) continue;

            Cell N = face.BackCell();
            Cell P = face.FrontCell();

            double xP[3], xN[3];
            P.Barycenter(xP);
            N.Barycenter(xN);

            double L[3];
            double diffXnorm2 = 0;
            for(unsigned i = 0; i < 3; i++)
                diffXnorm2 += (xP[i] - xN[i])*(xP[i] - xN[i]);
            for(unsigned i = 0; i < 3; i++)
                L[i] = (xP[i]-xN[i]) / diffXnorm2;

            double norF[3];
            face->UnitNormal(norF);

            double coef = L[0]*norF[0] + L[1]*norF[1] + L[2]*norF[2];

            face.Real(TCoeff) = -coef*face.Area();
        }
    }
    times[T_ASSEMBLE] += Timer() - t;
    mesh.ExchangeData(TCoeff, FACE);
}

void TwoPhaseFlow::copyTagReal(Tag Dest, Tag Src, ElementType mask)
{
    if(mask & CELL){
        for(auto icell = mesh.BeginCell(); icell != mesh.EndCell(); icell++){
            if(icell->GetStatus() == Element::Ghost) continue;

            icell->Real(Dest) = icell->Real(Src);
        }
    }
    if(mask & FACE){
        for(auto iface = mesh.BeginFace(); iface != mesh.EndFace(); iface++){
            if(iface->GetStatus() == Element::Ghost) continue;

            iface->Real(Dest) = iface->Real(Src);
        }
    }
    mesh.ExchangeData(Dest, mask);
}

void TwoPhaseFlow::initAutodiff()
{
    aut = new Automatizator;
    Automatizator::MakeCurrent(aut);

    INMOST_DATA_ENUM_TYPE XTagEntryIndex = 0;
    INMOST_DATA_ENUM_TYPE GasPressureTagEntryIndex = 0;
    INMOST_DATA_ENUM_TYPE PorosityTagEntryIndex = 0;

    XTagEntryIndex           = aut->RegisterTag(X, CELL);
    GasPressureTagEntryIndex = aut->RegisterTag(Pg, CELL);
    PorosityTagEntryIndex    = aut->RegisterTag(Phi, CELL);

    varX   = new dynamic_variable(*aut, XTagEntryIndex);
    varPg  = new dynamic_variable(*aut, GasPressureTagEntryIndex);
    varPhi = new dynamic_variable(*aut, PorosityTagEntryIndex);
    aut->EnumerateEntries();

    Residual R("2phase_full", aut->GetFirstIndex(), aut->GetLastIndex());
    R.Clear();
}

void TwoPhaseFlow::runSimulation()
{
    double t = Timer();
    mesh.Save("out.vtk");
    times[T_IO] += Timer() - t;
}

int main(int argc, char *argv[])
{
    if(argc != 3){
        std::cout << "Usage: twophase <param_file_path> <mesh_file_path>" << std::endl;
    }
    Solver::Initialize(&argc,&argv);
    Mesh::Initialize(&argc,&argv);

    TwoPhaseFlow Problem;
    Problem.readMesh(argv[2]);
    Problem.cleanMesh();
    Problem.initTags();
    Problem.computeTPFAcoeff();
    Problem.runSimulation();

    Mesh::Finalize();
    Solver::Finalize();
    return 0;
}
