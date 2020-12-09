#include <cstdio>
#include "header.h"

TwoPhaseFlow::TwoPhaseFlow()
    : aut(), times()
{
    ttt = Timer();
    setDefaultParams();
    mesh = &mesh_;
}

TwoPhaseFlow::~TwoPhaseFlow()
{
    if(aut != nullptr)
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
     Sl0     = 0.1;
     Pg0     = 1e6;
     phi0    = 0.16;

     dt      = 1e5;
     T       = 3e5;
     save_dir = ".";
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
        if(firstword == "dt")
            iss >> dt;
        if(firstword == "T")
            iss >> dt;
        if(firstword == "Sl0")
            iss >> Sl0;
        if(firstword == "Pg0")
            iss >> Pg0;
        if(firstword == "phi0")
            iss >> K0;
        if(firstword == "save_dir")
            iss >> save_dir;
        if(firstword == "problem_name")
            iss >> problem_name;
    }
    //std::cout << "Problem name is " << problem_name << std::endl;
    times[T_IO] += Timer() - t;
}

void TwoPhaseFlow::readMesh(std::string path)
{
    double t = Timer();
    mesh->Load(path);
    times[T_IO] += Timer()-t;

    t = Timer();
    Mesh::GeomParam param;
    param[ORIENTATION]  = FACE;
    param[MEASURE]      = FACE | CELL;
    param[BARYCENTER]   = FACE | CELL;
    param[NORMAL]       = FACE;
    mesh->PrepareGeometricData(param);
    mesh->AssignGlobalID(CELL|FACE);
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
        if(mesh->HaveTag(tagNames[i]))
            mesh->DeleteTag(mesh->GetTag(tagNames[i]));
    }
    times[T_INIT] += Timer() - t;

    //mesh->Save("out.vtk");
}

void TwoPhaseFlow::initTags()
{
    double t = Timer();
    Sl       = mesh->CreateTag("Liquid_Saturation",     DATA_REAL,    CELL, false, 1);
    Pl       = mesh->CreateTag("Liquid_Pressure",       DATA_REAL,    CELL, false, 1);
    Pc       = mesh->CreateTag("Capillary_Pressure",    DATA_REAL,    CELL, false, 1);
    Pg       = mesh->CreateTag("Gas_Pressure",          DATA_REAL,    CELL, false, 1);
    Pf_old   = mesh->CreateTag("Fluid_Pressure_Old",    DATA_REAL,    CELL, false, 1);
    X        = mesh->CreateTag("Primary_Variable",      DATA_REAL,    CELL, false, 1);
    Phi      = mesh->CreateTag("Porosity",              DATA_REAL,    CELL, false, 1);
    Phi_old  = mesh->CreateTag("Porosity_Old",          DATA_REAL,    CELL, false, 1);
    PV       = mesh->CreateTag("Primary_Variable_Type", DATA_INTEGER, CELL, false, 1);
    Sl_old   = mesh->CreateTag("Liquid_Saturation_Old", DATA_REAL,    CELL, false, 1);
    TCoeff   = mesh->CreateTag("TPFA_Coefficient",      DATA_REAL,    FACE, false, 1);
    Sltmp    = mesh->CreateTag("Sl_tmp",                DATA_REAL,    CELL, false, 1);
    Xtmp     = mesh->CreateTag("X_tmp",                 DATA_REAL,    CELL, false, 1);
    Pgtmp    = mesh->CreateTag("Pg_tmp",                DATA_REAL,    CELL, false, 1);
    Phitmp   = mesh->CreateTag("Phi_tmp",               DATA_REAL,    CELL, false, 1);

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
    for(auto iface = mesh->BeginFace(); iface != mesh->EndFace(); iface++){
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
    mesh->ExchangeData(TCoeff, FACE);
}

void TwoPhaseFlow::copyTagReal(Tag Dest, Tag Src, ElementType mask)
{
    if(mask & CELL){
        for(auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
            if(icell->GetStatus() == Element::Ghost) continue;

            icell->Real(Dest) = icell->Real(Src);
        }
    }
    if(mask & FACE){
        for(auto iface = mesh->BeginFace(); iface != mesh->EndFace(); iface++){
            if(iface->GetStatus() == Element::Ghost) continue;

            iface->Real(Dest) = iface->Real(Src);
        }
    }
    mesh->ExchangeData(Dest, mask);
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

    varX   = dynamic_variable(*aut, XTagEntryIndex);
    varPg  = dynamic_variable(*aut, GasPressureTagEntryIndex);
    varPhi = dynamic_variable(*aut, PorosityTagEntryIndex);
    aut->EnumerateEntries();

    R = Residual("2phase_full", aut->GetFirstIndex(), aut->GetLastIndex());
    R.Clear();
}

variable TwoPhaseFlow::get_Sl(variable Pcc)
{
    if(Pcc.GetValue() < 0.0)
        return 1.0;
    return pow(1.0 + pow(vg_a/rhol/g*Pcc, vg_n), -vg_m);
}

variable TwoPhaseFlow::get_Pc(variable S)
{
    if(S.GetValue() < 0.0 || S.GetValue() > 1.0){
        mesh->Save("err.vtk");
        std::cout << "Bad saturation " << S.GetValue() << std::endl;
        exit(1);
    }
    return rhol*9.81/vg_a * pow(pow(S,-1./vg_m) - 1., 1./vg_n);
}

void TwoPhaseFlow::assembleResidual()
{
    for(auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
        if(icell->GetStatus() == Element::Ghost) continue;
        Cell cellP = icell->getAsCell();

        // Values needed for accumulation part
        variable SP, Pf, PlP, massL, massG;
        double SPn, massL_old, massG_old;
        SPn = cellP.Real(Sl_old);

        double V = cellP.Volume();

        if(cellP.Integer(PV) == PV_PRES){
            PlP = varX(cellP);
            variable Pcc = varPg(cellP) - PlP;
            SP = get_Sl(Pcc);
        }
        else{
            // Saturation formulation, X is Sl
            SP = varX(cellP);
            PlP = varPg(cellP) - get_Pc(SP);
        }


        massL         = rhol *     SP   * varPhi(cellP);
        massG         = rhog * (1.-SP)  * varPhi(cellP);
        massL_old     = rhol *     SPn  * cellP.Real(Phi_old);
        massG_old     = rhol * (1.-SPn) * cellP.Real(Phi_old);

        Pf            = SP * PlP + (1.-SP) * varPg(cellP);

        R[varX.Index(cellP)]   = (massL - massL_old)/dt;// * V;
        R[varPg.Index(cellP)]  = (massG - massG_old)/dt;// * V;
        //R[varPg.Index(cellP)]  = varPg(cellP) - 1e6;

        R[varPhi.Index(cellP)] = (varPhi(cellP) - cellP.Real(Phi_old))/dt - c_phi*(Pf - cellP.Real(Pf_old))/dt;
        //R[varPhi.Index(cellP)] *= V;

        if(cellP->GetStatus() != Element::Ghost){
            auto faces = cellP->getFaces();
            // Loop over cell faces
            for(auto iface = faces.begin(); iface != faces.end(); iface++){
                Face face = iface->getAsFace();
                if(face.Boundary()){

                }
                else{ // Internal face
                    Cell cellN;
                    if(cellP == face->BackCell())
                        cellN = face->FrontCell();
                    else{
                        cellN = face->BackCell();
                    }
                    if(cellP == cellN){
                        std::cout << "P == N" << std::endl;
                        exit(1);
                    }

                    double coef = face.Real(TCoeff);


                    variable PlN, SN, Krl, Krg, ql, qg;
                    // Liquid pressure for cell N
                    if(cellN.Integer(PV) == PV_PRES){
                        PlN = varX(cellN);
                        variable Pcc = varPg(cellN) - PlN;
                        SN = get_Sl(Pcc);
                    }
                    else{ // X is Sl
                        SN = varX(cellN);
                        PlN = varPg(cellN) - get_Pc(SN);
                    }


                    Krl = 0.5 * (SP*SP + SN*SN);
                    Krg = 0.5 * ((1.-SP)*(1.-SP) + (1.-SN)*(1.-SN));

                    if(Krg.GetValue() < 1e-9)
                        Krg = 1e-9;

                    ql = -rhol*Krl*K0/mul * coef * (PlP - PlN);
                    qg = -rhog*Krg*K0/mug * coef * (varPg(cellP) - varPg(cellN));

                    R[varX.Index(cellP)] += ql/V;
                    R[varPg.Index(cellP)] += qg/V;
                }
            }
        }
    }
}

void TwoPhaseFlow::setInitialConditions()
{
    double t = Timer();
    for(auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
        if(icell->GetStatus() == Element::Ghost) continue;

        icell->Real(Pg) = Pg0;
        icell->Real(Sl) = Sl0;
        if(problem_name == "2phase_center"){
            double x[3];
            icell->Barycenter(x);
            double r = (x[0]-0.006)*(x[0]-0.006) + (x[1]-0.006)*(x[1]-0.006);
            r = sqrt(r);
            if(r < 0.012 / 6.)
                icell->Real(Sl) = 0.95;
        }
        icell->Real(Pl) = icell->Real(Pg) - (get_Pc(icell->Real(Sl))).GetValue();
        icell->Real(Phi) = phi0;
    }
    times[T_INIT] += Timer() - t;
}

void TwoPhaseFlow::makeTimeStep()
{

}

void TwoPhaseFlow::runSimulation()
{
    setInitialConditions();
    double t = Timer();
    mesh->Save(save_dir + "/sol0.vtk");
    times[T_IO] += Timer() - t;

    int nt = static_cast<int>(T/dt);
    for(int it = 1; it <= nt; it++){
        std::cout << "===== TIME STEP " << it << ", T = " << it*dt << " =====" << std::endl;
    }
}

int main(int argc, char *argv[])
{
    if(argc != 3){
        std::cout << "Usage: twophase <param_file_path> <mesh_file_path>" << std::endl;
    }
    Solver::Initialize(&argc,&argv);
    Mesh::Initialize(&argc,&argv);

    TwoPhaseFlow Problem;
    Problem.readParams(argv[1]);
    Problem.readMesh(argv[2]);
    Problem.cleanMesh();
    Problem.initTags();
    Problem.computeTPFAcoeff();
    Problem.runSimulation();

    Mesh::Finalize();
    Solver::Finalize();
    return 0;
}
