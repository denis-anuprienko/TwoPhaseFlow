#include <cstdio>
#include "header.h"

TwoPhaseFlow::TwoPhaseFlow()
    : aut(), times(), iterLinear(0), iterNewton(0)
{
    ttt = Timer();
    setDefaultParams();
    mesh = &mesh_;
}

TwoPhaseFlow::~TwoPhaseFlow()
{
    if(aut != nullptr)
        delete aut;

    std::cout << "Newton iterations: " << iterNewton << std::endl;
    std::cout << "Linear iterations: " << iterLinear << std::endl;

    printf("\n+=========================\n");
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
     Sl0_c   = 0.9;
     Pg0     = 1e6;
     phi0    = 0.16;

     dt      = 1e5;
     T       = 3e5;
     maxit   = 20;
     rtol    = 1e-6;
     atol    = 1e-9;
     save_dir = ".";
     solver_type = "inner_ilu2";
     w      = 1.0;
     inflowFluxL = 0.0;
     outflowPresL = 1e6;
     saveIntensity = 1;
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
            iss >> T;
        if(firstword == "Sl0")
            iss >> Sl0;
        if(firstword == "Sl0_c")
            iss >> Sl0_c;
        if(firstword == "Pg0")
            iss >> Pg0;
        if(firstword == "phi0")
            iss >> K0;
        if(firstword == "save_dir")
            iss >> save_dir;
        if(firstword == "problem_name")
            iss >> problem_name;
        if(firstword == "maxit")
            iss >> maxit;
        if(firstword == "rtol")
            iss >> rtol;
        if(firstword == "atol")
            iss >> atol;
        if(firstword == "solver_type")
            iss >> solver_type;
        if(firstword == "relax_param")
            iss >> w;
        if(firstword == "liquid_inflow_flux")
            iss >> inflowFluxL;
        if(firstword == "liquid_outflow_pressure")
            iss >> outflowPresL;
        if(firstword == "save_intensity")
            iss >> saveIntensity;
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
    Pf       = mesh->CreateTag("Fluid_Pressure",        DATA_REAL,    CELL, false, 1);
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
    BCtype   = mesh->CreateTag("BCtype",                DATA_INTEGER, FACE, FACE,  2);
    BCval    = mesh->CreateTag("BCval",                 DATA_REAL,    FACE, FACE,  2);

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

            if(face.Boundary()){
                Cell P = face.BackCell();

                double xP[3], xN[3];
                P.Barycenter(xP);
                face.Barycenter(xN);

                double L[3];
                double diffXnorm2 = 0;
                for(unsigned i = 0; i < 3; i++)
                    diffXnorm2 += (xP[i] - xN[i])*(xP[i] - xN[i]);
                for(unsigned i = 0; i < 3; i++)
                    L[i] = (xP[i]-xN[i]) / diffXnorm2;

                double norF[3];
                face.UnitNormal(norF);

                double coef = L[0]*norF[0] + L[1]*norF[1] + L[2]*norF[2];

                face.Real(TCoeff) = -coef*face.Area();
            }
            else{
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
                face.UnitNormal(norF);

                double coef = L[0]*norF[0] + L[1]*norF[1] + L[2]*norF[2];

                face.Real(TCoeff) = -coef*face.Area();
            }
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
    double t = Timer();
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
    times[T_INIT] += Timer() - t;
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
    double t = Timer();
    for(auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
        if(icell->GetStatus() == Element::Ghost) continue;
        Cell cellP = icell->getAsCell();

        // Values needed for accumulation part
        variable SP, PfP, PlP, massL, massG;
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

        PfP            = SP * PlP + (1.-SP) * varPg(cellP);

        R[varX.Index(cellP)]   = (massL - massL_old)/dt;// * V;
        R[varPg.Index(cellP)]  = (massG - massG_old)/dt;// * V;
        //R[varPg.Index(cellP)]  = varPg(cellP) - 1e6;

        R[varPhi.Index(cellP)] = (varPhi(cellP) - cellP.Real(Phi_old))/dt - c_phi*(PfP - cellP.Real(Pf_old))/dt;
        //R[varPhi.Index(cellP)] *= V;

        if(cellP->GetStatus() != Element::Ghost){
            auto faces = cellP->getFaces();
            // Loop over cell faces
            for(auto iface = faces.begin(); iface != faces.end(); iface++){
                Face face = iface->getAsFace();
                if(face.Boundary()){

                    // BC for liquid
                    int faceBCtypeL = face.IntegerArray(BCtype)[BCAT_L];
                    variable ql;
                    if(faceBCtypeL == BC_NEUM){
                        //std::cout << "Face with Neumann BC for liquid" << std::endl;
                        ql = face.RealArray(BCval)[BCAT_L];
                    }
                    else if(faceBCtypeL == BC_DIR){
                        double PlBC = face.RealArray(BCval)[BCAT_L];

                        variable Krl = SP*SP;

                        double coef = face.Real(TCoeff);

                        variable ql = -rhol*Krl*K0/mul * coef * (PlP - PlBC);
                    }
                    R[varX.Index(cellP)] -= ql / V;
                }
                else{ // Internal face
                    Cell cellN;
                    if(cellP == face->BackCell())
                        cellN = face->FrontCell();
                    else{
                        cellN = face->BackCell();
                    }
                    if(cellP == cellN)
                        exit(1);

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
        //R[varPg.Index(cellP)] = varPg(cellP) - Pg0;
    }
    times[T_ASSEMBLE] += Timer() - t;
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
                icell->Real(Sl) = Sl0_c;
        }
        double S = icell->Real(Sl);
        icell->Real(Pl) = icell->Real(Pg) - (get_Pc(S)).GetValue();
        icell->Real(Pf) = S * icell->Real(Pl) + (1.-S) * icell->Real(Pg);
        icell->Real(Phi) = phi0;
    }
    times[T_INIT] += Timer() - t;
}

void TwoPhaseFlow::setBoundaryConditions()
{
//    if(problem_name != "shale_test")
//        return;

    double t = Timer();
    for(auto iface = mesh->BeginFace(); iface != mesh->EndFace(); iface++){
        if(iface->GetStatus() == Element::Ghost) continue;

        Face face = iface->getAsFace();
        if(!face.Boundary())
            continue;
        double x[3];
        face.Barycenter(x);

        // Upper boundary - hardcoded
        if(fabs(x[2]-0.012) < 1e-7){
            //std::cout << "Boundary face " << face.GlobalID() << ", z = " << x[2] << std::endl;
            face.IntegerArray(BCtype)[BCAT_L] = BC_NEUM;
            face.RealArray(BCval)[BCAT_L] = inflowFluxL;
        }

        // Lower boundary - hardcoded
        if(fabs(x[2]-0.0) < 1e-7){
            //std::cout << "Boundary face " << face.GlobalID() << ", z = " << x[2] << std::endl;
            face.IntegerArray(BCtype)[BCAT_L] = BC_DIR;
            face.RealArray(BCval)[BCAT_L] = outflowPresL;
        }
    }
    times[T_INIT] += Timer() - t;
}

void TwoPhaseFlow::setPrimaryVariables()
{
    for(auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
        if(icell->GetStatus() == Element::Ghost) continue;

        // Hardcoded values for now...
        if(icell->Real(Sl) > 0.9){
            icell->Integer(PV) = PV_PRES;
            icell->Real(X) = icell->Real(Pl);
        }
        else{
            icell->Integer(PV) = PV_SAT;
            icell->Real(X) = icell->Real(Sl);
        }
    }
    mesh->ExchangeData(PV, CELL);
    mesh->ExchangeData(X, CELL);
}

void TwoPhaseFlow::makeTimeStep()
{
    // Save old values
    copyTagReal(Sl_old, Sl, CELL);
    copyTagReal(Pf_old, Pf, CELL);
    copyTagReal(Phi_old, Phi, CELL);

    double r2, r2_0;
    double t;

    Sparse::Vector sol("Newton_sol", aut->GetFirstIndex(), aut->GetLastIndex());

    std::cout << "Newton: maxit = " << maxit;
    std::cout << ", rtol = " << rtol << ", atol = " << atol << std::endl;
    std::cout << "Solver: " << solver_type << std::endl;
    bool converged = false;
    for(int iter = 0; iter < maxit; iter++){
        setPrimaryVariables();
        assembleResidual();

        r2 = R.Norm();
        if(iter == 0)
            r2_0 = r2;
        std::cout << " iter " << iter << ", |r|_2 = " << r2 << std::endl;

        if(r2 < atol || r2 < rtol*r2_0){
            converged = true;
            std::cout << "Converged" << std::endl;
            break;
        }

        t = Timer();
        S->SetMatrix(R.GetJacobian());
        times[T_PRECOND] += Timer() - t;

        t = Timer();
        bool solved = S->Solve(R.GetResidual(), sol);
        if(!solved){
            std::cout << "Linear solver failed: " << S->GetReason() << std::endl;
            exit(1);
        }
        iterLinear += S->Iterations();
        times[T_SOLVE] += Timer() - t;

        t = Timer();
        // Line search loop
        copyTagReal(Sltmp, Sl, CELL);
        copyTagReal(Pgtmp, Pg, CELL);
        copyTagReal(Phitmp, Phi, CELL);
        bool lsSuccess = false;
        w = 1.0;
        for(int ils = 0; ils < 7; ils++){
            bool gotBad = false;
            for(auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
                if(icell->GetStatus() == Element::Ghost) continue;

                Cell cell = icell->getAsCell();
                cell.Real(Pg) -= w*sol[varPg.Index(cell)];
                if(cell.Integer(PV) == PV_SAT){
                    double Snew = cell.Real(Sl) - w*sol[varX.Index(cell)];
                    if(Snew > 1.0 || Snew < 0.0){
                        std::cout << "Bad Sl = " << Snew << " at cell " << cell.GlobalID() << std::endl;
                        gotBad = true;
                        break;
                    }

                    cell.Real(Pc) = get_Pc(Snew).GetValue();
                    cell.Real(Pl) = cell.Real(Pg) - cell.Real(Pc);
                    cell.Real(Sl) = Snew;
                }
                else{
                    cell.Real(Pl) -= w*sol[varX.Index(cell)];
                    cell.Real(Pc) = cell.Real(Pg) - cell.Real(Pl);
                    cell.Real(Sl) = get_Sl(cell.Real(Pc)).GetValue();
                }
                cell.Real(Phi) -= sol[varPhi.Index(cell)];
                double S = cell.Real(Sl);
                cell.Real(Pf) = S*cell.Real(Pl) + (1.-S)*cell.Real(Pg);
            }
            if(gotBad){
                w *= 0.5;
                std::cout << "Decreasing w to " << w << std::endl;
                copyTagReal(Sl, Sltmp, CELL);
                copyTagReal(Pg, Pgtmp, CELL);
                copyTagReal(Phi, Phitmp, CELL);
            }
            else {
                lsSuccess = true;
                break;
            }
        }
        if(!lsSuccess){
            std::cout << "Line search failed" << std::endl;
            exit(1);
        }
        times[T_UPDATE] += Timer() - t;
        mesh->ExchangeData(Pl, CELL);
        mesh->ExchangeData(Pg, CELL);
        mesh->ExchangeData(Sl, CELL);
        mesh->ExchangeData(Phi, CELL);
        iterNewton++;
    }
    if(!converged){
        std::cout << "Newton failed" << std::endl;
        exit(1);
    }
}

void TwoPhaseFlow::runSimulation()
{
    setInitialConditions();
    setBoundaryConditions();
    double t = Timer();
    mesh->Save(save_dir + "/sol0.vtk");
    times[T_IO] += Timer() - t;

    std::ofstream out("P.txt");

    initAutodiff();

    S = new Solver(solver_type);
    S->SetParameter("absolute_tolerance","1e-15");
    S->SetParameter("relative_tolerance","1e-10");

    int nt = static_cast<int>(T/dt);
    for(int it = 1; it <= nt; it++){
        std::cout << std::endl;
        std::cout << "===== TIME STEP " << it << ", T = " << it*dt << " =====" << std::endl;
        makeTimeStep();
        out << mesh->CellByLocalID(0).Real(Pl) << std::endl;

        if(it%saveIntensity == 0){
            t = Timer();
            mesh->Save(save_dir + "/sol" + std::to_string(it/saveIntensity) + ".vtk");
            times[T_IO] += Timer() - t;
        }
    }

    delete S;
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
