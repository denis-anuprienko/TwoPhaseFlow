#include <cstdio>
#include "header.h"

#define V_ID(x, y, z) static_cast<unsigned long long>(((x-0)*(Ny+1)*(Nz+1) + (y-0)*(Nz+1) + (z-0)))

TwoPhaseFlow::TwoPhaseFlow()
    : aut(), times(), iterLinear(0), iterNewton(0), mass(0.0)
{
    ttt = Timer();
    setDefaultParams();
    mesh = new Mesh;
    rank = mesh->GetProcessorRank();
    if(mesh->GetProcessorsNumber() == 1)
        outpExt = ".vtk";
    else
        outpExt = ".pvtk";
}

TwoPhaseFlow::~TwoPhaseFlow()
{
    if(aut != nullptr)
        delete aut;
//    if(mesh != nullptr)
//        delete mesh;

    if(rank == 0){
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
     injectRadius = 1.0;
     outflowPresF = 1e6;
     saveIntensity = 1;
     loadMesh = false;
     Nx = Ny = Nz = 4;
     saveSol = true;
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
        if(firstword == "rhol")
            iss >> rhol;
        if(firstword == "rhog")
            iss >> rhog;
        if(firstword == "mul")
            iss >> mul;
        if(firstword == "mug")
            iss >> mug;
        if(firstword == "Sl0")
            iss >> Sl0;
        if(firstword == "Sl0_c")
            iss >> Sl0_c;
        if(firstword == "Pg0")
            iss >> Pg0;
        if(firstword == "phi0")
            iss >> phi0;
        if(firstword == "Pt")
            iss >> Pt;
        if(firstword == "gamma")
            iss >> gamma;
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
        if(firstword == "injection_radius")
            iss >> injectRadius;
        if(firstword == "fluid_outflow_pressure")
            iss >> outflowPresF;
        if(firstword == "save_intensity")
            iss >> saveIntensity;
        if(firstword == "load_mesh")
            iss >> loadMesh;
        if(firstword == "mesh_name")
            iss >> meshName;
        if(firstword == "Nx")
            iss >> Nx;
        if(firstword == "Ny")
            iss >> Ny;
        if(firstword == "Nz")
            iss >> Nz;
        if(firstword == "save_solution")
            iss >> saveSol;
    }
    //std::cout << "Pdnstr " << outflowPresL << std::endl;
    times[T_IO] += Timer() - t;
}

void TwoPhaseFlow::setMesh()
{
    if(rank == 0){
        if(loadMesh){
            readMesh(meshName);
            cleanMesh();
        }
        else
            createMesh();
        std::cout << "Mesh has " << mesh->NumberOfCells() << " cells" << std::endl;
        std::cout << "Mesh has " << mesh->NumberOfFaces() << " faces" << std::endl;
        std::cout << "Mesh has " << mesh->NumberOfNodes() << " nodes" << std::endl;
    }

    std::cout << "ready to partition mesh\n";

    Partitioner p(mesh);
    p.SetMethod(Partitioner::INNER_KMEANS,Partitioner::Partition);
    p.Evaluate();
    //printf("Proc %d ready to redistr.\n", rank);
    mesh->Redistribute();
    //printf("Proc %d did redistr.\n", rank);
    mesh->ReorderEmpty(CELL|FACE|EDGE|NODE);
    mesh->AssignGlobalID(CELL|FACE|EDGE|NODE);
//    MPI_Barrier(MPI_COMM_WORLD);
    mesh->ExchangeGhost(1, NODE);

    double t = Timer();
    Mesh::GeomParam param;
    param[ORIENTATION]  = FACE;
    param[MEASURE]      = FACE | CELL;
    param[BARYCENTER]   = FACE | CELL;
    param[NORMAL]       = FACE;
    mesh->PrepareGeometricData(param);
    mesh->AssignGlobalID(CELL|FACE);
    times[T_INIT] += Timer() - t;

    std::cout << "Finished partitioning" << std::endl;
}

void TwoPhaseFlow::readMesh(std::string path)
{
    double t = Timer();
    mesh->Load(path);
    times[T_IO] += Timer()-t;
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
    //tagNames.push_back("_BLK_0_Offset");
    //tagNames.push_back("_BLK_1_Offset");
    //tagNames.push_back("_BLK_2_Offset");
    tagNames.push_back("Primary_Variable");
    tagNames.push_back("Primary_Variable_Type");

    for(unsigned i = 0; i < tagNames.size(); i++){
        if(mesh->HaveTag(tagNames[i]))
            mesh->DeleteTag(mesh->GetTag(tagNames[i]));
    }
    times[T_INIT] += Timer() - t;

    //mesh->Save("out.vtk");
}

void TwoPhaseFlow::createMesh()
{
    if(problem_name != "2phase_center"){
        std::cout << "Trying to create mesh for unknown problem" << std::endl;
    }

    double L = 0.012;
    double dx = L/Nx, dy = L/Ny, dz = L/Nz;

    mesh = new Mesh;
    ElementArray<Node> nodes(mesh);
    for(int i = 0; i <= Nx; i++){
        for(int j = 0; j <= Ny; j++){
            for(int k = 0; k <= Nz; k++){
                double coords[3] = {i*dx, j*dy, k*dz};
                nodes.push_back(mesh->CreateNode(coords));
            }
        }
    }

    for(int i = 1; i <= Nx; i++){
        for(int j = 1; j <= Ny; j++){
            for(int k = 1; k <= Nz; k++){
                ElementArray<Node> verts(mesh);
                verts.push_back(nodes[V_ID(i-1, j-1, k-1)]);
                verts.push_back(nodes[V_ID(i-0, j-1, k-1)]);
                verts.push_back(nodes[V_ID(i-1, j-0, k-1)]);
                verts.push_back(nodes[V_ID(i-0, j-0, k-1)]);
                verts.push_back(nodes[V_ID(i-1, j-1, k-0)]);
                verts.push_back(nodes[V_ID(i-0, j-1, k-0)]);
                verts.push_back(nodes[V_ID(i-1, j-0, k-0)]);
                verts.push_back(nodes[V_ID(i-0, j-0, k-0)]);

                const INMOST_DATA_INTEGER_TYPE face_nodes[24] = {0,4,6,2, 1,3,7,5, 0,1,5,4, 2,6,7,3, 0,2,3,1, 4,5,7,6};
                const INMOST_DATA_INTEGER_TYPE num_nodes[6]   = {4,       4,       4,       4,       4,       4};

                mesh->CreateCell(verts,face_nodes,num_nodes,6); // Create the cubic cell in the mesh
            }
        }
    }
    mesh->ResolveShared();

    std::cout << "Created mesh with " << mesh->NumberOfCells() << " cells" << std::endl;
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
    Perm     = mesh->CreateTag("Permeability",          DATA_REAL,    CELL, false, 1);
    Phi      = mesh->CreateTag("Porosity",              DATA_REAL,    CELL, false, 1);
    Phi_old  = mesh->CreateTag("Porosity_Old",          DATA_REAL,    CELL, false, 1);
    PV       = mesh->CreateTag("Primary_Variable_Type", DATA_INTEGER, CELL, false, 1);
    Sl_old   = mesh->CreateTag("Liquid_Saturation_Old", DATA_REAL,    CELL, false, 1);
    TCoeff   = mesh->CreateTag("TPFA_Coefficient",      DATA_REAL,    FACE, false, 1);
    Grav     = mesh->CreateTag("Gravity",               DATA_REAL,    FACE, false, 1);
    Sltmp    = mesh->CreateTag("Sl_tmp",                DATA_REAL,    CELL, false, 1);
    Xtmp     = mesh->CreateTag("X_tmp",                 DATA_REAL,    CELL, false, 1);
    Pgtmp    = mesh->CreateTag("Pg_tmp",                DATA_REAL,    CELL, false, 1);
    Phitmp   = mesh->CreateTag("Phi_tmp",               DATA_REAL,    CELL, false, 1);
    BCtype   = mesh->CreateTag("BCtype",                DATA_INTEGER, FACE, FACE,  3);
    BCval    = mesh->CreateTag("BCval",                 DATA_REAL,    FACE, FACE,  3);

    // Some tags don't need to be printed
    X.SetPrint(false);
    TCoeff.SetPrint(false);
    Grav.SetPrint(false);
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

                if(std::isinf(coef) || std::isnan(coef)){
                    printf("Bad coef at face %d\n", face.LocalID());
                }
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
                if(std::isinf(coef) || std::isnan(coef)){
                    printf("Bad coef at face %d\n", face.LocalID());
                }
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
    //INMOST_DATA_ENUM_TYPE PorosityTagEntryIndex = 0;

    XTagEntryIndex           = aut->RegisterTag(X, CELL);
    GasPressureTagEntryIndex = aut->RegisterTag(Pg, CELL);
    //PorosityTagEntryIndex    = aut->RegisterTag(Phi, CELL);

    varX   = dynamic_variable(*aut, XTagEntryIndex);
    varPg  = dynamic_variable(*aut, GasPressureTagEntryIndex);
    //varPhi = dynamic_variable(*aut, PorosityTagEntryIndex);
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
        mesh->Save("err" + outpExt);
        std::cout << "Bad saturation " << S.GetValue() << std::endl;
        exit(1);
    }
    return rhol*9.81/vg_a * pow(pow(S,-1./vg_m) - 1., 1./vg_n);
}

void TwoPhaseFlow::get_Kr(variable S, variable &Krl, variable &Krg)
{
    double s = S.GetValue();

    if(s > 1.0 || s < 0.0 || std::isinf(s) || std::isnan(s)){
        std::cout << "Bad s = " << s << std::endl;
        exit(1);
    }

    Krl = pow(S, 2.0);
    //Krl = sqrt(S) * pow(1.-pow(1.-pow(S,1./vg_m),vg_m),2.0);

    Krg = pow(1.-S, 2.0);
    //Krg = pow(1.-S,2.0)*(1.-S*S);
}

void TwoPhaseFlow::assembleResidual()
{
    double t = Timer();
    for(auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
        if(icell->GetStatus() == Element::Ghost) continue;
        Cell cellP = icell->getAsCell();

        // Values needed for accumulation part
        variable SP, PfP, PlP, massL, massG, PhiP;
        double SPn, massL_old, massG_old, PhiPn;
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


        PfP           = SP * PlP + (1.-SP) * varPg(cellP);
        PhiPn         = cellP.Real(Phi_old);
        PhiP          = PhiPn + c_phi * (PfP - cellP.Real(Pf_old));

        massL         = rhol *     SP   * PhiP;//varPhi(cellP);
        massG         = rhog * (1.-SP)  * PhiP;//varPhi(cellP);
        massL_old     = rhol *     SPn  * cellP.Real(Phi_old);
        massG_old     = rhog * (1.-SPn) * cellP.Real(Phi_old);

        R[varX.Index(cellP)]   = (massL - massL_old)/dt;// * V;
        R[varPg.Index(cellP)]  = (massG - massG_old)/dt;// * V;
        //R[varPg.Index(cellP)]  = varPg(cellP) - 1e6;

        //R[varPhi.Index(cellP)] = (varPhi(cellP) - cellP.Real(Phi_old))/dt - c_phi*(PfP - cellP.Real(Pf_old))/dt;
        //R[varPhi.Index(cellP)] *= V;

        variable KrlP, KrgP;
        get_Kr(SP, KrlP, KrgP);

        double KP = cellP.Real(Perm);

        auto faces = cellP->getFaces();
        // Loop over cell faces
        for(auto iface = faces.begin(); iface != faces.end(); iface++){
            Face face = iface->getAsFace();
            if(face.Boundary()){
                int faceBCtypeL = face.IntegerArray(BCtype)[BCAT_L];
                int faceBCtypeG = face.IntegerArray(BCtype)[BCAT_G];
                int faceBCtypeF = face.IntegerArray(BCtype)[BCAT_F];
                variable ql = 0.0, qg = 0.0;

                // Liduid
                if(faceBCtypeL == BC_NEUM){
                    //std::cout << "Face with Neumann BC for liquid" << std::endl;
                    ql = face.Area()*face.RealArray(BCval)[BCAT_L];
                }
                else if(faceBCtypeL == BC_DIR){
                    std::cout << "Face with Dirichlet BC for liquid" << std::endl;
                    double PlBC = face.RealArray(BCval)[BCAT_L];
                    variable Krl, Krg;
                    get_Kr(SP, Krl, Krg);
                    variable Ke = exp(-gamma*(Pt - PfP - 0.1e6));

                    double coef = face.Real(TCoeff);

                    ql = -rhol*Krl*KP*Ke/mul * coef * (PlP - PlBC);
                }

                // Gas
//                if(faceBCtypeG == BC_NEUM){
//                    //std::cout << "Face with Neumann BC for gas" << std::endl;
//                    qg = face.Area()*face.RealArray(BCval)[BCAT_G];
//                }
//                else if(faceBCtypeG == BC_DIR){
//                    std::cout << "Face with Dirichlet BC for gas" << std::endl;
//                    double PgBC = face.RealArray(BCval)[BCAT_G];
//                    variable Krl, Krg;
//                    get_Kr(SP, Krl, Krg);
//                    variable Ke = exp(-gamma*(Pt - PfP - 0.1e6));

//                    double coef = face.Real(TCoeff);

//                    qg = -rhog*Krg*K0*Ke/mug * coef * (varPg(cellP) - PgBC);
//                }

                if(faceBCtypeF == BC_NEUM){
                    if(fabs(face.RealArray(BCval)[BCAT_F]) > 1e-15){
                        std::cout << "Face with Neumann BC for fluid" << std::endl;
                        exit(1);
                    }
                }
                else if(faceBCtypeF == BC_DIR){
                    //std::cout << "Face with Dirichlet BC for fluid" << std::endl;
                    double PBC = face.RealArray(BCval)[BCAT_F];

                    variable PlBC, PgBC;
                    // Pf = S*Pl + (1-S)*Pg = Pg + S*(Pl-Pg) = Pg - S*Pc
                    PgBC = PBC + SP*get_Pc(SP);
                    PlBC = (PBC - (1.-SP)*PgBC)/SP;
                    variable Krl, Krg;
                    get_Kr(SP, Krl, Krg);
                    variable Ke = exp(-gamma*(Pt - PfP - 0.1e6));

                    double coef = face.Real(TCoeff);

                    ql = -rhol*Krl*KP*Ke/mul * coef * (PlP - PlBC);
                    qg = -rhog*Krg*KP*Ke/mug * coef * (varPg(cellP) - PgBC);
                }

                R[varX.Index(cellP)] -= ql / V;
                R[varPg.Index(cellP)] -= qg / V;
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


                variable PlN, PfN, SN, Krl, Krg, Ke, ql, qg;

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
                PfN = SN*PlN + (1.-SN)*varPg(cellN);

                Ke = 0.5*(exp(-gamma*(Pt - PfP - 0.1e6))
                        + exp(-gamma*(Pt - PfN - 0.1e6)));

                //Ke = 1.0;

                variable KrlN, KrgN;
                get_Kr(SN, KrlN, KrgN);

                //Krl = 0.5 * (KrlP + KrlN);
                //Krg = 0.5 * (KrgP + KrgN);
                if(PlP.GetValue() > PlN.GetValue())
                    Krl = KrlP;
                else
                    Krl = KrlN;
                if(varPg(cellP).GetValue() > varPg(cellN).GetValue())
                    Krg = KrgP;
                else
                    Krg = KrgN;

                if(Krg.GetValue() < 1e-9)
                    Krg = 1e-9;

                double K, KN = cellN.Real(Perm);

                K = 2.0 * KP*KN / (KP + KN);

                ql = -rhol*Krl*K*Ke/mul * coef * (PlP - PlN);
                qg = -rhog*Krg*K*Ke/mug * coef * (varPg(cellP) - varPg(cellN));

                R[varX.Index(cellP)] += ql/V;
                R[varPg.Index(cellP)] += qg/V;
            }
        }

        //R[varPg.Index(cellP)] = varPg(cellP) - Pg0;
    }
    times[T_ASSEMBLE] += Timer() - t;
}

void TwoPhaseFlow::setInitialConditions()
{
    double t = Timer();
    mass = 0.0;
    srand(time(nullptr));
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
        icell->Real(Pc) = (get_Pc(S)).GetValue();
        icell->Real(Pl) = icell->Real(Pg) - icell->Real(Pc);
        icell->Real(Pf) = S * icell->Real(Pl) + (1.-S) * icell->Real(Pg);
        icell->Real(Phi) = phi0;

        mass += S * icell->Real(Phi) * icell->Volume();
    }

    if(problem_name == "spe"){
        Tag PORO = mesh->GetTag("PORO");
        Tag KK;
//        if(mesh->HaveTag("K"))
//            KK = mesh->GetTag("K");
//        if(mesh->HaveTag("Perm"))
//            KK = mesh->GetTag("K");
        if(mesh->HaveTag("Permeability_scalar"))
            KK = mesh->GetTag("Permeability_scalar");
        else {
            std::cout << "No 'Permeability_scalar' tag!\n";
            exit(-1);
        }

        double x = -200.0, y = 0.0;
        mass = 0.;

        for(auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
            if(icell->GetStatus() == Element::Ghost) continue;

            double c[3];
            icell->Barycenter(c);

            double r = (c[0]-x)*(c[0]-x) + (c[1]-y)*(c[1]-y);
            r = sqrt(r);


            icell->Real(Pg) = Pg0;
            if(r < 400.)
                icell->Real(Sl) = Sl0_c;
            else {
                icell->Real(Sl) = Sl0;
            }

            icell->Real(Perm) = 1e-10*icell->RealArray(KK)[0];
            icell->Real(Phi) = icell->Real(PORO);
            double S = icell->Real(Sl);
            icell->Real(Pl) = icell->Real(Pg) - (get_Pc(S)).GetValue();
            icell->Real(Pc) = icell->Real(Pg) - icell->Real(Pl);
            mass += S * icell->Real(Phi) * icell->Volume();
        }
    }

    mass = rhol * mesh->Integrate(mass);
    mesh->ExchangeData(Perm, CELL);
    mesh->ExchangeData(Pg, CELL);
    mesh->ExchangeData(Pl, CELL);
    mesh->ExchangeData(Pf, CELL);
    mesh->ExchangeData(Sl, CELL);
    mesh->ExchangeData(Phi, CELL);
    times[T_INIT] += Timer() - t;
}

void TwoPhaseFlow::setBoundaryConditions()
{
    if(problem_name != "shale_test")
        return;

    double t = Timer();
    ElementArray<Face> inflowFaces;
    double inflowArea = 0.0;
    for(auto iface = mesh->BeginFace(); iface != mesh->EndFace(); iface++){
        if(iface->GetStatus() == Element::Ghost) continue;

        Face face = iface->getAsFace();
        if(!face.Boundary())
            continue;
        double x[3];
        face.Barycenter(x);

        // Upper boundary - hardcoded
        if(fabs(x[2]-0.012) < 1e-7){
            double r = (x[0]-0.006)*(x[0]-0.006) + (x[1]-0.006)*(x[1]-0.006);
            r = sqrt(r);
            if(r <= injectRadius){
                inflowArea += face.Area();
                inflowFaces.push_back(face);
                inflowCells.push_back(face.BackCell());
            }
        }

        // Lower boundary - hardcoded
        if(fabs(x[2]-0.0) < 1e-7){
            double r = (x[0]-0.006)*(x[0]-0.006) + (x[1]-0.006)*(x[1]-0.006);
            r = sqrt(r);
            if(r <= injectRadius){
                //std::cout << "Bottom boundary face " << face.GlobalID() << ", z = " << x[2] << std::endl;
                face.IntegerArray(BCtype)[BCAT_F] = BC_DIR;
                face.RealArray(BCval)[BCAT_F] = outflowPresF;
            }
        }
    }

    for(auto iface = inflowFaces.begin(); iface != inflowFaces.end(); iface++){
        Face face = iface->getAsFace();
        //std::cout << "Top boundary face " << face.GlobalID() << ", z = " << x[2] << std::endl;
        face.IntegerArray(BCtype)[BCAT_L] = BC_NEUM;
        face.RealArray(BCval)[BCAT_L] = inflowFluxL/inflowArea;
    }
    mesh->ExchangeData(BCtype, FACE);
    mesh->ExchangeData(BCval, FACE);
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

void TwoPhaseFlow::countMass()
{
    double mass_new = 0.0;
    for(auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
        if(icell->GetStatus() == Element::Ghost) continue;

        mass_new += icell->Real(Sl) * icell->Real(Phi) * icell->Volume();
    }
    mass_new = rhol * mesh->Integrate(mass_new);
    std::cout << "Mass change is " << mass_new-mass;
    std::cout << " (" << (mass_new-mass)/mass*1e2 << "%)" << std::endl;
    mass = mass_new;
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

    if(rank == 0){
        std::cout << "Newton: maxit = " << maxit;
        std::cout << ", rtol = " << rtol << ", atol = " << atol << std::endl;
        std::cout << "Solver: " << solver_type << std::endl;
    }
    bool converged = false;
    for(int iter = 0; iter < maxit; iter++){
        setPrimaryVariables();
        assembleResidual();

        r2 = R.Norm();
        if(iter == 0)
            r2_0 = r2;
        if(rank == 0)
            std::cout << " iter " << iter << ", |r|_2 = " << r2 << std::endl;

        if(r2 < atol || r2 < rtol*r2_0){
            converged = true;
            if(rank == 0)
                std::cout << "Converged" << std::endl;
            break;
        }

        t = Timer();
        S->SetMatrix(R.GetJacobian());
        times[T_PRECOND] += Timer() - t;

        //std::cout << "System size: " << R.GetResidual().Size() << std::endl;

        t = Timer();
        bool solved = S->Solve(R.GetResidual(), sol);
        if(!solved){
            std::cout << "Linear solver failed: " << S->GetReason() << std::endl;
            std::cout << "Residual: " << S->Residual() << std::endl;
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
            int gotBad = 0;
            for(auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
                if(icell->GetStatus() == Element::Ghost) continue;

                Cell cell = icell->getAsCell();
                cell.Real(Pg) -= w*sol[varPg.Index(cell)];
                if(cell.Integer(PV) == PV_SAT){
                    double Snew = cell.Real(Sl) - w*sol[varX.Index(cell)];
                    if(Snew > 1.0 || Snew < 0.0){
                        std::cout << "    Bad Sl = " << Snew << " at cell " << cell.GlobalID() << std::endl;
                        gotBad = 1;
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
                //cell.Real(Phi) -= sol[varPhi.Index(cell)];
                double S = cell.Real(Sl);
                cell.Real(Pf) = S*cell.Real(Pl) + (1.-S)*cell.Real(Pg);
                cell.Real(Phi) = cell.Real(Phi_old) + c_phi*(cell.Real(Pf)-cell.Real(Pf_old));
            }
            gotBad = mesh->Integrate(gotBad);
            if(gotBad){
                w *= 0.25;
                if(rank == 0)
                    std::cout << "    Decreasing w to " << w << std::endl;
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
            mesh->Save("err" + outpExt);
            exit(1);
        }
        times[T_UPDATE] += Timer() - t;
        mesh->ExchangeData(Pl, CELL);
        mesh->ExchangeData(Pg, CELL);
        mesh->ExchangeData(Pc, CELL);
        mesh->ExchangeData(Sl, CELL);
        mesh->ExchangeData(Phi, CELL);
        iterNewton++;
    }
    if(!converged){
        std::cout << "Newton failed" << std::endl;
        mesh->Save("err" + outpExt);
        exit(1);
    }
}

void TwoPhaseFlow::runSimulation()
{
    setInitialConditions();
    setBoundaryConditions();
    double t = Timer();
    if(saveSol)
        mesh->Save(save_dir + "/sol0" + outpExt);
    times[T_IO] += Timer() - t;

    std::ofstream out("P.txt");

    initAutodiff();

    S = new Solver(solver_type);
    S->SetParameter("absolute_tolerance","1e-10");
    S->SetParameter("relative_tolerance","1e-6");
    S->SetParameter("maximum_iterations","2000");
    S->SetParameter("gmres_substeps","0");
    S->SetParameter("drop_tolerance","1e-2");

    int nt = static_cast<int>(T/dt);
    for(int it = 1; it <= nt; it++){
        if(rank == 0){
            std::cout << std::endl;
            std::cout << "===== TIME STEP " << it << ", T = " << it*dt << " =====" << std::endl;
        }
        makeTimeStep();
        countMass();

        if(it%saveIntensity == 0){
            t = Timer();
            if(saveSol)
                mesh->Save(save_dir + "/sol" + std::to_string(it/saveIntensity) + outpExt);

            if(inflowCells.size() > 0){
                double avPin = 0.0;
                for(auto icell = inflowCells.begin(); icell != inflowCells.end(); icell++)
                    avPin += icell->Real(Pl);
                out << avPin/inflowCells.size() << std::endl;
            }
            times[T_IO] += Timer() - t;
        }
    }

    delete S;
}

int main(int argc, char *argv[])
{
    if(argc != 2){
        std::cout << "Usage: twophase <param_file_path>" << std::endl;
        exit(1);
    }
    Solver::Initialize(&argc,&argv);
    Mesh::Initialize(&argc,&argv);
    Partitioner::Initialize(&argc,&argv);

    TwoPhaseFlow Problem;
    Problem.readParams(argv[1]);
    Problem.setMesh();
    Problem.initTags();
    Problem.computeTPFAcoeff();

    Problem.runSimulation();

    Solver::Finalize();
    Mesh::Finalize();
    Partitioner::Finalize();
    return 0;
}
