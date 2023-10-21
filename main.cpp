#include <cstdio>
#include "header.h"

//#define V_ID(x, y, z) static_cast<unsigned long long>(((x-0)*(Ny+1)*(Nz+1) + (y-0)*(Nz+1) + (z-0)))
#define V_ID(x, y, z) ((x-localstart[0])*(localsize[1]+1)*(localsize[2]+1) + (y-localstart[1])*(localsize[2]+1) + (z-localstart[2]))

#define Rsample 0.0125
#define Lsample 0.01

double inflowArea;
double dtDesir, dtnew;

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
        outpExt = ".pvtu";
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
        printf("| T_mesh     = %lf\n", times[T_MESHGEN]);
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
     rhol0   = 7.0e2;
     rhog    = 7.0e2;
     g       = 9.81;
     c_f     = 1./22e9;
     c_phi   = 9e-3*1e-6;
     c_w     = 1e-6;
     Pt      = 15e6;
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
        dtDesir = dt;
        if(firstword == "T")
            iss >> T;
        if(firstword == "rhol")
            iss >> rhol0;
        if(firstword == "c_w")
            iss >> c_w;
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

    if(rank == 0)
        printf("mu/k/rho RATIO = %e\n", mul/K0/rhol0);
}

void TwoPhaseFlow::setMesh()
{
    if(loadMesh){
        if(rank == 0){
            readMesh(meshName);
            cleanMesh();
        }
    }
    else
        createMesh(); // always in parallel mode!

    int nc = mesh->Integrate(mesh->NumberOfCells());
    int nf = mesh->Integrate(mesh->NumberOfFaces());
    int nn = mesh->Integrate(mesh->NumberOfNodes());
    if(rank == 0){
        std::cout << "Mesh has " << nc << " cells" << std::endl;
        std::cout << "Mesh has " << nf << " faces" << std::endl;
        std::cout << "Mesh has " << nn << " nodes" << std::endl;
    }


    //std::cout << "ready to partition mesh\n";

    Partitioner p(mesh);
    p.SetMethod(Partitioner::INNER_KMEANS, Partitioner::Repartition);
    p.Evaluate();
    //printf("Proc %d ready to redistr.\n", rank);
    mesh->Redistribute();
    //printf("Proc %d did redistr.\n", rank);
    mesh->ReorderEmpty(CELL|FACE|EDGE|NODE);
    mesh->AssignGlobalID(CELL|FACE|EDGE|NODE);
//    MPI_Barrier(MPI_COMM_WORLD);
    mesh->ExchangeGhost(1, NODE);

    std::cout << "Process " << rank << " has " << mesh->NumberOfCells() << " cells" << std::endl;

    double t = Timer();
    Mesh::GeomParam param;
    param[ORIENTATION]  = FACE;
    param[MEASURE]      = FACE | CELL;
    param[BARYCENTER]   = FACE | CELL;
    param[NORMAL]       = FACE;
    mesh->PrepareGeometricData(param);
    mesh->AssignGlobalID(CELL|FACE);
    times[T_INIT] += Timer() - t;

    if(rank == 0)
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
    double t = Timer();
    if(problem_name != "2phase_center"){
        if(rank == 0)
            std::cout << "Trying to create mesh for unknown problem" << std::endl;
    }

    double dx = 2*Rsample/Nx, dy = 2*Rsample/Ny, dz = Lsample/Nz;

    mesh = new Mesh;

    int procs_per_axis[3] = {1,1,1};
    int sizes[3] = {Nx,Ny,Nz};

    {
        int divsize = mesh->GetProcessorsNumber();
        std::vector<int> divs;
        while( divsize > 1 )
        {
            for(int k = 2; k <= divsize; k++)
                if( divsize % k == 0 )
                {
                    divs.push_back(k);
                    divsize /= k;
                    break;
                }
        }
        int elements_per_procs[3] = {sizes[0],sizes[1],sizes[2]};
        for(std::vector<int>::reverse_iterator it = divs.rbegin(); it != divs.rend(); it++)
        {
            int * max = std::max_element(elements_per_procs+0,elements_per_procs+3);
            procs_per_axis[max-elements_per_procs] *= *it;
            (*max) /= *it;
        }
    }

    //rank = proc_coords[2] * procs_per_axis[0] *procs_per_axis[1] + proc_coords[1] * procs_per_axis[0] + proc_coords[0];
    int proc_coords[3] = {rank % procs_per_axis[0] , rank / procs_per_axis[0] % procs_per_axis[1], rank / (procs_per_axis[0] *procs_per_axis[1]) };

    int localsize[3], localstart[3], localend[3];
    int avgsize[3] =
            {
                    (int)ceil((double)sizes[0]/procs_per_axis[0]),
                    (int)ceil((double)sizes[1]/procs_per_axis[1]),
                    (int)ceil((double)sizes[2]/procs_per_axis[2])
            };

    for(int j = 0; j < 3; j++)
    {
        localstart[j] = avgsize[j] * proc_coords[j];
        if( proc_coords[j] == procs_per_axis[j] - 1 )
            localsize[j] = sizes[j] - avgsize[j] * (procs_per_axis[j]-1);
        else localsize[j] = avgsize[j];
        localend[j] = localstart[j] + localsize[j];
    }

    ElementArray<Node> nodes(mesh);
//    for(int i = 0; i <= Nx; i++){
//        for(int j = 0; j <= Ny; j++){
//            for(int k = 0; k <= Nz; k++){
    for(int i = localstart[0]; i <= localend[0]; i++){
        for(int j = localstart[1]; j <= localend[1]; j++){
            for(int k = localstart[2]; k <= localend[2]; k++){
                double coords[3] = {i*dx, j*dy, k*dz};
                nodes.push_back(mesh->CreateNode(coords));
            }
        }
    }

//    for(int i = 1; i <= Nx; i++){
//        for(int j = 1; j <= Ny; j++){
//            for(int k = 1; k <= Nz; k++){
    for(int i = localstart[0]+1; i <= localend[0]; i++){
        for(int j = localstart[1]+1; j <= localend[1]; j++){
            for(int k = localstart[2]+1; k <= localend[2]; k++){
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

                double avX[3] = {0.,0.,0.};
                for(auto iv = verts.begin(); iv != verts.end(); iv++){
                    double x[3];
                    iv->Barycenter(x);
                    for(int i = 0; i < 3; i++)
                        avX[i] += x[i]/verts.size();
                }
                double r = 0.0;
                for(int i = 0; i < 2; i++)
                    r += (avX[i]-Rsample)*(avX[i]-Rsample);
                //printf("r = %lf\n", sqrt(r));
                if(sqrt(r) < Rsample)
                    mesh->CreateCell(verts,face_nodes,num_nodes,6); // Create the cubic cell in the mesh
            }
        }
    }
    //mesh->ResolveShared();

    //std::cout << "Created mesh with " << mesh->NumberOfCells() << " cells" << std::endl;
    times[T_MESHGEN] += Timer() - t;
}

void TwoPhaseFlow::initTags()
{
    double t = Timer();
	Sl        = mesh->CreateTag("Liquid_Saturation",     DATA_REAL,    CELL, false, 1);
	Pl        = mesh->CreateTag("Liquid_Pressure",       DATA_REAL,    CELL, false, 1);
	Pc        = mesh->CreateTag("Capillary_Pressure",    DATA_REAL,    CELL, false, 1);
	Pg        = mesh->CreateTag("Gas_Pressure",          DATA_REAL,    CELL, false, 1);
	Pf        = mesh->CreateTag("Fluid_Pressure",        DATA_REAL,    CELL, false, 1);
	Pf_old    = mesh->CreateTag("Fluid_Pressure_Old",    DATA_REAL,    CELL, false, 1);
	Pl_old    = mesh->CreateTag("Liquid_Pressure_Old",   DATA_REAL,    CELL, false, 1);
	Pg_old    = mesh->CreateTag("Gas_Pressure_Old",      DATA_REAL,    CELL, false, 1);
	X         = mesh->CreateTag("Primary_Variable",      DATA_REAL,    CELL, false, 1);
	Perm      = mesh->CreateTag("Permeability",          DATA_REAL,    CELL, false, 1);
	Phi       = mesh->CreateTag("Porosity",              DATA_REAL,    CELL, false, 1);
	Phi_old   = mesh->CreateTag("Porosity_Old",          DATA_REAL,    CELL, false, 1);
	PV        = mesh->CreateTag("Primary_Variable_Type", DATA_INTEGER, CELL, false, 1);
	Sl_old    = mesh->CreateTag("Liquid_Saturation_Old", DATA_REAL,    CELL, false, 1);
	TCoeff    = mesh->CreateTag("TPFA_Coefficient",      DATA_REAL,    FACE, false, 1);
	Grav      = mesh->CreateTag("Gravity",               DATA_REAL,    FACE, false, 1);
	fluxFaceL = mesh->CreateTag("fluxFaceL",             DATA_REAL,    FACE, false, 1);
	fluxFaceG = mesh->CreateTag("fluxFaceG",             DATA_REAL,    FACE, false, 1);
	Sltmp     = mesh->CreateTag("Sl_tmp",                DATA_REAL,    CELL, false, 1);
	Xtmp      = mesh->CreateTag("X_tmp",                 DATA_REAL,    CELL, false, 1);
	Pgtmp     = mesh->CreateTag("Pg_tmp",                DATA_REAL,    CELL, false, 1);
	Phitmp    = mesh->CreateTag("Phi_tmp",               DATA_REAL,    CELL, false, 1);
	BCtype    = mesh->CreateTag("BCtype",                DATA_INTEGER, FACE, FACE,  3);
	BCval     = mesh->CreateTag("BCval",                 DATA_REAL,    FACE, FACE,  3);

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
    return pow(1.0 + pow(vg_a/rhol0/g*Pcc, vg_n), -vg_m);
}

variable TwoPhaseFlow::get_Poro(variable PfP, double PfPn, double PhiPn)
{
    //return PhiPn + c_phi * (PfP - PfPn);
    return PhiPn / (1. - c_phi /* /PfP */ * (PfP - PfPn));
    //return phi0 * pow(PfP/0.1e6, -c_phi);
}

double TwoPhaseFlow::get_Poro(double PfP, double PfPn, double PhiPn)
{
    return get_Poro(variable(PfP), PfPn, PhiPn).GetValue();
}

variable TwoPhaseFlow::get_Ke(variable PfP)
{
    return exp(-gamma*(Pt - PfP - 0.1e6));
    //return pow(PfP/0.1e6, -c_phi);
}

variable TwoPhaseFlow::get_Pc(variable S)
{
    if(S.GetValue() < 0.0 || S.GetValue() > 1.0){
        mesh->Save("err" + outpExt);
        std::cout << "Bad saturation " << S.GetValue() << std::endl;
        exit(1);
    }
    return rhol0*9.81/vg_a * pow(pow(S,-1./vg_m) - 1., 1./vg_n);
}

void TwoPhaseFlow::get_Kr(variable S, variable &Krl, variable &Krg)
{
    double s = S.GetValue();

    if(s > 1.0 || s < 0.0 || std::isinf(s) || std::isnan(s)){
        std::cout << "Bad s = " << s << std::endl;
        exit(1);
    }

    Krl = pow(S, 2.0);
    double Krmin = 1e-10;
    if(Krl.GetValue() < Krmin)
        Krl = Krmin;
    //Krl = sqrt(S) * pow(1.-pow(1.-pow(S,1./vg_m),vg_m),2.0);

    Krg = pow(1.-S, 2.0);
    //Krg = pow(1.-S,2.0)*(1.-S*S);
}

variable TwoPhaseFlow::get_rhol(variable pl)
{
    double Pl0 = Pg0 - get_Pc(Sl0).GetValue();
    return rhol0 + c_w * (pl - Pl0);
}

void TwoPhaseFlow::assembleResidual()
{
    double t = Timer();
    for(auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
        if(icell->GetStatus() == Element::Ghost) continue;
        Cell cellP = icell->getAsCell();

        // Values needed for accumulation part
        variable SP, PfP, PlP, massL, massG, PhiP, rholP;
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

        rholP = get_rhol(PlP);

        PfP           = SP * PlP + (1.-SP) * varPg(cellP);
        PhiPn         = cellP.Real(Phi_old);
        PhiP          = get_Poro(PfP, cellP.Real(Pf_old), cellP.Real(Phi_old));

        massL         = rholP *     SP   * PhiP;//varPhi(cellP);
        massG         = rhog * (1.-SP)  * PhiP;//varPhi(cellP);
        massL_old     = rholP.GetValue() * SPn  * cellP.Real(Phi_old);
        massG_old     = rhog * (1.-SPn) * cellP.Real(Phi_old);

        R[varX.Index(cellP)]   = (massL - massL_old)/dt;// * V;
        R[varPg.Index(cellP)]  = (massG - massG_old)/dt;// * V;
        //R[varPg.Index(cellP)]  = varPg(cellP) - Pg0;

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
					face.Real(fluxFaceL) = ql.GetValue();
                }
                else if(faceBCtypeL == BC_DIR){
                    std::cout << "Face with Dirichlet BC for liquid" << std::endl;
                    double PlBC = face.RealArray(BCval)[BCAT_L];
                    variable Krl, Krg;
                    get_Kr(SP, Krl, Krg);
                    //variable Ke = exp(-gamma*(Pt - PfP - 0.1e6));
                    variable Ke = get_Ke(PfP);

                    double coef = face.Real(TCoeff);

                    ql = -get_rhol(PlBC)*Krl*KP*Ke/mul * coef * (PlP - PlBC);
					face.Real(fluxFaceL) = ql.GetValue();
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
                    //SP = 1.0;//
                    PgBC = PBC + SP*get_Pc(SP);
                    PlBC = (PBC - (1.-SP)*PgBC)/SP;
                    variable Krl, Krg;
                    get_Kr(SP, Krl, Krg);
                    variable Ke = exp(-gamma*(Pt - PfP - 0.1e6));

                    double coef = face.Real(TCoeff);

                    ql = -get_rhol(PlBC)*Krl*KP*Ke/mul * coef * (PlP - PlBC);
                    qg = -rhog*Krg*KP*Ke/mug * coef * (varPg(cellP) - PgBC);

					face.Real(fluxFaceL) = ql.GetValue();
					face.Real(fluxFaceG) = qg.GetValue();
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

//                if(Krg.GetValue() < 1e-9)
//                    Krg = 1e-9;

                double K, KN = cellN.Real(Perm);

                K = 2.0 * KP*KN / (KP + KN);

                variable rhol =
                        (PlP-PlN).GetValue() > 0 ?
                            rholP
                          : get_rhol(PlN);

                ql = -rhol*Krl*K*Ke/mul * coef * (PlP - PlN);
                qg = -rhog*Krg*K*Ke/mug * coef * (varPg(cellP) - varPg(cellN));

				face.Real(fluxFaceL) = ql.GetValue();
				face.Real(fluxFaceG) = qg.GetValue();

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

    Tag Ref = mesh->CreateTag("REF", DATA_REAL, CELL, NONE, 1);

    for(auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
        if(icell->GetStatus() == Element::Ghost) continue;


        double x[3];
        icell->Barycenter(x);
        //icell->Real(Ref) = mul/rhol/K0*inflowFluxL/inflowArea*(x[2]-Lsample) + outflowPresF;

        //if(!loadMesh)
            icell->Real(Perm) = K0;// * (0.5 + rand()/(double)RAND_MAX);

        icell->Real(Pg) = Pg0;
        icell->Real(Sl) = Sl0;
        if(problem_name == "2phase_center"){
            double x[3];
            icell->Barycenter(x);
            double r = (x[0]-Rsample)*(x[0]-Rsample) + (x[1]-Rsample)*(x[1]-Rsample);
            r = sqrt(r);
            if(r < Rsample / 6.)
                icell->Real(Sl) = Sl0_c;
        }
        double S = icell->Real(Sl);
        icell->Real(Pc) = (get_Pc(S)).GetValue();
        icell->Real(Pl) = 8e6;//icell->Real(Pg) - icell->Real(Pc);
        icell->Real(Pg) = icell->Real(Pl) + (get_Pc(S)).GetValue();
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
//            KK = mesh->GetTag("Perm");
//        if(mesh->HaveTag("Permeability"))
//            KK = mesh->GetTag("Permeability");
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
            if(r < 200.)
                icell->Real(Sl) = Sl0_c;
            else {
                icell->Real(Sl) = Sl0;
            }

            icell->Real(Perm) = 1e-10*icell->RealArray(KK)[0];
            icell->Real(Phi) = icell->Real(PORO);
            double S = icell->Real(Sl);
            icell->Real(Pl) = icell->Real(Pg) - (get_Pc(S)).GetValue();
            icell->Real(Pc) = icell->Real(Pg) - icell->Real(Pl);
            mass += get_rhol(icell->Real(Pl)).GetValue() * S * icell->Real(Phi) * icell->Volume();
        }
    }

    mass = mesh->Integrate(mass);
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
    inflowFaces.SetMeshLink(mesh);
    inflowCells.SetMeshLink(mesh);
    outflowFaces.SetMeshLink(mesh);
    inflowArea = 0.0;
    iclsize = 0.0;
    Tag ref = mesh->GetTag("REF");
    for(auto iface = mesh->BeginFace(); iface != mesh->EndFace(); iface++){
        Face face = iface->getAsFace();
        if(!face.Boundary())
            continue;
        double x[3];
        face.Barycenter(x);

        // Upper boundary - hardcoded
        if(fabs(x[2]-Lsample) < 1e-7){
            double r = (x[0]-Rsample)*(x[0]-Rsample) + (x[1]-Rsample)*(x[1]-Rsample);
            r = sqrt(r);
            //printf("r = %lf\n", r);
            if(r <= injectRadius){
                inflowFaces.push_back(face);
                inflowCells.push_back(face.BackCell());
                if(face.BackCell().GetStatus() != Element::Ghost){
                    iclsize += 1.0;
                    inflowArea += face.Area();
                }
                face.BackCell().Real(ref) = -1e20;
            }
        }

        // Lower boundary - hardcoded
        if(fabs(x[2]-0.0) < 1e-7){
            double r = (x[0]-Rsample)*(x[0]-Rsample) + (x[1]-Rsample)*(x[1]-Rsample);
            r = sqrt(r);
            if(r <= injectRadius){
                //std::cout << "Bottom boundary face " << face.GlobalID() << ", z = " << x[2] << std::endl;
                face.IntegerArray(BCtype)[BCAT_F] = BC_DIR;
                face.RealArray(BCval)[BCAT_F] = outflowPresF;
                outflowFaces.push_back(face);
            }
        }
    }

    iclsize = mesh->Integrate(iclsize);//static_cast<double>(inflowCells.size()));
    inflowArea = mesh->Integrate(inflowArea);
    if(rank == 0){
        std::cout << "Number of inflow cells " << iclsize << std::endl;
        std::cout << "Inflow area: " << inflowArea << std::endl;
    }


    for(auto iface = inflowFaces.begin(); iface != inflowFaces.end(); iface++){
        if(iface->GetStatus() == Element::Ghost)
            continue;

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
        double rhol = (get_rhol(icell->Real(Pl))).GetValue();
        mass_new += rhol * icell->Real(Sl) * icell->Real(Phi) * icell->Volume();
    }
    mass_new = mesh->Integrate(mass_new);
    if(rank == 0){
        std::cout << "Mass change is " << mass_new-mass;
        std::cout << " (" << (mass_new-mass)/mass*1e2 << "%)" << std::endl;
    }
    mass = mass_new;
}

void TwoPhaseFlow::computeFluxes()
{
    Tag Lflux, Gflux;
    if(mesh->HaveTag("Liquid_Flux"))
        Lflux = mesh->GetTag("Liquid_Flux");
    else
        Lflux = mesh->CreateTag("Liquid_Flux", DATA_REAL, CELL, NONE, 3);
    if(mesh->HaveTag("Gas_Flux"))
        Gflux = mesh->GetTag("Gas_Flux");
    else
        Gflux = mesh->CreateTag("Gas_Flux", DATA_REAL, CELL, NONE, 3);

    for(auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
        if(icell->GetStatus() == Element::Ghost)
            continue;
        Cell c = icell->getAsCell();
        auto faces = icell->getFaces();
        unsigned m = static_cast<unsigned>(faces.size());
        double cfluxL[3] = {0., 0., 0.};
        double cfluxG[3] = {0., 0., 0.};
        for(unsigned k = 0; k < m; k++){
            double nor[3];
            faces[k].UnitNormal(nor);
            for(int i = 0; i < 3; i++){
				cfluxL[i] += nor[i] * faces[k].Real(fluxFaceL);
				cfluxG[i] += nor[i] * faces[k].Real(fluxFaceG);
            }
        }
        for(unsigned i = 0; i < 3; i++){
            icell->RealArray(Lflux)[i] = cfluxL[i]/m;
            icell->RealArray(Gflux)[i] = cfluxG[i]/m;
        }

    }
    mesh->ExchangeData(Lflux, CELL);
    mesh->ExchangeData(Gflux, CELL);
}

bool TwoPhaseFlow::makeTimeStep()
{

    double r2, r2_0;
    double t;

    Sparse::Vector sol("Newton_sol", aut->GetFirstIndex(), aut->GetLastIndex());

    if(rank == 0){
        std::cout << "Newton: maxit = " << maxit;
        std::cout << ", rtol = " << rtol << ", atol = " << atol << std::endl;
        std::cout << "Solver: " << solver_type << std::endl;
    }
    bool converged = false;
    int iter;
    for(iter = 0; iter < maxit; iter++){
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
            if(rank == 0){
                std::cout << "Linear solver failed: " << S->GetReason() << std::endl;
                std::cout << "Residual: " << S->Residual() << std::endl;
            }
            //exit(1);
            return false;
        }
        if(rank == 0)
            std::cout << "  Lin.it: " << S->Iterations() << std::endl;
        iterLinear += S->Iterations();
        times[T_SOLVE] += Timer() - t;

        t = Timer();
        // Line search loop
        copyTagReal(Sltmp, Sl, CELL);
        copyTagReal(Pgtmp, Pg, CELL);
        copyTagReal(Phitmp, Phi, CELL);
        bool lsSuccess = false;
        w = 1.0;
        for(int ils = 0; ils < 20; ils++){
            int gotBad = 0;
            for(auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
                if(icell->GetStatus() == Element::Ghost) continue;

                Cell cell = icell->getAsCell();
                cell.Real(Pg) -= w*sol[varPg.Index(cell)];
                if(cell.Integer(PV) == PV_SAT){
                    double Snew = cell.Real(Sl) - w*sol[varX.Index(cell)];
                    if(Snew > 1.0 || Snew < 0.0){
                        //std::cout << "    Bad Sl = " << Snew << " at cell " << cell.GlobalID() << std::endl;
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
                cell.Real(Phi) = get_Poro(cell.Real(Pf), cell.Real(Pf_old), cell.Real(Phi_old));//cell.Real(Phi_old) + c_phi/cell.Real(Pf)*(cell.Real(Pf)-cell.Real(Pf_old));
                if(cell.Real(Phi) > 1e5){
                    //printf("Bad porosity %e at cell %d\n", cell.Real(Phi), cell.LocalID());
                    return false;
                }
            }

            //std::cout << "Here" << std::endl;
            gotBad = mesh->Integrate(gotBad);
            if(gotBad){
                w *= 0.05;

                if(w < 1e-8){
                    lsSuccess = false;
                    break;
                }

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
            if(rank == 0)
                std::cout << "Line search failed" << std::endl;
            //mesh->Save("err" + outpExt);
            //exit(1);
            return false;
        }
        times[T_UPDATE] += Timer() - t;
        mesh->ExchangeData(Pl, CELL);
        mesh->ExchangeData(Pg, CELL);
        mesh->ExchangeData(Pc, CELL);
        mesh->ExchangeData(Pf, CELL);
        mesh->ExchangeData(Sl, CELL);
        mesh->ExchangeData(Phi, CELL);
        iterNewton++;
    }
    if(!converged){
        if(rank == 0)
            std::cout << "Newton failed" << std::endl;
//        mesh->Save("err" + outpExt);
//        exit(1);
        return false;
    }
    if(iter < 7)
        dtnew = dt * 1.5;
    else
        dtnew = dt;
    return true;
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
    std::ofstream out1("t.txt");
    std::ofstream out2("q.txt");

    initAutodiff();

    S = new Solver(solver_type);
    S->SetParameter("absolute_tolerance","1e-10");
    S->SetParameter("relative_tolerance","1e-6");
    S->SetParameter("maximum_iterations","100");
    S->SetParameter("gmres_substeps","0");
    S->SetParameter("condition_estimation","1");
    S->SetParameter("schwartz_overlap","3");
    S->SetParameter("drop_tolerance","1e-1");

    int nt = static_cast<int>(T/dt);
    dt = 1e-1;
    double tbeg = 0.0;
    for(int it = 1; ; it++){
        if(rank == 0){
            std::cout << std::endl;
            std::cout << "===== TIME STEP " << it << ", T_beg = " << tbeg << " =====" << std::endl;
        }
        // Save old values
        copyTagReal(Sl_old, Sl, CELL);
        copyTagReal(Pf_old, Pf, CELL);
        copyTagReal(Pl_old, Pl, CELL);
        copyTagReal(Pg_old, Pg, CELL);
        copyTagReal(Phi_old, Phi, CELL);
        bool found = false;
        for(dt = std::min(dtDesir, std::min(T-tbeg, dt)); dt >= 1e-10; dt *= 0.5){
            if(rank == 0)
                std::cout << "Trying dt = " << dt << std::endl;
            bool success = makeTimeStep();
            if(success){
                found = true;
                break;
            }
            else {
                copyTagReal(Sl, Sl_old, CELL);
                copyTagReal(Pf, Pf_old, CELL);
                copyTagReal(Pl, Pl_old, CELL);
                copyTagReal(Pg, Pg_old, CELL);
                copyTagReal(Phi, Phi_old, CELL);
            }
        }
        if(!found){
            std::cout << "Failed to find dt" << std::endl;
            exit(2);
        }
        countMass();
        computeFluxes();


        t = Timer();
        if(it%saveIntensity == 0 && saveSol){
            mesh->Save(save_dir + "/sol" + std::to_string(it/saveIntensity) + outpExt);
        }
        //printf("proc %d of %d, size = %lf\n", rank, mesh->GetProcessorsNumber(), iclsize);
        //MPI_Barrier(INMOST_MPI_COMM_WORLD);
        if(iclsize > 0){
            double avPin = 0.0;
            for(auto icell = inflowCells.begin(); icell != inflowCells.end(); icell++){
                if(icell->GetStatus() != Element::Ghost)
                avPin += icell->Real(Pf);
            }
            avPin = mesh->Integrate(avPin);
            double Q = 0.0; // outflow flux
            for(auto iface = outflowFaces.begin(); iface != outflowFaces.end(); iface++){
                if(iface->GetStatus() == Element::Ghost)
                    continue;
                Face f = iface->getAsFace();
                Cell c = f.BackCell();
                variable Krl, Krg;
                get_Kr(c.Real(Sl), Krl,Krg);
                variable Ke = get_Ke(c.Real(Pf));
                double k = c.Real(Perm);
                double t = f.Real(TCoeff);
                double bcL = f.RealArray(BCval)[BCAT_L];
                double bcG = f.RealArray(BCval)[BCAT_G];
                variable rhol = get_rhol(c.Real(Pl));
                double fluxL = -(k * rhol * Krl / mul * t * (bcL - c.Real(Pl))).GetValue();
                double fluxG = -(k * rhol * Krl / mul * t * (bcG - c.Real(Pg))).GetValue();
                Q += fluxL + fluxG;
            }
            Q = mesh->Integrate(Q);
            if(rank == 0){
                printf("avPin = %e\n", avPin/iclsize);
                out << 1e-6*avPin/iclsize << std::endl;
                out1 << tbeg+dt << std::endl;// << "\t" << 1e-6*avPin/iclsize << std::endl;
                out2 << Q << std::endl;
            }
        }
        times[T_IO] += Timer() - t;

        tbeg += dt;
        if(tbeg > T - 1e-5)
            break;

        dt = dtnew;
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
