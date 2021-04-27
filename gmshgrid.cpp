#include "header.h"


void interpPerm(Mesh *m){
    Mesh *mp = new Mesh("");
    mp->Load("refperm.vtk");

    Tag K = m->CreateTag("Permeability", DATA_REAL, CELL, NONE, 1);
    Tag PermRef = mp->GetTag("Permeability");

    for(auto icell = m->BeginCell(); icell != m->EndCell(); icell++){
        double xi[3];
        printf("Processing cell %d of %d\n",icell->DataLocalID(), m->NumberOfCells());
        icell->Barycenter(xi);
        double rmin = 1e20;
        double prm;
        for(auto jcell = mp->BeginCell(); jcell != mp->EndCell(); jcell++){
            double xj[3];
            jcell->Barycenter(xj);
            double r = (xi[0]-xj[0])*(xi[0]-xj[0]) + (xi[1]-xj[1])*(xi[1]-xj[1]) + (xi[2]-xj[2])*(xi[2]-xj[2]);
            r = sqrt(r);
            if(r < rmin){
                rmin = r;
                prm = jcell->Real(PermRef);
            }
        }
        icell->Real(K) = prm;
    }
}

int main(int argc, char *argv[])
{
    if(argc != 3){
        std::cout << "Usage: gmshgrid <gmsh_mesh_file> <out_mesh_file>" << std::endl;
        exit(1);
    }

    Mesh::Initialize(&argc, &argv);
    Partitioner::Initialize(&argc, &argv);

    //Mesh *mesh = new Mesh("mesh");
    Mesh *m    = new Mesh("gmsh_processed");

    m->Load(argv[1]);

    std::vector<std::string> tagNames;
    tagNames.push_back("K");
    tagNames.push_back("Permeability_scalar");
    tagNames.push_back("Permeability");
    tagNames.push_back("PORO");

    m->AssignGlobalID(CELL|FACE);
    //interpPerm(m);
//    Tag K = m->CreateTag("Permeability", DATA_REAL, CELL, NONE, 1);
//    srand(time(nullptr));
//    for(auto icell = m->BeginCell(); icell != m->EndCell(); icell++){
//        icell->Real(K) = 1e-18 * (0.5 + rand()/(double)RAND_MAX);
//    }


    std::cout << "Generated mesh with " << m->NumberOfCells() << " cells" << std::endl;

    m->Save(argv[2]);

    delete m;

    Mesh::Finalize();
    Partitioner::Finalize();
}
