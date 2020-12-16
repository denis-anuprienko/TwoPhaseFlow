#include "header.h"

int main(int argc, char *argv[])
{
    if(argc != 3){
        std::cout << "Usage: gridrefine <mesh_file> <refinement_level>" << std::endl;
        exit(1);
    }
    std::cout << "Only cubic-like grids are supported" << std::endl;

    Mesh::Initialize(&argc, &argv);
    Partitioner::Initialize(&argc, &argv);

    Mesh *mesh = new Mesh("mesh");

    int rank = mesh->GetProcessorRank();
    std::cout << "rank = " << rank << std::endl;
    MPI_Barrier(INMOST_MPI_COMM_WORLD);

    if(rank == 0)
        mesh->Load(argv[1]);

    Partitioner p(mesh);
    p.SetMethod(Partitioner::INNER_KMEANS, Partitioner::Partition);
    p.Evaluate();

    mesh->Redistribute();
    mesh->ExchangeGhost(1, NODE);

    Tag T = mesh->CreateTag("My", DATA_REAL, CELL, false, 1);
    T.SetPrint(true);
    for(auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
        if(icell->GetStatus() == Element::Ghost)
            continue;
        double x[3];
        icell->Barycenter(x);
        icell->Real(T) = x[0]+x[1]+x[2];
    }
    mesh->Save("parallel.pvtk");

    delete mesh;

    Mesh::Finalize();
    Partitioner::Finalize();
}
