#include "header.h"

int main(int argc, char *argv[])
{
    if(argc != 3){
        std::cout << "Usage: gridrefine <mesh_file> <refinement_level>" << std::endl;
        exit(1);
    }
    std::cout << "Only cubic-like grids are supported" << std::endl;

    Mesh m, mesh;
    m.Load(argv[1]);

    int refLev = atoi(argv[2]);

    for(int iref = 0; iref < refLev; iref++){
        std::cout << std::endl << "Refinement level " << iref << std::endl;

        for(auto inode = m.BeginNode(); inode != m.EndNode(); inode++){
            Storage::real_array crds = inode->Coords();
            double coords[3] = {crds[0], crds[1], crds[2]};
            mesh.CreateNode(coords);
        }

        for(auto icell = m.BeginCell(); icell != m.EndCell(); icell++){
            //mesh.CreateCell(icell->getFaces(), icell->getNodes());
            std::cout << "Added cell " << icell->DataLocalID() << std::endl;
        }
    }
    mesh.Save("new.vtk");
}
