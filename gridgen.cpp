#include "header.h"

#define V_ID(x, y, z) static_cast<unsigned long long>(((x-0)*(Ny+1)*(Nz+1) + (y-0)*(Nz+1) + (z-0)))

int main(int argc, char *argv[])
{
    Mesh::Initialize(&argc, &argv);

    if(argc != 5){
        std::cout << "Usage: gridgen <mesh_name> <Nx> <Ny> <Nz>" << std::endl;
        exit(1);
    }

    int Nx = atoi(argv[2]);
    int Ny = atoi(argv[3]);
    int Nz = atoi(argv[4]);

    double L = 0.012;
    double dx = L/Nx, dy = L/Ny, dz = L/Nz;

    Mesh *mesh = new Mesh;
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

                double avX[3] = {0.,0.,0.};
                for(auto iv = verts.begin(); iv != verts.end(); iv++){
                    double x[3];
                    iv->Barycenter(x);
                    for(int i = 0; i < 3; i++)
                        avX[i] += x[i]/verts.size();
                }
                double r = 0.0;
                for(int i = 0; i < 2; i++)
                    r += (avX[i]-0.006)*(avX[i]-0.006);
                //printf("r = %lf\n", sqrt(r));
                if(sqrt(r) < 0.006)
                    mesh->CreateCell(verts,face_nodes,num_nodes,6); // Create the cubic cell in the mesh
            }
        }
    }
    std::cout << "Created mesh with " << mesh->NumberOfCells() << " cells" << std::endl;
    mesh->Save(argv[1]);

    Mesh::Finalize();
}
