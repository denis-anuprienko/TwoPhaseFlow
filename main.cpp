#include <cstdio>
#include "header.h"

TwoPhaseFlow::TwoPhaseFlow()
{
    setDefaultParams();
}

TwoPhaseFlow::~TwoPhaseFlow()
{

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
}

void TwoPhaseFlow::readMesh(std::string path)
{
    mesh.Load(path);
}

void TwoPhaseFlow::cleanMesh()
{
    // Tags that are likely to be on mesh
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

    //mesh.Save("out.vtk");
}

int main(int argc, char *argv[])
{
    if(argc != 3){
        std::cout << "Usage: twophase <param_file_path> <mesh_file_path>" << std::endl;
    }
    TwoPhaseFlow Problem;
    Problem.readMesh(argv[2]);
    Problem.cleanMesh();
    return 0;
}
