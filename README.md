# TwoPhaseFlow
This project contains code for finite volume modeling of two-phase fluid flow in deformable porous medium.
Currently many parameters are hardcoded since the code is mainly used for one problem: flow in a cylinder. However, the code
supports arbitrary domains and polyhedral meshes.

INMOST (https://github.com/INMOST-DEV/INMOST) is needed for this project.

The following executables are produced:

- twophase: main executable that performs simulation
- gridgen: executable that creates grid (currently: meshes with cubic cells that approximate cylinder domain)
- gridrefine: executable that takes a mesh and refines it, generally producing mesh with 8x cells
- gmshgrid: executable that takes externally produced mesh (from generators like Gmsh) and prepares it