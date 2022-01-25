# Design-of-Metamaterials-for-Optics

This repository contains all computational resources that were used for computations in the aforementioned publication. All computations are done in MATLAB using FELICITY (Finite ELement Implementation and Computational Interface Tool for You) which is freely available at https://www.mathworks.com/matlabcentral/fileexchange/31141-felicity/.


## Dependencies

1. The [FELICITY](https://www.mathworks.com/matlabcentral/fileexchange/31141-felicity/) package needs to be installed in the current MATLAB directory.
2. The [AGMG](http://agmg.eu) needs to be downloaded and placed in the current MATLAB directory.
3. The  Matlab function `licols.m` which is freely available [HERE](https://www.mathworks.com/matlabcentral/fileexchange/77437-extract-linearly-independent-subset-of-matrix-columns) is used in the Permeability_Function.m.
4. The Matlab function `load_gmsh2.m` which is freely available [HERE](https://github.com/cycheung/gmsh/blob/master/utils/converters/matlab/load_gmsh2.m) is used in both  Permeability_Function.m and Permittivity_Function.m. 

## Folder Information

- **Permeability_Function**: This folder contains  the Matlab function `Permeability_Function.m` that computes the effective magnetic permeability tensor by solving the relevant eigenvalue problem. It also contains all the linear and bilinear forms written in FELICITY to assemble matrices needed in the process of executing `Permeability_Function.m`. 

- **Permittivity_Function**: This folder contains  the Matlab function `Permittivity_Function.m` that computes the effective dielectric permittivity tensor by solving the relevant eigenvalue problem. It also contains all the linear and bilinear forms written in FELICITY to assemble matrices needed in the process of executing `Permittivity_Function.m`. 

- **Mesh Files for Different Geometries**: This folder contains the Mesh files generated using GMSH. Currently, both `Permeability_Function.m` and `Permittivity_Function.m` are set up to use the mesh file `con0.3(0.12)`. If  another geometry is chosen to  be used  from this folder, the corresponding parts of `Permeability_Function.m` and `Permittivity_Function.m` including `All_Mesh4.m` need to be manually adjusted. 
