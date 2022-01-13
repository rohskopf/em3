# Easy Modular Molecular Mechanics (EM3) Program

A plethora of powerful open-source atomistic simulation programs are available to the public, including LAMMPS, GROMACS, and many others. With the current availability of powerful atomistic simulation programs, the development of another might seem futile. EM3, however, possesses a goal that separates it from the rest of the currently available programs. Unlike existing programs possessing goals centered around broad performance across all materials/molecules/systems, the goal of EM3 is more educational and serves as a simple platform for studyign MD simulation *performance*. This educational goal centers around the simple (easy) objective-oriented (modular) C++ source code of EM3, which aside from being powerful and efficient, seeks to bridge the gap between less powerful but easily understood MATLAB and Python codes scattered about the internet and more powerful but complicated larger programs like LAMMPS and GROMACS (also written in C++ and C, respectively). 

### Installation

In a Linux environment, build with:

    cd src
    make clean
    make

This creates an executable called `em3`.

### Using EM3

The manual is located here: [https://github.com/rohskopf/em3/blob/master/Manual.pdf](https://github.com/rohskopf/em3/blob/master/Manual.pdf)

There are two main input files: (1) `INPUT` and (2) `CONFIGS`. These files must be placed in the directory that the simulation occurs, and then execute the `em3` exectuable in this directory.

Example `INPUT` and `CONFIG` files for each of the cases discussed in the manual given in the 
`/examples` directory. 
