/*
 input.h

 Copyright (c) 2018 Andrew Rohskopf

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <vector>
#include <string>
#include "mpi.h"

#include <iostream>
#include <new>
#include <cstdlib>
#include "pointers.h"

using namespace std;

namespace EM3_NS
{
  class Input: protected Pointers
  {
  public:
    Input(class EM3 *);
    ~Input();

    FILE * fh_debug; // Debug file handle

    // Declare member functions

    void readinput(); // function to read INPUT file
    void readconfig(); // function to read CONFIG file
    void initialize(); // function to initialize velocities and convert inputs

    // readinput() variables

    int nsteps; // number of timesteps 
    int nequil; // number of NVT equilibration timesteps
    int nout; // output data this many timesteps
    int neighmax; // maximum number of neighbors per atom
    double rc; // cutoff
    int neighcount;
    int newton;
    int offset; // whether or not to calculate energy offset
    double m; // mass
    double temp;
    double dt; // timestep
    double epsilon;
    double sigma;
   

    // readconfig() variables

    int natoms;
    double box[3];
    double **x; // positions
    double **x0; // original positions to store

    // initilize() variables

    double tau; // LJ time unit
    double velocity; // LJ velocity unit
    double force; // LJ force unit
    double pressure; // LJ pressure unit
    double temperature; // LJ temperature unit
    double **v; // velocities
    double volume; // box volume
    

  };
}

