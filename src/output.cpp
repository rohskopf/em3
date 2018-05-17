/*
 output.cpp

 Copyright (c) 2018 Andrew Rohskopf

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include "mpi.h"
#include <math.h>       /* sqrt */
#include <random>

#include "output.h"
#include "potential.h"
#include "input.h"
#include "neighbor.h"
#include "memory.h"

using namespace std;

using namespace EM3_NS;

Output::Output(EM3 *em3) : Pointers(em3) {

    fh_xyz = fopen("dump.xyz", "w");
    fh_forces = fopen("forces.output", "w");


}

Output::~Output() 
{
    fclose(fh_xyz);
    fclose(fh_forces);

};

void Output::write_xyz()
{


    fprintf(fh_xyz, "%d\n", input->natoms);
    fprintf(fh_xyz, "Atoms. Timestep: %d\n", em3->t);
    for (int i=0; i<input->natoms; i++){
        fprintf(fh_xyz, "1 %f %f %f\n", neighbor->x[i][0], neighbor->x[i][1], neighbor->x[i][2]);
    }


}

void Output::write_forces()
{

    for (int i=0; i<input->natoms; i++){
        fprintf(fh_forces, "%f %f %f\n", potential->f[i][0], potential->f[i][1], potential->f[i][2]);
    }


}
