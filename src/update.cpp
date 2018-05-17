/*
 update.cpp

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

#include "update.h"
#include "potential.h"
#include "neighbor.h"
#include "memory.h"
#include "input.h"
#include "compute.h"

using namespace std;

using namespace EM3_NS;

Update::Update(EM3 *em3) : Pointers(em3) {

    dt = input->dt;
    natoms = input->natoms;

    memory->allocate(v, natoms, 3);

    for (int i=0; i<natoms; i++){
        for (int j=0; j<3; j++){
            v[i][j] = input->v[i][j];
        }
    }

    //fh_debug = fopen("DEBUG_update", "w");
}

Update::~Update() 
{

    memory->deallocate(v);

    //fclose(fh_debug);

};

void Update::integrate()
{

    // In LJ units, the accelerations are simply the forces

    double **a = potential->f;
    double **x = neighbor->x;

    for (int i=0; i<natoms; i++){

        v[i][0] += 0.5*dt*a[i][0];
        v[i][1] += 0.5*dt*a[i][1];
        v[i][2] += 0.5*dt*a[i][2];

        x[i][0] += dt*v[i][0];
        x[i][1] += dt*v[i][1];
        x[i][2] += dt*v[i][2];

    }

    // Update neighborlist

    em3->neighbor->generate();

    // Calculate forces with these new positions

    em3->potential->calculate();

    a = potential->f;

    // Perform final step of Verlet algorithm

    for (int i=0; i<natoms; i++){
        v[i][0] += + 0.5*dt*a[i][0];
        v[i][1] += + 0.5*dt*a[i][1];
        v[i][2] += + 0.5*dt*a[i][2];
    }

    /*
    Rescale the velocities if the timestep is less than the requested equilibration time
    */

    if (em3->t < input->nequil){
        compute->compute_ke();
        double temp_init = input->temp;
        double temp_current = compute->temp;
        for (int i=0; i<natoms; i++){
            v[i][0] = v[i][0]*sqrt(temp_init/temp_current);
            v[i][1] = v[i][1]*sqrt(temp_init/temp_current);
            v[i][2] = v[i][2]*sqrt(temp_init/temp_current);
        }
    }

}
