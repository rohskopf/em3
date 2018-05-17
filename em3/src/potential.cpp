/*
 potential.cpp

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

#include "potential.h"
#include "neighbor.h"
#include "memory.h"
#include "input.h"
#include "compute.h"

using namespace std;

using namespace EM3_NS;

Potential::Potential(EM3 *em3) : Pointers(em3) {

    natoms = input->natoms;

    memory->allocate(f, natoms, 3);

    for (int i=0; i<natoms; i++){
        for (int j=0; j<3; j++){
            f[i][j] = 0.0;
        }
    }

    //fh_debug = fopen("DEBUG_potential", "w");
}

Potential::~Potential() 
{

    memory->deallocate(f);

    //fclose(fh_debug);

};

void Potential::calculate()
{

    
    double **x = neighbor->x;
    int **neighlist = neighbor->neighlist;
    int **neightags = neighbor->neightags;
    int *numneigh = neighbor->numneigh;
    double rc = input->rc;
    double cutsq = rc*rc;

    /* Compute offset, if requested. */
    double offset;
    if (input->offset == 0){
        double rinv2_off = (1.0/rc)*(1.0/rc);
        double rinv6_off = rinv2_off*rinv2_off*rinv2_off;
        offset = 4.0*rinv6_off*(rinv6_off - 1.0);
    }
    else if (input->offset == 1){
        offset = 0.0;
    }
    
    int i,j, itag, jtag;
    double xij,yij,zij,rsq,rij;

    double xtmp, ytmp,ztmp;

    double rinv2, rinv6, uij, fpair;

    int nearintx, nearinty, nearintz;
    double xlength = input->box[0];
    double ylength = input->box[1];
    double zlength = input->box[2];
    double boxinvx = 1.0/xlength;
    double boxinvy = 1.0/ylength;
    double boxinvz = 1.0/zlength;

    pe = 0.0; // potential energy

    for (int i=0; i<natoms; i++){
        for (int j=0; j<3; j++){
            f[i][j] = 0.0;
        }
    }

    p_xx = 0.0; // total stress tensor component
    p_yy = 0.0; // total stress tensor component
    p_zz = 0.0; // total stress tensor component
    double p_xx_i, p_yy_i, p_zz_i; // stress tensor components for atom i

    for (int ii=0; ii<natoms; ii++){

        xtmp = x[ii][0];
        ytmp = x[ii][1];
        ztmp = x[ii][2];

        p_xx_i = 0.0;
        p_yy_i = 0.0;
        p_zz_i = 0.0;
   
        for (int jj=0; jj<numneigh[ii]; jj++){

            j = neighlist[ii][jj];
            jtag = neightags[ii][jj];

            xij=xtmp-x[j][0];
            yij=ytmp-x[j][1];
            zij=ztmp-x[j][2];

            rsq = xij*xij + yij*yij + zij*zij;

            if (rsq < cutsq){
            
                rij = sqrt(rsq);

                rinv2 = 1.0/rsq;
                rinv6 = rinv2*rinv2*rinv2;

                /*
                If we have a full neighbor list, the LJ functional form is half the original form.
                Apply this change to energy and forces, depending on the request neighbor list.
                */

                if (input->neighcount == 0){
                    uij = 4.0*rinv6*(rinv6 - 1.0) - offset;
    
                }
                else if (input->neighcount == 1){
                    uij = 2.0*rinv6*(rinv6 - 1.0) - offset;
                }
                pe += uij;

                if (input->neighcount == 0){
                    fpair = 24.0*rinv6*rinv2*(2.0*rinv6 - 1.0);
                }
                else if (input->neighcount == 1){
                    fpair = 12.0*rinv6*rinv2*(2.0*rinv6 - 1.0);
                }

                f[ii][0] += xij*fpair;
                f[ii][1] += yij*fpair;
                f[ii][2] += zij*fpair;

                // Apply Newton's 2nd law

                if (input->newton == 0){
                    f[jtag][0] -= xij*fpair;
                    f[jtag][1] -= yij*fpair;
                    f[jtag][2] -= zij*fpair;
                }

                // Calculate pressure tensor 

                p_xx_i += xij*xij*fpair;
                p_yy_i += yij*yij*fpair;
                p_zz_i += zij*zij*fpair;

            } // if (rsq < cutsq){

        } // for (int jj=0; jj<numneigh[ii]; jj++){

        p_xx += p_xx_i;
        p_yy += p_yy_i;
        p_zz += p_zz_i;

    } // for (int ii=0; ii<natoms; ii++){

}
