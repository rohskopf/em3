/*
 neighbor.cpp

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

#include "neighbor.h"
#include "memory.h"
#include "input.h"

using namespace std;

using namespace EM3_NS;

Neighbor::Neighbor(EM3 *em3) : Pointers(em3) {


    natoms = em3->input->natoms;
    neighmax = em3->input->neighmax;

    memory->allocate(numneigh, natoms);
    memory->allocate(neighlist, natoms, neighmax);
    memory->allocate(neightags, natoms, neighmax);
    memory->allocate(x, natoms*neighmax, 3);
    
    for (int i=0; i<natoms; i++){
        numneigh[i] = 0;
    }

    for (int i=0; i<natoms; i++){
        for (int j=0; j<neighmax; j++){
            neighlist[i][j] = 0;
        }
    }

    for (int i=0; i<natoms; i++){
        for (int j=0; j<neighmax; j++){
            neightags[i][j] = 0;
        }
    }

    for (int i=0; i<natoms*neighmax; i++){
        for (int j=0; j<3; j++){
            x[i][j] = 0.0;
        }
    }

    for (int i=0; i<natoms; i++){
        for (int j=0; j<3; j++){
            x[i][j] = input->x[i][j];
        }
    }

    //fh_debug = fopen("DEBUG_neighbor", "w");

}

Neighbor::~Neighbor() 
{

    memory->deallocate(numneigh);
    memory->deallocate(neighlist);
    memory->deallocate(neightags);
    memory->deallocate(x);

    //fclose(fh_debug);

};

void Neighbor::generate()
{


    double cutsq = input->rc*em3->input->rc;
    double xlength = input->box[0];
    double ylength = input->box[1];
    double zlength = input->box[2];
    double boxinvx = 1.0/xlength;
    double boxinvy = 1.0/ylength;
    double boxinvz = 1.0/zlength;

    int neighcount;
    int imagecounter = natoms; // used to store new images in the list of positions
    int nearintx, nearinty, nearintz;
    int nearintx2, nearinty2, nearintz2;
    double xi,yi,zi,xj,yj,zj;
    double xjp, yjp, zjp; // Image positions
    double xij,yij,zij, rsq, rij;

    double xtmp, ytmp, ztmp;

    /* If a half list is requested: */
    if (input->neighcount == 0){

        for (int i=0; i<natoms-1; i++){

            neighcount = 0;
            xtmp = x[i][0];
            ytmp = x[i][1];
            ztmp = x[i][2];

            for (int j=i+1; j < natoms; j++){
                if (j != i)
                {

                    xij = xtmp - x[j][0];
                    yij = ytmp - x[j][1];
                    zij = ztmp - x[j][2];

                    nearintx = std::round(xij*boxinvx);
                    nearinty = std::round(yij*boxinvy);
                    nearintz = std::round(zij*boxinvz);

                    xij = xij - xlength*nearintx;
                    yij = yij - ylength*nearinty;
                    zij = zij - zlength*nearintz;

                    rsq = xij*xij + yij*yij + zij*zij;

                    if (rsq < cutsq){

                        rij = sqrt(rsq);

                        /*
                        Need to add neighbor index to neighlist and position to x.
                        If neighbor is periodic, a new index is created along with a new position.
                        */

                        neightags[i][neighcount] = j;

                        if (nearintx == 0.0 && nearinty == 0.0 && nearintz == 0.0){

                            neighlist[i][neighcount] = j;

                        }
                        
                        else if (nearintx != 0.0 || nearinty != 0.0 || nearintz != 0.0){

                            // In this case, we must add this image of j to our list (x)

                            xjp = x[i][0] - xij;
                            yjp = x[i][1] - yij;
                            zjp = x[i][2] - zij;

                            x[imagecounter][0] = xjp;
                            x[imagecounter][1] = yjp;
                            x[imagecounter][2] = zjp;

                            neighlist[i][neighcount] = imagecounter;

                            imagecounter++;

                        } // if (nearintx == 0 && nearinty == 0 && nearintz == 0){

                        neighcount++;
                
                    } // if (rsq < cutsq){

                } // if (j!=i)

            } // for (int j=0; j<natoms; j++){

            numneigh[i] = neighcount;

        } // for (int i=0; i<natoms; i++){

    } // if (input->neighcount == 0){


    /* If a full list is requested: */
    else if (input->neighcount == 1){

        for (int i=0; i<natoms; i++){

            neighcount = 0;
            xtmp = x[i][0];
            ytmp = x[i][1];
            ztmp = x[i][2];

            for (int j=0; j < natoms; j++){
                if (j != i)
                {

                    xij = xtmp - x[j][0];
                    yij = ytmp - x[j][1];
                    zij = ztmp - x[j][2];

                    nearintx = std::round(xij*boxinvx);
                    nearinty = std::round(yij*boxinvy);
                    nearintz = std::round(zij*boxinvz);

                    xij = xij - xlength*nearintx;
                    yij = yij - ylength*nearinty;
                    zij = zij - zlength*nearintz;

                    rsq = xij*xij + yij*yij + zij*zij;

                    if (rsq < cutsq){

                        rij = sqrt(rsq);

                        /*
                        Need to add neighbor index to neighlist and position to x.
                        If neighbor is periodic, a new index is created along with a new position.
                        */

                        neightags[i][neighcount] = j;

                        if (nearintx == 0.0 && nearinty == 0.0 && nearintz == 0.0){
                            neighlist[i][neighcount] = j;
                        }
                        
                        if (nearintx != 0.0 || nearinty != 0.0 || nearintz != 0.0){

                            // In this case, we must add this image of j to our list (x)

                            xjp = x[i][0] - xij;
                            yjp = x[i][1] - yij;
                            zjp = x[i][2] - zij;

                            x[imagecounter][0] = xjp;
                            x[imagecounter][1] = yjp;
                            x[imagecounter][2] = zjp;

                            neighlist[i][neighcount] = imagecounter;

                            imagecounter++;

                        } // if (nearintx == 0 && nearinty == 0 && nearintz == 0){

                        neighcount++;
                
                    } // if (rsq < cutsq){

                } // if (j!=i)

            } // for (int j=0; j<natoms; j++){

            numneigh[i] = neighcount;

        } // for (int i=0; i<natoms; i++){

    } // if (input->neighcount == 1){



    // Find average number of neighbors
    neighavg =0.0;
    for (int i=0; i<natoms; i++){
        neighavg += numneigh[i];
    }
    neighavg = neighavg/natoms;


}
