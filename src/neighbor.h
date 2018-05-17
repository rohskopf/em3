/*
 neighbor.h

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
  class Neighbor: protected Pointers
  {
  public:
    Neighbor(class EM3 *);
    ~Neighbor();

    FILE * fh_debug; // Debug file handle

    // Declare member functions

    void generate();

    // Variables

    double **x; // positions of all atoms, including image neighbors
    int *numneigh; // Number of neighbors for every atom i
    int **neighlist; // contains neighbor indices j for every atom i (used to index positions from x)
    int **neightags; // tags (IDs) of neighbors which correspond to original atoms
    int natoms;
    int neighmax;
    double neighavg; // average number of neighbors
    

  };
}

