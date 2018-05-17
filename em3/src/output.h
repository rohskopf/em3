/*
 output.h

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
  class Output: protected Pointers
  {
  public:
    Output(class EM3 *);
    ~Output();

    FILE * fh_xyz; // positions output in .xyz format
    FILE * fh_forces; // forces

    void write_xyz(); // function to write .xyz format position file
    void write_forces(); // function to write forces
    

  };
}

