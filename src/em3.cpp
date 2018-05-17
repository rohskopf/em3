/*
 em3.cpp

 Copyright (c) 2018 Andrew Rohskopf

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include "mpi.h"

#include "memory.h"
#include "timer.h"
#include "input.h"
#include "neighbor.h"
#include "potential.h"
#include "update.h"
#include "output.h"
#include "compute.h"

using namespace std;

using namespace EM3_NS;

EM3::EM3(int narg, char **arg)
{

    //fh_debug = fopen("DEBUG_em3", "w");

    /************************** Set up MPI settings **************************/

    int color,key,global,local;
    MPI_Comm comm;

    // Split the communicators so that multiple instances can be run
    MPI_Comm_rank(MPI_COMM_WORLD, &global);
    color = global / 1; // Change "1" to 2 in order to use 2 procs per instance, etc..
    key = global; 
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &comm);
    MPI_Comm_rank(comm,&local);

    //  Get the number of processes.
    procs = MPI::COMM_WORLD.Get_size ( ); //Get_size gets number of processes (np) in communicator group
    //  Get the individual process ID.
    rank = MPI::COMM_WORLD.Get_rank ( ); // Get_rank gets the rank of the calling process in the communicator

    /************************** Initial Screen Output **************************/
    if (rank == 0)
    {
        std::cout << " +-----------------------------------------------------------------+" << std::endl;
        std::cout << " +                            EM3 0.0                              +" << std::endl;
        std::cout << " +-----------------------------------------------------------------+" << std::endl;
        std::cout << " Running on " << procs << " procs" << std::endl;
    }
  
    timer = new Timer(this);

    if (rank == 0) std::cout << " Job started at " << timer->DateAndTime() << std::endl;

    /************************** Proceed with Program **************************/

    // Initialize system

    input = new Input(this);
    
    initialize();

    int natoms = input->natoms;
    int neighmax = input->neighmax;

    // Declare output pointer

    output = new Output(this);

    // Dynamically allocated pointers

    create();

    // Generate neighbor list

    neighbor->generate();

    // Calculate potential and forces

    potential->calculate();

    // Perform MD

    printf("Running for %d steps.\n", input->nsteps);
    printf("Step PE KE E T P MSD AverageNeigh\n");
    compute->compute_ke();
    compute->compute_pressure();
    compute->compute_msd();
    compute->compute_com();

    printf("%d %f %f %f %f %f %f %f\n", 0, potential->pe/natoms, compute->ke/natoms, compute->etot/natoms, compute->temp,\
            compute->pressure, compute->msd, neighbor->neighavg);
    output->write_xyz();
    //output->write_forces();

    for (t=1; t < input->nsteps+1;t++){

        // Compute KE for temperature rescaling

        compute->compute_ke();

        // Update the timestep

        update->integrate();

        // Write output files

        if (t % input->nout == 0){

            compute->compute_pressure();
            compute->compute_msd();
            compute->compute_com();

            printf("%d %f %f %f %f %f %f %f\n", t, potential->pe/natoms, compute->ke/natoms, compute->etot/natoms, compute->temp,\
                    compute->pressure, compute->msd, neighbor->neighavg);
            output->write_xyz();
            //output->write_forces();
        }
    }
    

    // Delete dynamically allocated pointers

    finalize();

    if (rank == 0) std::cout << std::endl << " Job finished at " 
        << timer->DateAndTime() << std::endl;
    if (rank == 0) timer->print_elapsed();
}

void EM3::create()
{
    memory = new Memory(this);
    neighbor = new Neighbor(this);
    potential = new Potential(this);
    update = new Update(this);
    compute = new Compute(this);

}

void EM3::initialize()
{
    input->readinput();
    input->readconfig();
    input->initialize();
}

EM3::~EM3()
{
    delete timer;
    delete input;
    delete output;

    //fclose(fh_debug);

}

void EM3::finalize()
{

    delete memory;
    delete neighbor;
    delete potential;
    delete update;
    delete compute;


}

