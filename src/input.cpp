/*
 input.cpp

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

#include "input.h"
#include "memory.h"
#include "compute.h"

using namespace std;

using namespace EM3_NS;

Input::Input(EM3 *em3) : Pointers(em3) {
    //fh_debug = fopen("DEBUG_input", "w");
}

Input::~Input() 
{

    em3->memory->deallocate(x);
    em3->memory->deallocate(v);
    em3->memory->deallocate(x0);

    //fclose(fh_debug);

};

void Input::readinput()
{
    /* Read INPUT file */

    string line;

    // Declare scalar inputs
    double value;

    // Open INPUT file
    ifstream INPUT("INPUT");
    // Ignore the first line
    getline(INPUT, line); 
    string characters;
    // Get input variables
    for (int i=1; i<=13; i++)
    {
        getline(INPUT, line);
        switch (i)
        {
            case 1:{ 
                stringstream ss(line);
                ss >> characters >> nsteps;
            }
            case 2:{ 
                stringstream ss(line);
                ss >> characters >> nequil;
            }
            case 3:{ 
                stringstream ss(line);
                ss >> characters >> nout;
            }
            case 4:{ 
                stringstream ss(line);
                ss >> characters >> neighmax;
            }
            case 5:{ 
                stringstream ss(line);
                ss >> characters >> rc;
            }
            case 6:{
                stringstream ss(line);
                ss >> characters >> neighcount;
            }
            case 7:{
                stringstream ss(line);
                ss >> characters >> newton;
            }
            case 8:{
                stringstream ss(line);
                ss >> characters >> offset;
            }
            case 9:{ 
                stringstream ss(line);
                ss >> characters >> m;
            }
            case 10:{ 
                stringstream ss(line);
                ss >> characters >> temp;
            }
            case 11:{ 
                stringstream ss(line);
                ss >> characters >> dt;
            }
            case 12:{
                stringstream ss(line);
                ss >> characters >> epsilon;
            }
            case 13:{
                stringstream ss(line);
                ss >> characters >> sigma;
            }

  
        } // switch (i)

    } // for (int i=1..)

    INPUT.close();

}

void Input::readconfig()
{
    /* Read CONFIG file */

    ifstream config("CONFIG");
    string line;

    getline(config, line);
    stringstream ss(line);
    ss >> natoms;

    getline(config, line);
    stringstream ss2(line);
    ss2 >> box[0] >> box[1] >> box[2];
    //fprintf(fh_debug, "Box: %f %f %f\n", box[0], box[1], box[2]);

    printf("Reading %d atoms in box (%.2f, %.2f, %.2f)\n", natoms, box[0],box[1],box[2]);

    em3->memory->allocate(x, natoms*neighmax, 3); 
    em3->memory->allocate(x0, natoms, 3); 
    for (int i=0; i<natoms*neighmax; i++){
        for (int j=0; j<3; j++){
            x[i][j] = 0.0;
        }
    }
    for (int i=0; i<natoms; i++){
        for (int j=0; j<3; j++){
            x0[i][j] = 0.0;
        }
    }

    int type,tag;
    for (int i=0; i<natoms; i++){

        getline(config, line);
        stringstream ss3(line);
        ss3 >> type >> tag >> x[i][0] >> x[i][1] >> x[i][2];
        
    }
    config.close();

    // Store original positions

    for (int i=0; i<natoms; i++){
        for (int j=0; j<3; j++){
            x0[i][j] = x[i][j];
        }
    }

}

void Input::initialize()
{

    double k = 1.38064852e-23; // Boltzmann constant [J/K]

    // Convert inputs to SI units before converting to LJ units

    dt = dt*1e-15;
    

    // Convert all quantities to LJ units

    tau = sigma*sqrt(m/epsilon);
    velocity = sqrt(epsilon/m);
    force = epsilon/sigma;
    pressure = epsilon/(sigma*sigma*sigma);
    temperature = epsilon/k;
    printf("Units measured in:\n");
    printf("  length: %e m\n", sigma);
    printf("  energy: %e J\n", epsilon);
    printf("  mass: %e kg\n", m);
    printf("  time: %e s\n", tau );
    printf("  velocity: %e m/s\n", velocity );
    printf("  force: %e N\n", force );
    printf("  pressure: %e N/m^2\n", pressure );
    printf("  temperature: %f K\n", temperature );

    volume = box[0]*box[1]*box[2];

    temp = temp/temperature;

    dt = dt/tau;
    
    printf("Tau: %e\n", tau);
    printf("Cutoff [sigma]: %f\n", rc);
    printf("Temperature [k/epsilon]: %f\n", temp);
    printf("Box [sigma]: %f %f %f\n", box[0],box[1],box[2]);
    printf("Timestep [tau]: %e\n", dt);

    // Initialize MB velocities

    em3->memory->allocate(v, natoms, 3); 
    for (int i=0; i<natoms; i++){
        for (int j=0; j<3; j++){
            v[i][j] = 0.0;
        }
    }

    // random device class instance, source of 'true' randomness for initializing random seed
    std::random_device rd; 

    // Mersenne twister PRNG, initialized with seed from previous random device instance
    std::mt19937 gen(rd()); 
    
    // instance of class std::normal_distribution with specific mean and stddev
    std::normal_distribution<float> d(0.0, 1.0); 

    double normal;
    for (int i=0; i<natoms; i++){
        for (int j=0; j<3; j++){

            // get random number with normal distribution using gen as random source
            normal = d(gen); 
            //v[i][j] = sqrt(temp)*normal;
            v[i][j] = normal;
        }
    }

    // Calculate velocity center of mass

    double vx_com = 0.0;
    double vy_com = 0.0;
    double vz_com = 0.0;
    for (int i=0; i<natoms; i++){
        vx_com += v[i][0];
        vy_com += v[i][1];
        vz_com += v[i][2];
    }
    vx_com = vx_com/(natoms);
    vy_com = vy_com/(natoms);
    vz_com = vz_com/(natoms);

    // Subtract velocity COM from all components
    for (int i=0; i<natoms; i++){

        v[i][0] -= vx_com;
        v[i][1] -= vy_com;
        v[i][2] -= vz_com;

    }

    // Calculate velocity center of mass for sanity check
    /*
    for (int i=0; i<natoms; i++){
        vx_com += v[i][0];
        vy_com += v[i][1];
        vz_com += v[i][2];
    }
    vx_com = vx_com/(natoms);
    vy_com = vy_com/(natoms);
    vz_com = vz_com/(natoms);
    fprintf(fh_debug, "%f %f %f\n", vx_com, vy_com, vz_com);
    */

    // Scale each velocity component
    for (int i=0; i<natoms; i++){
        for (int j=0; j<3; j++){
            v[i][j] = sqrt(temp)*v[i][j];
        }
    }

    
    // Kinetic energy is half times the sum of squared speeds in LJ units

    double ke=0.0;
    for (int i=0; i<natoms; i++){

        ke += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2])/2.0;

    }

    double temp_actual = ke*(2.0/3.0)/natoms;

    // Recale the velocities
    for (int i=0; i<natoms; i++){
        v[i][0] = v[i][0]*sqrt(temp/temp_actual);
        v[i][1] = v[i][1]*sqrt(temp/temp_actual);
        v[i][2] = v[i][2]*sqrt(temp/temp_actual);
    }

}


