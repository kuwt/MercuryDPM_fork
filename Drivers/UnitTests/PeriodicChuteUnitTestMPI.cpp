//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include <iostream>
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include <cstdlib> //needed to be cygwin compatible (system not found)
#include "Boundaries/PeriodicBoundary.h"
#include <chrono>

class MpiPeriodicBoundaryUnitTest : public Mercury3D
{
public:

    void setupInitialConditions() override {
        //Define wall
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0, 0, -1),Vec3D(0,0,0));
        wallHandler.copyAndAddObject(w0);

        //Add periodic boundary
        boundaryHandler.clear();
        PeriodicBoundary b0;
        b0.set(Vec3D(0.0, 1.0, 0.0), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(b0);
        b0.set(Vec3D(1.0, 0.0, 0.0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(b0);
    
        //Add a particles
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        p0.setRadius(1.0);
        double np = 1000;
        int fail = 0;
        Vec3D position;
        for (int i = 0; i < np; i++) 
        {
            position = Vec3D(   random.getRandomNumber(getXMin() + 2 , getXMax() - 2),
                                random.getRandomNumber(getYMin() + 2, getYMax() - 2),
                                random.getRandomNumber(1.9,2.1));

            p0.setPosition(position);
            p0.setVelocity({0.0, 0.0, 0.0});
            if (checkParticleForInteraction(p0))
            {
                particleHandler.copyAndAddObject(p0);
                fail = 0;
            }
            else
            {
                fail++;
            }
            
            if (fail > 5000)
            {
                break;
            }
        }

        logger(INFO,"Number of particles: %",particleHandler.getSize());
  }

    void actionsAfterTimeStep() override {
    }       
 
};

int main(int argc, char *argv[])
{

    // Start measuring elapsed time
    std::chrono::time_point<std::chrono::system_clock> startClock, endClock;
    startClock = std::chrono::system_clock::now();

    //Create the problem
    MpiPeriodicBoundaryUnitTest problem;
    problem.setName("PeriodicChuteUnitTestMPI");

    //Create a species
    LinearViscoelasticFrictionSpecies* species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    species->setDensity(1.0);

    //Set domain accordingly (domain boundaries are not walls!)
    problem.setXMin(0.0);
    problem.setXMax(80.0);
    problem.setYMin(0.0);
    problem.setYMax(80.0);
    problem.setZMin(0.0);
    problem.setZMax(40.0);

    //specify particle properties
    species->setDensity(6.0/constants::pi);

    //specify body forces
    double angle = 27.0;
    double g = 1.0;
    double val = g*std::sin(angle/360.0*2.0*constants::pi);
    problem.setGravity(Vec3D(val, -val, -g*std::cos(angle/360.0*2.0*constants::pi)));

    //Set the number of domains for parallel decomposition
    problem.setNumberOfDomains({2,2,1});

    //specify contact properties
    //normal forces
    species->setStiffness(1e5);
    species->setDissipation(0.1);
    //tangential (sliding) forces
    species->setSlidingFrictionCoefficient(0.5);
    species->setSlidingStiffness(1.2e4);
    species->setSlidingDissipation(0.16);
    //tangential (rolling) torques
    species->setRollingFrictionCoefficient(0.2);
    species->setRollingStiffness(1.2e4);
    species->setRollingDissipation(6.3e-2);
    //normal (torsion/spin) torques
    species->setTorsionFrictionCoefficient(0.1);
    species->setTorsionStiffness(1.2e4);
    species->setSlidingDissipation(6.3e-2);
 
    //set other simulation parameters
    Mdouble MinParticleMass = species->getDensity()*4.0 / 3.0 * constants::pi ;    
    Mdouble tc = species->getCollisionTime(MinParticleMass);
    problem.setTimeStep(tc / 50.0);
    //problem.setTimeMax(100.0);//run until 3.0 to see full simulation
    problem.setTimeMax(40.0);//run until 3.0 to see full simulation
    problem.setSaveCount(200); //used to be 500

    //Set output to paraview
    problem.setParticlesWriteVTK(true);
    
    problem.solve(argc, argv);

    // Measure elapsed time
    endClock = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = endClock - startClock;
    logger(INFO, "Elapsed time for solving the PDE: % s", elapsed_seconds.count());


    return 0;
}

