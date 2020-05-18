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
#include <iostream>
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include <cstdlib> //needed to be cygwin compatible (system not found)
#include "Boundaries/PeriodicBoundary.h"

class MpiPeriodicBoundaryUnitTest : public Mercury3D
{
public:

    void setupInitialConditions() override {
        //Add a particle
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(1.0);
    
        //Test 1: Particle moves through the corner of a periodic box
        p0.setPosition(Vec3D(7.0,7.0,99.5));
        p0.setVelocity(Vec3D(5.0, 0.0, 0.0));

        //Test 2: Particle moving through periodic boundaries while in an mpi zone
        //p0.setPosition(Vec3D(5.0,5.0,50.5));
        //p0.setVelocity(Vec3D(5.0, 0.0, 0.0));

        //Test 3: Particle moving through mpi boundaries while in a periodic zone
        //p0.setPosition(Vec3D(7.0,7.0,60.0));
        //p0.setVelocity(Vec3D(0.0, 0.0, -20.0));
        
        //Test 4: Diagonal test moving periodic and mpi boundary simultaniously
        //p0.setPosition(Vec3D(7.9,7.9,50.1));
        //p0.setVelocity(Vec3D(5.0, 5.0, -5.0));
        
        //Test 5: Tests what happens if a particle is located exactly on a periodic and mpi boundary
        //p0.setPosition(Vec3D(8.0,8.0,50.0));
        //p0.setVelocity(Vec3D(0.0, 0.0, 0.0)); 

        //Test 6: Perform an accuracy test that proves that my approach is justified
        //p0.setPosition(Vec3D(7.0, 7.0, 0.0000000000000001));
        //p0.setVelocity(Vec3D(0.0, 0.0, -0.00000000000001));

        //Test 7: Interaction to local domain test
        //p0.setPosition(Vec3D(4.9, 4.0, 97.99));
        //p0.setVelocity(Vec3D(-2.0, 0.0, 2.5));
        //particleHandler.copyAndAddObject(p0);
        //p0.setVelocity(Vec3D(2.0, 0.0, 2.5));
        //p0.setPosition(Vec3D(2.9001, 4.0, 97.99));


        particleHandler.copyAndAddObject(p0);

        //Add periodic boundary
        boundaryHandler.clear();
        PeriodicBoundary b0;
        b0.set(Vec3D(0.0, 0.0, 1.0), getZMin(), getZMax());
        boundaryHandler.copyAndAddObject(b0);
        b0.set(Vec3D(1.0, 0.0, 0.0), -8.0, 8.0);
        boundaryHandler.copyAndAddObject(b0);
        b0.set(Vec3D(0.0, 1.0, 0.0), -8.0, 8.0);
        boundaryHandler.copyAndAddObject(b0);
    }

    

};

int main(int argc, char *argv[])
{
    //Create the problem
    MpiPeriodicBoundaryUnitTest problem;
    problem.setName("PeriodicBoundaryUnitTestMPI");

    //Create a species
    LinearViscoelasticFrictionSpecies* species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    species->setDensity(1.0);

    //Set domain accordingly (domain boundaries are not walls!)
    problem.setXMin(-10.0);
    problem.setXMax(10.0);
    problem.setYMin(-10.0);
    problem.setYMax(10.0);
    problem.setZMin(0.0);
    problem.setZMax(100.0);

    //specify particle properties
    species->setDensity(6.0/constants::pi);

    //specify body forces
    problem.setGravity(Vec3D(0.0, 0.0, 0.0));

    //Set the number of domains for parallel decomposition
    problem.setNumberOfDomains({2,1,1});

    //specify contact properties
    //normal forces
    species->setStiffness(1e5);
    species->setDissipation(0.63);
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
    //Mdouble tc = species->getCollisionTime(MinParticleMass);
    Mdouble tc = 0.02;
    logger(INFO,"tc: %",tc);
    problem.setTimeStep(tc / 50.0);
    problem.setTimeMax(0.2);
    problem.setSaveCount(25); //used to be 500

    //Set output to paraview
    problem.setParticlesWriteVTK(true);
    
    problem.solve(argc, argv);

    return 0;
}

