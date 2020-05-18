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
#include "Boundaries/SubcriticalMaserBoundaryTEST.h"
#include "Walls/IntersectionOfWalls.h"


class MpiMaserChuteTest : public Mercury3D
{
public:

    void setupInitialConditions() override {
        //Add a particle
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));


        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 8; j++)
            {
                for (int k = 0; k < 10; k++)
                {
                    p0.setPosition(Vec3D(1.0+i*2.0, 1.0+j*2.0,1.0+k*2.0));
                    p0.setRadius(random.getRandomNumber(0.9, 1.0));
                    p0.setVelocity(Vec3D(3.0, random.getRandomNumber(-0.1,0.1), 0.0));
                    particleHandler.copyAndAddObject(p0);
                }
            }
        }

/*
        p0.setRadius(1.0);
        p0.setPosition(Vec3D(1.0,1.0,1.0));
        particleHandler.copyAndAddObject(p0);
*/

        //Add maser boundary
        boundaryHandler.clear();
        SubcriticalMaserBoundaryTEST m0;
        m0.set(Vec3D(1.0, 0.0, 0.0), 0.0, 10.0);
        m0.setActivationTime(10.0);
        boundaryHandler.copyAndAddObject(m0);

        //Add periodic boundary at side
        PeriodicBoundary b0;
        b0.set(Vec3D(0.0, 1.0, 0.0), 0.0, 20.0);
        boundaryHandler.copyAndAddObject(b0);

        //Add bottom wall
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0, 0, -1), Vec3D(0.0, 0, 0));
        wallHandler.copyAndAddObject(w0);
    }
};

int main(int argc, char *argv[])
{
    //Create the problem
    MpiMaserChuteTest problem;
    problem.setName("MaserChuteTestMPI");

    //Create a species
    LinearViscoelasticFrictionSpecies* species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    species->setDensity(1.0);

    //Set domain accordingly (domain boundaries are not walls!)
    problem.setXMin(-10.0);
    problem.setXMax(150.0);
    problem.setYMin(-10.0);
    problem.setYMax(100.0);
    problem.setZMin(-10.0);
    problem.setZMax(100.0);

    //specify particle properties
    species->setDensity(6.0/constants::pi);

    //specify body forces
    problem.setGravity(Vec3D(1.0*0.37460659341, 0.0, -1.0*0.92718385456));

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
    problem.setTimeMax(52.4);
    problem.setSaveCount(250); //used to be 500

    //Set output to paraview
    problem.setParticlesWriteVTK(true);
    
    problem.solve(argc, argv);

    return 0;
}

