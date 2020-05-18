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
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/IntersectionOfWalls.h"

/*!
 * This test tests what happens when two particles that have an interaction across a periodic boundary pass into another
 * MPI domain together. Currently there is an ugly solution to make sure we do not have a segfault in
 * Domain::processReceivedInteractionData, where the interaction data is deleted in this case. This is a known bug, but
 * it happens so rarely that we accept this temporary fix.
 */
class PeriodicBounaryEnteringMPIDomainTest : public Mercury3D
{
public:
    
    void setupInitialConditions() override
    {
        //Add a particle
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(0.5);
        
        //Test 1: Particle moves through the corner of a periodic box
        p0.setPosition(Vec3D(8, 9.51, 0));
        p0.setVelocity(Vec3D(1, 0, 0));
        particleHandler.copyAndAddObject(p0);
        p0.setPosition(Vec3D(8, 0.5, 0));
        p0.setVelocity(Vec3D(1, 0, 0));
        particleHandler.copyAndAddObject(p0);
        
        //Add periodic boundary
        boundaryHandler.clear();
        PeriodicBoundary b0;
        b0.set(Vec3D(0.0, 1.0, 0), 0, 10);
        boundaryHandler.copyAndAddObject(b0);
        
        IntersectionOfWalls wall;
        wall.addObject(Vec3D(0, 1, 0), Vec3D(0, 0.95, 0));
        wall.addObject(Vec3D(0, 0, 1), Vec3D(0,0,0));
        wall.addObject(Vec3D(0, -1, 0), Vec3D(0, 9.05, 0));
        wallHandler.copyAndAddObject(wall);
    }
};

int main(int argc, char* argv[])
{
    //Create the problem
    PeriodicBounaryEnteringMPIDomainTest problem;
    problem.setName("PeriodicBounaryEnteringMPIDomainTest");
    
    //Create a species
    LinearViscoelasticFrictionSpecies* species = problem.speciesHandler.copyAndAddObject(
            LinearViscoelasticFrictionSpecies());
    species->setDensity(1.0);
    
    //Set domain accordingly (domain boundaries are not walls!)
    problem.setXMin(0);
    problem.setXMax(20);
    problem.setYMin(-10.0);
    problem.setYMax(10.0);
    problem.setZMin(0.0);
    problem.setZMax(100.0);
    
    //specify particle properties
    species->setDensity(6.0 / constants::pi);
    
    //specify body forces
    problem.setGravity(Vec3D(0.0, 0.0, 0.0));
    
    //Set the number of domains for parallel decomposition
    problem.setNumberOfDomains({2, 1, 1});
    
    //specify contact properties
    //normal forces
    species->setStiffness(1e5);
    species->setDissipation(0.63);
    //tangential (sliding) forces
    species->setSlidingFrictionCoefficient(0);
    
    //set other simulation parameters
    Mdouble MinParticleMass = species->getDensity() * 4.0 / 3.0 * constants::pi;
    //Mdouble tc = species->getCollisionTime(MinParticleMass);
    Mdouble tc = 0.02;
    logger(INFO, "tc: %", tc);
    problem.setTimeStep(tc / 50.0);
    problem.setTimeMax(3);
    problem.setSaveCount(500);
    
    //Set output to paraview
    problem.setParticlesWriteVTK(true);
    
    problem.solve(argc, argv);
    
    return 0;
}

