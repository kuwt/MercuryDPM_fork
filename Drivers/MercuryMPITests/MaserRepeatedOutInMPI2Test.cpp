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


#include "DPMBase.h"
#include "Boundaries/SubcriticalMaserBoundaryTEST.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include "Species/LinearViscoelasticSpecies.h"


class MaserRepeatedOutInMPI2Test : public DPMBase
{
public:
    
    MaserRepeatedOutInMPI2Test()
    {
        setName("MaserRepeatedOutInMPI2Test");
    
        //set species properties: make the restitution coefficient very high so that it bounces back with enough energy
        LinearViscoelasticSpecies species;
        species.setDensity(6.0 / constants::pi);
        species.setCollisionTimeAndRestitutionCoefficient(0.005, 0.95, 1);
        speciesHandler.copyAndAddObject(species);
    
        //set time and file properties
        setTimeStep(species.getCollisionTime(1) / 50.0);
        setTimeMax(15);
        setSaveCount(5000);
        setParticlesWriteVTK(true);
        setWallsWriteVTK(FileType::ONE_FILE);
    
        //set domain size
        setMin({0,-1,-1});
        setMax({50, 1, 1});
    
        //for testing purposes, just set the gravity to 0.
        setGravity({0,0,0});
    }
    
    void setupInitialConditions() override {
        //Start with just one particle, which moves out, in, out of maser-boundary
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getLastObject());
        p0.setPosition({19,3,0});
        p0.setVelocity({1,1,0});
        p0.setRadius(0.5);
        particleHandler.copyAndAddObject(p0);
        
        //set the maser boundary
        SubcriticalMaserBoundaryTEST* b0 = boundaryHandler.copyAndAddObject(SubcriticalMaserBoundaryTEST());
        b0->set(Vec3D(1.0, 0.0, 0.0), 0.0, 20.0);
        b0->setActivationTime(0);
        
        PeriodicBoundary b;
        b.set(Vec3D(0,1,0), 0, 5);
        boundaryHandler.copyAndAddObject(b);
        
        //add two walls for the particle to bounce back from: an infinite wall on the right side, and a "double wall" in
        //the middle of the maser domain.
        InfiniteWall w;
        w.setSpecies(speciesHandler.getObject(0));
        w.setPosition(Vec3D(22, 0, 0));
        w.setNormal(Vec3D(1,0,0));
        wallHandler.copyAndAddObject(w);
        
        IntersectionOfWalls w1;
        w1.addObject(Vec3D(1,0,0), Vec3D(15, 0, 0));
        w1.addObject(Vec3D(0, 1, 0), Vec3D(0,-1,0));
        w1.addObject(Vec3D(-1, 0, 0), Vec3D(18, 0, 0));
        wallHandler.copyAndAddObject(w1);
    }
    
    void actionsAfterTimeStep() override
    {
        ///remove the walls to check if the new particles flow out & copy correctly
        if (getTime() > 7.5)
        {
            wallHandler.clear();
        }
    }
    
    void actionsAfterSolve() override
    {
        if (PROCESSOR_ID == 0)
        {
            logger(INFO, "I'm processor 0");
            logger(INFO, "Number of particles on processor 0: %", particleHandler.getNumberOfRealObjectsLocal());
        }
        if (PROCESSOR_ID == 1)
        {
            logger(INFO, "I'm processor 1");
            logger(INFO, "Number of particles on processor 1: %", particleHandler.getNumberOfRealObjectsLocal());
        }
    }
    
};

int main()
{
    MaserRepeatedOutInMPI2Test maserTest;
    maserTest.setNumberOfDomains({2,1,1});
    maserTest.solve();
}
