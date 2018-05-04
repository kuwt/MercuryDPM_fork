//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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


#include <Walls/IntersectionOfWalls.h>
#include "DPMBase.h"
#include "Boundaries/SubcriticalMaserBoundaryTEST.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticSpecies.h"


class MaserRepeatedOutInUnitTest : public DPMBase
{
public:
    
    MaserRepeatedOutInUnitTest()
    {
        setName("MaserRepeatedOutInUnitTest");
    
        //set species properties: some standard values
        LinearViscoelasticSpecies species;
        species.setDensity(6.0 / constants::pi);
        species.setCollisionTimeAndRestitutionCoefficient(0.005, 0.95, 1);
        speciesHandler.copyAndAddObject(species);
    
        //set time and file properties
        setTimeStep(species.getCollisionTime(1) / 50.0);
        setTimeMax(7.5);
        setSaveCount(1000);
        setParticlesWriteVTK(true);
    
        //set domain size
        setMin({0,-1,-1});
        setMax({50, 1, 1});
    
        //for testing purposes, just set the gravity to 0.
        setGravity({0,0,0});
    }
    
    void setupInitialConditions()
    {
        //Check if particle is copied correctly when moving
        BaseParticle p0;
        p0.setSpecies(speciesHandler.getLastObject());
        p0.setPosition({19,0,0});
        p0.setVelocity({1,0,0});
        p0.setRadius(0.5);
        particleHandler.copyAndAddObject(p0);
        
        //set the maser boundary
        SubcriticalMaserBoundaryTEST* b0 = boundaryHandler.copyAndAddObject(SubcriticalMaserBoundaryTEST());
        b0->set(Vec3D(1.0, 0.0, 0.0), 0.0, 20.0);
        b0->setActivationTime(0);
        
        //PeriodicBoundary b;
        //b.set(Vec3D(0,1,0), -.1, 1);
        
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
        
        setParticlesWriteVTK(true);
        setWallsWriteVTK(FileType::ONE_FILE);
    }
    
    void actionsAfterTimeStep() override
    {
        if (getTime() > 15)
        {
            wallHandler.clear();
        }
    }
    
};

int main()
{
    MaserRepeatedOutInUnitTest maserTest;
    maserTest.setNumberOfDomains({2,1,1});
    maserTest.solve();
}
