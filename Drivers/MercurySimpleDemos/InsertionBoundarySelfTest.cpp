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

#include <iostream>
#include "Mercury3D.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Particles/BaseParticle.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Walls/InfiniteWall.h"


class InsertionBoundarySelfTest : public Mercury3D
{
public:

    void setupInitialConditions()
    {
        setName("InsertionBoundarySelfTest");
        setSystemDimensions(3);
        setGravity(Vec3D(0.0,0.0,0.0));
        setTimeStep(1e-4);
        dataFile.setSaveCount(50);
        setTimeMax(0.5);
        setHGridMaxLevels(2);

        setXMin(0.0);
        setYMin(0.0);
        setZMin(0.0);
        setXMax(1.0);
        setYMax(0.01);
        setZMax(1.0);
        
        LinearViscoelasticSpecies species;
        species.setDensity(2000);
        species.setStiffness(10000);
        speciesHandler.copyAndAddObject(species);


        BaseParticle insertionBoundaryParticle;
        insertionBoundaryParticle.setSpecies(speciesHandler.getObject(0));

        //CubeInsertionBoundary::set(BaseParticle* particleToCopy, int maxFailed, Vec3D posMin, Vec3D posMax, Vec3D velMin, Vec3D velMax, double radMin, double radMax)
        CubeInsertionBoundary insertionBoundary;
        insertionBoundary.set(&insertionBoundaryParticle,1,getMin(),getMax(),Vec3D(0,0,0),Vec3D(0,0,0),0.025,0.05);
        //insertionBoundary.checkBoundaryBeforeTimeStep(this);
        boundaryHandler.copyAndAddObject(insertionBoundary);

//        InfiniteWall bottomWall;
//        bottomWall.setSpecies(speciesHandler.getObject(0));
//        bottomWall.set(Vec3D(1,1,1), 0.2*getMin()+0.8*getMax());
//        wallHandler.copyAndAddObject(bottomWall);

        PeriodicBoundary periodicBoundary;
        periodicBoundary.set(Vec3D(1,0,0),0,1);
        boundaryHandler.copyAndAddObject(periodicBoundary);
        periodicBoundary.set(Vec3D(0,0,1),0,1);
        boundaryHandler.copyAndAddObject(periodicBoundary);
    }

    void printTime() const override
    {
        logger(INFO,"t=%, tMax=%, N=%", getTime(),getTimeMax(), particleHandler.getSize());
    }


    CubeInsertionBoundary* insertionBoundary;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    logger(INFO,"Simple box for creating particles");

    InsertionBoundarySelfTest insertionBoundary_problem;
    insertionBoundary_problem.solve();
}
