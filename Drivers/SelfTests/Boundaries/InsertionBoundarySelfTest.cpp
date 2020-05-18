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

#include <iostream>
#include "Mercury3D.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Walls/InfiniteWall.h"


class InsertionBoundarySelfTest : public Mercury3D
{
public:

    void setupInitialConditions() override {
        setName("InsertionBoundarySelfTest");
        setSystemDimensions(3);
        setGravity(Vec3D(0, 0, 0));
        setTimeStep(1e-4);
        dataFile.setSaveCount(10);
        setTimeMax(2e-2);
        setHGridMaxLevels(2);

        setMin(Vec3D(0, 0, 0));
        setMax(Vec3D(1, 1, 1));

        LinearViscoelasticSpecies species;
        species.setDensity(2000);
        species.setStiffness(10000);
        speciesHandler.copyAndAddObject(species);


        SphericalParticle insertionBoundaryParticle;
        insertionBoundaryParticle.setSpecies(speciesHandler.getObject(0));

        CubeInsertionBoundary insertionBoundary;
        insertionBoundary.set(&insertionBoundaryParticle,1,getMin(),getMax(),Vec3D(1,0,0),Vec3D(1,0,0),0.025,0.05);
        boundaryHandler.copyAndAddObject(insertionBoundary);

    }

    void printTime() const override
    {
        logger(INFO,"t=%, tMax=%, N=%", getTime(),getTimeMax(), particleHandler.getSize());
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    logger(INFO,"Simple box for creating particles");

    InsertionBoundarySelfTest insertionBoundary_problem;
    insertionBoundary_problem.solve();
}
