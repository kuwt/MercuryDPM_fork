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
#include "Boundaries/PolydisperseInsertionBoundary.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Walls/InfiniteWall.h"

class PolydisperseInsertionBoundarySelfTest : public Mercury3D
{
public:

    void setupInitialConditions() override {
        setName("PolydisperseInsertionBoundarySelfTest");
        setSystemDimensions(2);
        setGravity(Vec3D(0.0,0.0,0.0));
        setTimeStep(1e-3);
        dataFile.setSaveCount(10);
        setTimeMax(5e-1);
        setHGridMaxLevels(2);

        setXMin(0.0);
        setYMin(0.0);
        setZMin(0.0);
        setXMax(1.0);
        setYMax(1.0);
        setZMax(1.0);
        
        double radA = 0.02;
        double radB = 0.04;

        double tc = 2e-2;
        double rest = 0.6;

        LinearViscoelasticSpecies specA;
        specA.setDensity(8);
        specA.setCollisionTimeAndRestitutionCoefficient(tc, rest, 
                4./3. * constants::pi * pow(radA, 3));
        speciesHandler.copyAndAddObject(specA);

        LinearViscoelasticSpecies specB;
        specB.setDensity(1);
        specB.setCollisionTimeAndRestitutionCoefficient(tc, rest, 
                4./3. * constants::pi * pow(radB, 3));
        speciesHandler.copyAndAddObject(specB);

        auto genA = new SphericalParticle();
        genA->setSpecies(speciesHandler.getObject(0));
        genA->setRadius(radA);
        auto genB = new SphericalParticle();
        genB->setSpecies(speciesHandler.getObject(1));
        genB->setRadius(radB);

        auto insb = boundaryHandler.copyAndAddObject(new PolydisperseInsertionBoundary());
        insb->setGeometry(1,
                Vec3D( getXMin(), getYMin(), getZMin() ),
                Vec3D( getXMax(), getYMax(), 0 ),
                Vec3D(-1,-1,0),Vec3D(1,1,0));
        insb->addGenerandum(genA, 1, 0.2);
        insb->addGenerandum(genB, 1, 0.4);

        InfiniteWall wall;
        wall.setSpecies(speciesHandler.getObject(0));
        /*
        wall.set(Vec3D(0,0,-1), getMin());
        wallHandler.copyAndAddObject(wall);
        */
        wall.set(Vec3D(-1,0,0), getMin());
        wallHandler.copyAndAddObject(wall);
        wall.set(Vec3D(+1,0,0), getMax());
        wallHandler.copyAndAddObject(wall);
        wall.set(Vec3D(0,-1,0), getMin());
        wallHandler.copyAndAddObject(wall);
        wall.set(Vec3D(0,+1,0), getMax());
        wallHandler.copyAndAddObject(wall);
    }

    void printTime() const override
    {
        logger(INFO,"t=%, tMax=%, N=%", 
                getTime(), getTimeMax(), particleHandler.getSize()
                );
    }

    void actionsAfterTimeStep() override {
        if (particleHandler.getSize() > 300 && boundaryHandler.getSize() > 0)
        {
            logger(INFO, "particleHandler.getSize() = %, boundaryHandler.getSize() = %",
                    particleHandler.getSize(), boundaryHandler.getSize());
            boundaryHandler.clear();
        }
    }
};

int main(int argc , char *argv[] )
{
    logger(INFO,"PolydisperseInsertionBoundarySelfTest test.");

    PolydisperseInsertionBoundarySelfTest problem;
    problem.solve(argc, argv);
}
