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
#include "Boundaries/ChuteInsertionBoundary.h"
#include "Boundaries/HopperInsertionBoundary.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Walls/InfiniteWall.h"


class InsertionBoundarySelfTest : public Mercury3D
{
public:

    InsertionBoundarySelfTest()
    {
        /* JMFT: When using MPI, these need to be set here, or in main(),
         * *before* we call solve().  Hence, we cannot put them in
         * setupInitialConditions(). */
        setMin(Vec3D(0, 0, 0));
        setMax(Vec3D(1, 1, 1));

        /* JMFT: Species needs to be defined *before* calling solve() */
        LinearViscoelasticSpecies species;
        species.setDensity(2000);
        species.setStiffness(10000);
        speciesHandler.copyAndAddObject(species);
    }

    void setupInitialConditions() override {
        setName("InsertionBoundarySelfTestMPI");
        setSystemDimensions(3);
        setGravity(Vec3D(0.0,0.0,0.0));
        setTimeStep(1e-4);
        dataFile.setSaveCount(1);
        setTimeMax(1);
        setHGridMaxLevels(2);

        BaseParticle* insertionBoundaryParticle = new SphericalParticle;
        insertionBoundaryParticle->setSpecies(speciesHandler.getObject(0));

        //CubeInsertionBoundary
        auto insertionBoundary = new CubeInsertionBoundary; //delete is done in boundaryHandler
        boundaryHandler.addObject(insertionBoundary);
        insertionBoundary->set(insertionBoundaryParticle,1,getMin(),getMax(),Vec3D(1,0,0),Vec3D(1,0,0),0.025,0.05);

        /*
        //ChuteInsertionBoundary
        auto insertionBoundary = new ChuteInsertionBoundary;
        boundaryHandler.addObject(insertionBoundary);
        insertionBoundary->set(insertionBoundaryParticle,10,Vec3D(-0.5,-0.5,-0.5),Vec3D(0.5,0.5,0.5),0.1,0.2,0.1,1.0,0.2);
        */

        //HopperInsertionBoundary
        /*
        auto insertionBoundary = new HopperInsertionBoundary;
        boundaryHandler.addObject(insertionBoundary);
        double ExitHeight = 8.0;
        double ExitLength = ExitHeight;
        double hopperAngle = 45.0;
        double hopperLength = 3.0*ExitLength;
        double chuteAngle = 25;
        double lift = 0.0; 
        double fillPercentage = 50.0;
        insertionBoundary->set(insertionBoundaryParticle,100,0.0,0.0,0.1,0.1,
                                chuteAngle, 0.0, false, 2, hopperAngle, hopperLength, ExitLength,
                                ExitHeight, lift, fillPercentage);
        */

        delete insertionBoundaryParticle;
                
        InfiniteWall bottomWall;
        bottomWall.setSpecies(speciesHandler.getObject(0));
        const Vec3D normal = Vec3D(std::pow(1.0/3.0,0.5), std::pow(1.0/3.0,0.5), std::pow(1.0/3.0,0.5));
        bottomWall.set(normal, normal);
        wallHandler.copyAndAddObject(bottomWall);
    }

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    logger(INFO,"Simple box for creating particles");

    InsertionBoundarySelfTest insertionBoundary_problem;
    // insertionBoundary_problem.setParticlesWriteVTK(true);
    insertionBoundary_problem.setNumberOfDomains({2,1,1}); //For cube
    //insertionBoundary_problem.setNumberOfDomains({1,2,1});   //For chute
    //insertionBoundary_problem.setNumberOfDomains({1,2,1});   //For hopper

    insertionBoundary_problem.solve();
}
