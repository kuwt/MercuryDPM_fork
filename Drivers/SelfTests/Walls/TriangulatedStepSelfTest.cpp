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
#include "Walls/TriangulatedWall.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Logger.h"

/*!
 * \brief Tests the implementation of TriangulatedWall.
 * \details
 *  <img src="Walls/triangulatedStepSelfTest.png" height="250px">
 */
class TriangulatedStepSelfTest : public Mercury3D
{
public:
    TriangulatedStepSelfTest ()
    {
        setName("TriangulatedStepSelfTest");
        setMin({0,0,0});
        setMax({3,1,2});

        setGravity(Vec3D(0,0,-10));
        setTimeStep(1e-4);
        setTimeMax(1.7);
        //setTimeMax(getTimeStep());
        setSaveCount(10);
        setXBallsAdditionalArguments("-v0 -solidf");
        //fStatFile.setFileType(FileType::NO_FILE);
        restartFile.setSaveCount(1e5);

        auto s = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        s->setDensity(1);
        s->setCollisionTimeAndRestitutionCoefficient(50.0*getTimeStep(),0.01,s->getMassFromRadius(0.5));
        s->setSlidingStiffness(0.00001*2./7.*s->getStiffness());
        s->setSlidingDissipation(0);
        s->setSlidingFrictionCoefficient(1e-2);
    }

    void setupInitialConditions() override
    {
        //introduce Wall
        TriangulatedWall w("TriangulatedStep.vtk",speciesHandler.getLastObject());
        w.write(std::cout);
        wallHandler.copyAndAddObject(w);

        //introduce single particle
        SphericalParticle p;
        p.setSpecies(speciesHandler.getLastObject());
        p.setRadius(0.5);
        p.setPosition({0.95,0.5,1.5});
        p.setVelocity({1,0,0});
        particleHandler.copyAndAddObject(p);
    }
};

/*!
 * \brief Tests the implementation of TriangulatedWall.
 * \details
 *  <img src="Walls/triangulatedStepSelfTest.png" height="250px">
 */
class TriangulatedStepWallSelfTest : public Mercury3D
{
public:
    TriangulatedStepWallSelfTest ()
    {
        setName("TriangulatedStepWallSelfTest");
        setMin({-.5,0,-.5});
        setMax({3.5,1,3.5});

        setGravity(Vec3D(0,0,0));
        setTimeStep(1e-12);
        setTimeMax(getTimeStep());
        setSaveCount(2);
        setXBallsAdditionalArguments("-v0 -solid -rred 15");
        fStatFile.setFileType(FileType::NO_FILE);
    }

    void setupInitialConditions() override
    {
        auto s = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        s->setDensity(1);
        s->setStiffness(1);
        auto ws = speciesHandler.copyAndAddObject(s);
        s->setStiffness(0);

        //introduce Wall
        TriangulatedWall w("TriangulatedStepWall.vtk",ws);
        w.write(std::cout);
        wallHandler.copyAndAddObject(w);

        //introduce single particle
        SphericalParticle p;
        p.setSpecies(s);
        p.setRadius(0.4);
        Mdouble h = 0.05*2;

        for (Mdouble x = getXMin(); x<getXMax(); x+=h){
            for (Mdouble z = getZMin(); z<getZMax(); z+=h)
            {
                Mdouble distance;
                Vec3D normal;
                p.setPosition({x,.5,z});
                if (w.getDistanceAndNormal(p,distance,normal))
                    particleHandler.copyAndAddObject(p);
            }
        }
        logger(INFO,"Inserted % particles",particleHandler.getNumberOfObjects());
    }
};


int main()
{
    //input file without wall at the end.
    helpers::writeToFile("TriangulatedStep.vtk","# vtk DataFile Version 2.0\n"
     "Step\n"
     "ASCII\n"
     "DATASET UNSTRUCTURED_GRID\n"
     "POINTS 8 double\n"
     "0 0 1\n"
     "1 0 1\n"
     "2 0 0.9\n"
     "3 0 0.9\n"
     "0 1 1\n"
     "1 1 1\n"
     "2 1 0.9\n"
     "3 1 0.9\n"
     "\n"
     "CELLS 6 24\n"
     "3 0 1 4\n"
     "3 4 1 5\n"
     "3 5 1 2\n"
     "3 5 2 6\n"
     "3 6 2 3\n"
     "3 6 3 7\n"
     "\n"
     "CELL_TYPES 6\n"
     "5\n"
     "5\n"
     "5\n"
     "5\n"
     "5\n"
     "5");

    //input file with a wall at the end
    helpers::writeToFile("TriangulatedStepWall.vtk","# vtk DataFile Version 2.0\n"
     "Step\n"
     "ASCII\n"
     "DATASET UNSTRUCTURED_GRID\n"
     "POINTS 10 double\n"
     "0 0 2\n"
     "1 0 2\n"
     "2 0 1\n"
     "3 0 1\n"
     "3 0 3\n"
     "0 1 2\n"
     "1 1 2\n"
     "2 1 1\n"
     "3 1 1\n"
     "3 1 3\n"
     "\n"
     "CELLS 8 32\n"
     "3 0 1 5\n"
     "3 5 1 6\n"
     "3 6 1 2\n"
     "3 6 2 7\n"
     "3 7 2 3\n"
     "3 7 3 8\n"
     "3 8 3 4\n"
     "3 8 4 9\n"
     "\n"
     "CELL_TYPES 8\n"
     "5\n"
     "5\n"
     "5\n"
     "5\n"
     "5\n"
     "5\n"
     "5\n"
     "5");

//    //input file with a wall at the end
//    helpers::writeToFile("TriangulatedStepWall.vtk","# vtk DataFile Version 2.0\n"
//     "Step\n"
//     "ASCII\n"
//     "DATASET UNSTRUCTURED_GRID\n"
//     "POINTS 6 double\n"
//     "0 0 .5\n"
//     ".5 0 .5\n"
//     "10 0 10\n"
//     "0 1 .5\n"
//     ".5 1 .5\n"
//     "10 1 10\n"
//     "\n"
//     "CELLS 4 16\n"
//     "3 0 1 3\n"
//     "3 3 1 4\n"
//     "3 4 1 2\n"
//     "3 4 2 5\n"
//     "\n"
//     "CELL_TYPES 4\n"
//     "5\n"
//     "5\n"
//     "5\n"
//     "5");

    //plot location of contact point
    helpers::writeToFile("TriangulatedStepSelfTest.gnu","p 'TriangulatedStepSelfTest.fstat' u 4:6 w lp");

    //TriangulatedStepSelfTest t;
    //t.solve();
    TriangulatedStepWallSelfTest t;
    t.solve();
    return 0;
}
