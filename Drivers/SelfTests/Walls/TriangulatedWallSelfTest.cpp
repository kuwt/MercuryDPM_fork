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
#include "Species/LinearViscoelasticSpecies.h"
#include "Logger.h"

/*!
 * \brief Tests the implementation of TriangulatedWall.
 * \details
 *  <img src="Walls/triangulatedWallsSelfTest.png" height="250px">
 */
class TriangulatedWallSelfTest : public Mercury3D
{
public:
    TriangulatedWallSelfTest ()
    {
        setName("TriangulatedWallSelfTest");
        setMin({-.2,-.2,-.2});
        setMax({1.2,1.2,1.2});

        setGravity(10*Vec3D(0,0,-1));
        setTimeStep(1e-4);
        setTimeMax(1e4*getTimeStep());
        //setTimeMax(getTimeStep());
        setSaveCount(100);
        setXBallsAdditionalArguments("-v0 -solidf");
        fStatFile.setFileType(FileType::NO_FILE);
        restartFile.setSaveCount(1e5);

        auto s = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        s->setDensity(1);
        s->setCollisionTimeAndRestitutionCoefficient(20.0*getTimeStep(),0.1,s->getMassFromRadius(radius));
    }

    void setupInitialConditions() override
    {
        InfiniteWall i;
        i.setSpecies(speciesHandler.getLastObject());
        i.set({0,0,-1},getMin());
        wallHandler.copyAndAddObject(i);
        i.set({0,-1,0},getMin());
        wallHandler.copyAndAddObject(i);
        i.set({-1,0,0},getMin());
        wallHandler.copyAndAddObject(i);
        i.set({0,1,0},getMax());
        wallHandler.copyAndAddObject(i);
        i.set({1,0,0},getMax());
        wallHandler.copyAndAddObject(i);

        TriangulatedWall w("TriangulatedWall.vtk",speciesHandler.getLastObject());
        //TriangulatedWall w("TriangulatedWallInverted.vtk",speciesHandler.getLastObject());
        w.write(std::cout);
        wallHandler.copyAndAddObject(w);

        //introduce particles
        SphericalParticle p;
        p.setSpecies(speciesHandler.getLastObject());

        //orderly insertion
//        Mdouble dx = 1.0/(n-1);
//        p.setRadius(0.49*dx);
//        for (unsigned i=0; i<n; i++)
//        for (unsigned j=0; j<n; j++)
//        for (unsigned k=0; k<n; k++)
//        {
//            p.setPosition({i*dx-1e-8,j*dx-1e-8,k*dx-1e-8});
//            particleHandler.copyAndAddObject(p);
//        }

        //random insertion
//        p.setRadius(0.1);
//        for (unsigned i=0; i<100; i++)
//        {
//            Vec3D r;
//            r.X = random.getRandomNumber(0,1);
//            r.Y = random.getRandomNumber(1-r.X,1);
//            r.Z = random.getRandomNumber(1-r.X-r.Y,2);
//            p.setPosition(r);
//            if (checkParticleForInteraction(p))
//                particleHandler.copyAndAddObject(p);
//        }

        //random insertion
        for (unsigned i=0; i<300; i++)
        {
            Vec3D r;
            r.X = random.getRandomNumber(getXMin(),getXMax());
            r.Y = random.getRandomNumber(getYMin(),getYMax());
            r.Z = random.getRandomNumber(getZMax(),getZMax()+1.0);
            p.setPosition(r);
            p.setRadius(random.getRandomNumber(radius,1.1*radius));
            if (checkParticleForInteraction(p))
                particleHandler.copyAndAddObject(p);
        }

        //single-particle insertion
//        p.setRadius(0.1);
//        p.setPosition({-.01,-.01,-.01}); //should create three contacts
//        particleHandler.copyAndAddObject(p);

        logger(INFO,"Inserted % particles",particleHandler.getNumberOfObjects());
    }

    Mdouble radius = 0.08;
};

int main()
{
    //write input files
    helpers::writeToFile("TriangulatedWallSimple.vtk","# vtk DataFile Version 2.0\n"
     "Defines a domain outside of a tetrahedron\n"
     "ASCII\n"
     "DATASET UNSTRUCTURED_GRID\n"
     "POINTS 4 double\n"
     "0 0 0\n"
     "1 0 0\n"
     "0 1 0\n"
     "0 0 1\n"
     "\n"
     "CELLS 4 16\n"
     "3 0 1 2\n"
     "3 0 2 3\n"
     "3 0 3 1\n"
     "3 1 3 2\n"
     "\n"
     "CELL_TYPES 4\n"
     "5\n"
     "5\n"
     "5\n"
     "5");

    helpers::writeToFile("TriangulatedWallInverted.vtk","# vtk DataFile Version 2.0\n"
     "Defines a domain inside of a tetrahedron\n"
     "ASCII\n"
     "DATASET UNSTRUCTURED_GRID\n"
     "POINTS 4 double\n"
     "0 0 0\n"
     "1 0 0\n"
     "0 1 0\n"
     "0 0 1\n"
     "\n"
     "CELLS 4 16\n"
     "3 0 2 1\n"
     "3 0 3 2\n"
     "3 0 1 3\n"
     "3 1 2 3\n"
     "\n"
     "CELL_TYPES 4\n"
     "5\n"
     "5\n"
     "5\n"
     "5");

    helpers::writeToFile("TriangulatedWall.vtk","# vtk DataFile Version 2.0\n"
     "Cube Normals (cross of first two vec's) into the wall\n"
     "ASCII\n"
     "DATASET UNSTRUCTURED_GRID\n"
     "POINTS 8 double\n"
     "0 0 0\n"
     "1 0 0\n"
     "1 1 0\n"
     "0 1 0\n"
     "0 0 1\n"
     "1 0 1\n"
     "1 1 1\n"
     "0 1 1\n"
     "\n"
     "CELLS 8 32\n"
     "3 0 1 2\n"
     "3 0 7 4\n"
     "3 2 6 7\n"
     "3 0 2 7\n"
     "3 1 0 4\n"
     "3 1 4 7\n"
     "3 6 2 1\n"
     "3 1 7 6\n"
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

    logger(INFO,"To get VTK output, uncomment the relevant lines in the main function.");


    TriangulatedWallSelfTest t;
    //uncomment to get VTK output
    //t.setParticlesWriteVTK(true);
    //t.setWallsWriteVTK(FileType::ONE_FILE);
    t.solve();
    return 0;
}
