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


    Mercury3D dpm;
    dpm.setName("TriangleWallCubeSelfTest");
    dpm.setMin({-.2,-.2,-.2});
    dpm.setMax({1.2,1.2,1.2});

    dpm.setGravity(10*Vec3D(0,0,-1));
    dpm.setTimeStep(1e-4);
    dpm.setTimeMax(1e4*dpm.getTimeStep());
    //setTimeMax(getTimeStep());
    dpm.setSaveCount(100);
    dpm.setXBallsAdditionalArguments("-v0 -solidf");
    dpm.fStatFile.setFileType(FileType::NO_FILE);
    dpm.restartFile.setSaveCount(1e5);

    Mdouble radius = 0.08;
    auto s = dpm.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    s->setDensity(1);
    s->setCollisionTimeAndRestitutionCoefficient(20.0*dpm.getTimeStep(),0.1,s->getMassFromRadius(radius));

    //create domain walls
    InfiniteWall i;
    i.setSpecies(dpm.speciesHandler.getLastObject());
    i.set({0,0,-1},dpm.getMin());
    dpm.wallHandler.copyAndAddObject(i);
    i.set({0,-1,0},dpm.getMin());
    dpm.wallHandler.copyAndAddObject(i);
    i.set({-1,0,0},dpm.getMin());
    dpm.wallHandler.copyAndAddObject(i);
    i.set({0,1,0},dpm.getMax());
    dpm.wallHandler.copyAndAddObject(i);
    i.set({1,0,0},dpm.getMax());
    dpm.wallHandler.copyAndAddObject(i);

    //read in walls
    dpm.wallHandler.readTriangleWall("TriangulatedWall.vtk",s);

    //introduce particles
    SphericalParticle p(s);

    //random insertion
    for (unsigned i=0; i<300; i++)
    {
        Vec3D r;
        r.X = dpm.random.getRandomNumber(dpm.getXMin(),dpm.getXMax());
        r.Y = dpm.random.getRandomNumber(dpm.getYMin(),dpm.getYMax());
        r.Z = dpm.random.getRandomNumber(dpm.getZMax(),dpm.getZMax()+1.0);
        p.setPosition(r);
        p.setRadius(dpm.random.getRandomNumber(radius,1.1*radius));
        if (dpm.checkParticleForInteraction(p))
            dpm.particleHandler.copyAndAddObject(p);
    }

    logger(INFO,"Inserted % particles",dpm.particleHandler.getNumberOfObjects());
    //uncomment to get VTK output
    dpm.setParticlesWriteVTK(true);
    dpm.setWallsWriteVTK(FileType::ONE_FILE);
    dpm.solve();
    return 0;
}
