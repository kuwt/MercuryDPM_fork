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

#include <Boundaries/CubeInsertionBoundary.h>
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
    Mercury3D dpm;
    dpm.setName("STLReader");
    dpm.removeOldFiles();
    dpm.setMax({250e-3,40e-3,10e-3});
    dpm.setMin(-dpm.getMax());

    dpm.setGravity(9.8*Vec3D(0,0,-1));
    dpm.setTimeStep(5e-4);
    dpm.setTimeMax(1);
    //setTimeMax(getTimeStep());
    dpm.setSaveCount(10);
    dpm.setXBallsAdditionalArguments("-v0 -solidf");
    dpm.fStatFile.setFileType(FileType::NO_FILE);
    dpm.restartFile.setSaveCount(1e5);

    Mdouble radius = 5e-3;
    auto species = dpm.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    species->setDensity(200);
    species->setCollisionTimeAndRestitutionCoefficient(20.0*dpm.getTimeStep(),0.1,species->getMassFromRadius(radius));

     //read in walls
    dpm.wallHandler.readTriangleWall("grossesK_binaer.stl",species,0.001);

    //introduce particles
    SphericalParticle p(species);
    p.setRadius(radius);
    p.setPosition({0,0,5e-3});
    //dpm.particleHandler.copyAndAddObject(p);

    dpm.particleHandler.copyAndAddObject(p);
    CubeInsertionBoundary insertionBoundary;
    Vec3D pMin = {-.22,-.02+radius,radius};
    Vec3D pMax = {+.22,+.02+radius,.016-radius};
    insertionBoundary.set(&p, 100, pMin, pMax, {0,0,0}, {0,0,0}, radius, radius);
    //insertionBoundary.set(&p, 100, {0,5e-3,0}, {0,5e-3,0}, {0,0,0}, {0,0,0}, radius, radius);
    insertionBoundary.checkBoundaryBeforeTimeStep(&dpm);
    logger(INFO,"Inserted % particles",dpm.particleHandler.getNumberOfObjects());

    //uncomment to get VTK output
    dpm.setParticlesWriteVTK(true);
    //dpm.setWallsWriteVTK(FileType::ONE_FILE);
    dpm.solve();
    return 0;
}
