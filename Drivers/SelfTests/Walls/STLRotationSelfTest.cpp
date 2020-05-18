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

#include <Mercury3D.h>
#include <Species/LinearViscoelasticSpecies.h>
#include <Walls/TriangleWall.h>

int main(int argc, char* argv[])
{
    // Problem setup
    Mercury3D problem;
    // Set name of output files
    problem.setName("STLRotationSelfTest");
    // Set domain size
    problem.setMax({450,50,50});
    problem.setMin({-450,-50,-50});
    // Set time step, final time and how often to output
    problem.setTimeStep(0.01);
    problem.setTimeMax(0.1);
    problem.setSaveCount(1);
    // Turn off most output files except data and wall-vtu files
    problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.eneFile.setFileType(FileType::NO_FILE);
    problem.restartFile.setFileType(FileType::NO_FILE);
    problem.setWallsWriteVTK(true);
    // Introduce an material (no properties set, as no collisions happen here)
    problem.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    // introduce an outer wall
    problem.wallHandler.readTriangleWall("Casing.stl",problem.speciesHandler.getObject(0));
    // introduce an inner, rotating wall
    Mdouble scaleFactor = 1;
    Vec3D shift = {0,0,0};
    Vec3D velocity = {0,0,0};
    Vec3D angularVelocity = {2.0*constants::pi,0,0};
    problem.wallHandler.readTriangleWall("Screw.stl",problem.speciesHandler.getObject(0),
            scaleFactor,shift,velocity,angularVelocity);
    // start solving in time
    problem.solve();
    logger.assert_always(problem.wallHandler.getSize()==4316,"Didn't read the right number of walls");
    logger(INFO,"Load %Wall_*.vtu in paraview to see the wall geometry",problem.getName());
    return 0;
}
