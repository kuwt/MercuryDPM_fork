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

#include "Nurbs/NurbsSurface.h"
#include "Walls/NurbsWall.h"
#include "Math/Helpers.h"
#include <sstream>
#include <DPMBase.h>
#include <Interactions/NormalForceInteractions/LinearViscoelasticInteraction.h>
#include <Species/LinearViscoelasticSpecies.h>

int main()
{
    //define quarter circle as nurbs surface
    std::vector<double> knotsU = {0,0,0,1,1,1};
    std::vector<double> knotsV = {0,0,1,1};
    std::vector<std::vector<Vec3D>> controlPoints = {{{1,0,0},{1,0,2}},{{1,1,0},{1,1,2}},{{0,1,0},{0,1,2}}} ;
    std::vector<std::vector<Mdouble>> weights = {{1,1},{1,1},{2,2}};
    NurbsSurface nurbsSurface(knotsU,knotsV,controlPoints,weights);

    //create a basic dpm class
    DPMBase dpm;
    dpm.setName("NurbsSurfaceUnitTest");
    auto s = dpm.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());

    //define nurbs wall and add to dpm class
    NurbsWall nurbsWall;
    nurbsWall.setSpecies(s);
    nurbsWall.set(nurbsSurface);
    auto w = dpm.wallHandler.copyAndAddObject(nurbsWall);

    SphericalParticle particle;
    particle.setSpecies(s);
    particle.setRadius(.99);
    for (double x=-.5; x<=4.5; x+=.125) {
        for (double z = -.5; z <= 2.5; z += .5) {
            particle.setPosition({x, 2-x, z});
            dpm.particleHandler.copyAndAddObject(particle);
        }
    }

    //test if contact can be found
    Vec3D normal;
    double distance;
    std::stringstream ss;
    for (auto p : dpm.particleHandler) {
        if (w->getDistanceAndNormal(*p, distance, normal)) {
            p->setPosition(distance * normal + p->getPosition());
            ss << p->getPosition() << '\t' << distance << '\t' << distance * normal + p->getPosition() << '\n';
        } else {
            p->setPosition({0,0,0});
        }
    }
    helpers::writeToFile("NurbsWallUnitTest.txt",ss.str());

    //write vtk file
    dpm.setName("NurbsWallUnitTest");
    dpm.setWallsWriteVTK(FileType::ONE_FILE);
    dpm.setParticlesWriteVTK(true);
    dpm.forceWriteOutputFiles();

    return 0;
}
