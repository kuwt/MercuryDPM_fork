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
#include "Walls/InfiniteWall.h"
#include "Math/Helpers.h"
#include <sstream>
#include <DPMBase.h>
#include <Interactions/NormalForceInteractions/LinearViscoelasticInteraction.h>
#include <Species/LinearViscoelasticSpecies.h>

/**
 * Tests if contact with nurbs surface is correctly implemented.
 * The setup consists of five particles and two walls, on of which is a NurbsWall.
 * The particles fall downwards by gravity, filling the space between the two walls.
 */
int main()
{
    //define quarter circle as nurbs surface
    std::vector<double> knotsU = {0,0,0,1,1,1};
    std::vector<double> knotsV = {0,0,1,1};
    std::vector<std::vector<Vec3D>> controlPoints = {{{1,0,0},{1,0,2}},{{1,1,0},{1,1,2}},{{0,1,0},{0,1,2}}};
    std::vector<std::vector<Mdouble>> weights = {{1,1},{1,1},{2,2}};
    NurbsSurface nurbsSurface(knotsU,knotsV,controlPoints,weights);
    //nurbsSurface.flipOrientation();

    //create a basic dpm class
    DPMBase dpm;
    dpm.setName("NurbsWallSelfTest");
    dpm.setGravity({0,0,-1});
    auto s = dpm.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    s->setCollisionTimeAndRestitutionCoefficient(.05,.1,s->getMassFromRadius(0.1));

    //define nurbs wall and add to dpm class
    NurbsWall nurbsWall;
    nurbsWall.setSpecies(s);
    nurbsWall.set(nurbsSurface);
    nurbsWall.setOrientationViaNormal({0,0,1});
    nurbsWall.setPosition({1,0,0});
    dpm.wallHandler.copyAndAddObject(nurbsWall);

    InfiniteWall infiniteWall;
    infiniteWall.setSpecies(s);
    infiniteWall.set({0,1,0},{1,1,1});
    dpm.wallHandler.copyAndAddObject(infiniteWall);

    SphericalParticle particle;
    particle.setSpecies(s);
    particle.setRadius(.1);
    for (double x=0; x<=0; x+=.2) {
        for (double y = .1; y <= .9; y += .2) {
            for (double z = 1.1; z <= 1.1; z += .2) {
                particle.setPosition({x, y, z});
                dpm.particleHandler.copyAndAddObject(particle);
            }
        }
    }

    dpm.setWallsWriteVTK(FileType::ONE_FILE);
    dpm.setParticlesWriteVTK(true);
    dpm.setTimeStep(0.001);
    dpm.setTimeMax(2.0);
    dpm.setSaveCount(100);
    dpm.setDomain({-1,0,0},{1,1,1});
    dpm.solve();

    return 0;
}
