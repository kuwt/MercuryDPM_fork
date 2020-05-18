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

    //output of the unit test: check if evaluation works
    Vec3D val = nurbsSurface.evaluate(0.5,0.5);
    helpers::check(val,{.6,.8,1},std::numeric_limits<double>::epsilon(),"Nurbs evaluation failed");

    //now demonstrate how nurbs wall can be defined:

    //create a basic dpm class
    DPMBase dpm;
    dpm.setName("NurbsSurfaceUnitTest");
    auto s = dpm.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    SphericalParticle p;
    p.setSpecies(s);
    p.setRadius(1);
    p.setPosition({1,1,.2});
    dpm.particleHandler.copyAndAddObject(p);

    //define nurbs wall and add to dpm class
    NurbsWall nurbsWall;
    nurbsWall.setSpecies(s);
    nurbsWall.set(nurbsSurface);
    auto w = dpm.wallHandler.copyAndAddObject(nurbsWall);

    //write vtk file
    dpm.setWallsWriteVTK(FileType::ONE_FILE);
    dpm.forceWriteOutputFiles();

    //test if contact can be found
    Vec3D normal;
    double distance;
    w->getDistanceAndNormal(p,distance,normal);
    helpers::check(normal,{-sqrt(.5),-sqrt(.5),0},4*std::numeric_limits<double>::epsilon(),"Nurbs contact evaluation");

    return 0;
}
