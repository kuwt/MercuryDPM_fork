//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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
#include "Walls/LevelSetWall.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticSpecies.h"

int main()
{
    //create an instance of Mercury3D
    Mercury3D dpm;
    //create species
    LinearViscoelasticSpecies* s = dpm.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    Mdouble radius = 0.04;
    Mdouble collisionTime = 5e-2;
    s->setCollisionTimeAndRestitutionCoefficient(collisionTime,0.5,s->getMassFromRadius(radius));
    //add stl walls
    unsigned g0 = dpm.wallHandler.readTriangleWall("wgcp_ac_inlet_3.dem_adsorber.stl.stl"  ,s);
    unsigned g1 = dpm.wallHandler.readTriangleWall("wgcp_ac_inlet_3.dem_boxsection.stl.stl",s);
    unsigned g2 = dpm.wallHandler.readTriangleWall("wgcp_ac_inlet_3.dem_floor.stl.stl"     ,s);
    unsigned g3 = dpm.wallHandler.readTriangleWall("wgcp_ac_inlet_3.dem_I-type-1.stl.stl"  ,s);
    unsigned g4 = dpm.wallHandler.readTriangleWall("wgcp_ac_inlet_3.dem_I-type-2.stl.stl"  ,s);
    unsigned g5 = dpm.wallHandler.readTriangleWall("wgcp_ac_inlet_3.dem_adsorber.stl.stl"  ,s);
    unsigned g6 = dpm.wallHandler.readTriangleWall("wgcp_ac_inlet_3.dem_rings.stl.stl"     ,s);
//    //shift them
//    for (auto& w: dpm.wallHandler) {
//        unsigned g = w->getGroupId();
//        if (g==g1) w->move({10,0,0});
//        else if (g==g2) w->move({20,0,0});
//        else if (g==g3) w->move({30,0,0});
//        else if (g==g4) w->move({40,0,0});
//        else if (g==g5) w->move({50,0,0});
//        else if (g==g6) w->move({60,0,0});
//    }
    // define values required for each simulation
    dpm.wallHandler.setWriteVTK(true);
    dpm.setTimeMax(1);
    dpm.setTimeStep(1);
    dpm.setName("WGCP");
    dpm.setMax(1,1,1);
    dpm.solve();
    return 0;
}
