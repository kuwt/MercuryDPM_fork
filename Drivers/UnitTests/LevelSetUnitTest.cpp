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
#include "Walls/LevelSetWall.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticSpecies.h"

int main() {

    DPMBase dpm;
    dpm.setName("LevelSetUnitTest");
    dpm.removeOldFiles();
    LinearViscoelasticSpecies* s = dpm.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    Mdouble r = 0.04;
    s->setCollisionTimeAndRestitutionCoefficient(0.5,0.5,s->getMassFromRadius(r));
    logger(INFO,"k % gam %",s->getStiffness(),s->getDissipation());
    LevelSetWall* w = dpm.wallHandler.copyAndAddObject(LevelSetWall(LevelSetWall::Shape::Diamond, 0.8, s));
    //w->setPosition({0,0,0}); //change this to check the transation properties
    w->setVelocity(Vec3D(-1,-1,0)*0.01);
    w->setAngularVelocity(Vec3D(0,0,1)*-0.05);
    dpm.wallHandler.copyAndAddObject(InfiniteWall({-1,0,0},{-1,0,0},s));
    dpm.wallHandler.copyAndAddObject(InfiniteWall({ 1,0,0},{ 0,0,0},s));
    dpm.wallHandler.copyAndAddObject(InfiniteWall({0,-1,0},{0,-1,0},s));
    dpm.wallHandler.copyAndAddObject(InfiniteWall({0, 1,0},{0, 0,0},s));
    dpm.wallHandler.copyAndAddObject(InfiniteWall({0,0,-1},{0,0,-r},s));
    dpm.wallHandler.copyAndAddObject(InfiniteWall({0,0, 1},{0,0, r},s));
    for (unsigned i=1; i<dpm.wallHandler.getNumberOfObjects(); ++i) {
        dpm.wallHandler.getObject(i)->setVTKVisibility(false);
    }
    dpm.setTimeMax(20);
    dpm.setTimeStep(0.01);
    dpm.setDomain({-1,-1,-1},{1,1,1});
    dpm.setGravity(Vec3D(1,1,0)*0.5*s->getStiffness());
    SphericalParticle p(s);
    Mdouble d;
    Vec3D n;
    p.setRadius(r);
    for (Mdouble x=-1+r; x<=0; x+=2*r)
        for (Mdouble y=-1+r; y<=0; y+=2*r)
            for (Mdouble z=0; z<=0; z+=2*r) {
                p.setPosition({x,y,z});
                if (!w->getDistanceAndNormal(p,d,n)) {
                    dpm.particleHandler.copyAndAddObject(p);
                }
            }
    logger(INFO,"#particles %",dpm.particleHandler.getNumberOfObjects());
    dpm.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    dpm.setParticlesWriteVTK(true);
    dpm.solve();
    w->writeToFile(10,0.0);
    return 0;
}
