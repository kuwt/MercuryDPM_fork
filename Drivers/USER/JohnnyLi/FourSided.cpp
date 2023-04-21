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

int main()
{
    //create an instance of Mercury3D
    Mercury3D dpm;
    //set the name of the output files
    dpm.setName("FourSided");
    //set a simple spring-dashpot contact law
    LinearViscoelasticSpecies* s = dpm.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    Mdouble radius = 0.04;
    Mdouble collisionTime = 5e-2;
    s->setCollisionTimeAndRestitutionCoefficient(collisionTime,0.5,s->getMassFromRadius(radius));
    logger(INFO,"k % gam %",s->getStiffness(),s->getDissipation());
    // create a level set wall of type fourSided
    LevelSetWall* w = dpm.wallHandler.copyAndAddObject(LevelSetWall(LevelSetWall::Shape::FourSided, 0.8, s));
    // give the wall angular velocity
    w->setAngularVelocity(Vec3D(0,0,1)*-0.05);
    // add walls on the domain boundary
    dpm.wallHandler.copyAndAddObject(InfiniteWall({-1,0,0},{-1,0,0},s));
    dpm.wallHandler.copyAndAddObject(InfiniteWall({ 1,0,0},{ 0,0,0},s));
    dpm.wallHandler.copyAndAddObject(InfiniteWall({0,-1,0},{0,-1,0},s));
    dpm.wallHandler.copyAndAddObject(InfiniteWall({0, 1,0},{0, 0,0},s));
    dpm.wallHandler.copyAndAddObject(InfiniteWall({0,0,-1},{0,0,-radius},s));
    dpm.wallHandler.copyAndAddObject(InfiniteWall({0,0, 1},{0,0, radius},s));
    // make the walls at teh domain boundary invisible
    for (unsigned i=1; i<dpm.wallHandler.getNumberOfObjects(); ++i) {
        dpm.wallHandler.getObject(i)->setVTKVisibility(false);
    }
    // Set time step, max time, domain size, savecount
    dpm.setTimeMax(20);
    dpm.setTimeStep(collisionTime/15);
    dpm.setDomain({-1,-1,-1},{1,1,1});
    dpm.setSaveCount(dpm.getTimeMax()/dpm.getTimeStep()/100);
    // add particles in a limited region
    SphericalParticle p(s);
    Mdouble d;
    Vec3D n;
    p.setRadius(radius);
    for (Mdouble x=-1+radius; x<=0; x+=2*radius)
        for (Mdouble y=-1+radius; y<=0; y+=2*radius)
            for (Mdouble z=0; z<=0; z+=2*radius) {
                p.setPosition({x,y,z});
                if (!w->getDistanceAndNormal(p,d,n)) {
                    dpm.particleHandler.copyAndAddObject(p);
                }
            }
    logger(INFO,"#particles %",dpm.particleHandler.getNumberOfObjects());
    // turn on vtk output
    dpm.wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);
    dpm.setParticlesWriteVTK(true);
    // start the simulation
    dpm.solve();
    return 0;
}
