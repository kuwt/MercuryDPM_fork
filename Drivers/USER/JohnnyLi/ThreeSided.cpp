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
#include "Walls/LevelSetWall.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticSpecies.h"

int main()
{
    //create an instance of Mercury3D
    Mercury3D dpm;
    //set the name of the output files
    dpm.setName("ThreeSided");
    //set a simple spring-dashpot contact law
    LinearViscoelasticSpecies* s = dpm.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    Mdouble radius = 10;
    Mdouble collisionTime = 5e-2;
    s->setCollisionTimeAndRestitutionCoefficient(collisionTime,0.5,s->getMassFromRadius(radius));
    logger(INFO,"k % gam %",s->getStiffness(),s->getDissipation());
    // Set time step, max time, domain size, savecount
    Vec3D centerOfRotation = {2240,1480,50};
    Vec3D halfWidth = {120,110,50};
    dpm.setTimeMax(20);
    dpm.setTimeStep(collisionTime/15);
    dpm.setDomain(centerOfRotation-halfWidth,centerOfRotation+halfWidth);
    dpm.setSaveCount(dpm.getTimeMax()/dpm.getTimeStep()/100);
    dpm.setGravity({0,-1,0});
    // create a level set wall of type fourSided
    // give the wall angular velocity
    dpm.wallHandler.readTriangleWall("3S.stl",s,1,centerOfRotation,{0,0,0},{0,0,-0.05});
    // add walls on the domain boundary
    dpm.wallHandler.copyAndAddObject(InfiniteWall({-1,0,0},dpm.getMin(),s));
    dpm.wallHandler.copyAndAddObject(InfiniteWall({ 1,0,0},dpm.getMax(),s));
    dpm.wallHandler.copyAndAddObject(InfiniteWall({0,-1,0},dpm.getMin(),s));
    //dpm.wallHandler.copyAndAddObject(InfiniteWall({0, 1,0},dpm.getMax(),s));
    dpm.wallHandler.copyAndAddObject(InfiniteWall({0,0,-1},dpm.getMin(),s));
    dpm.wallHandler.copyAndAddObject(InfiniteWall({0,0, 1},dpm.getMax(),s));
    // make the walls at teh domain boundary invisible
    //for (unsigned i=dpm.wallHandler.getNumberOfObjects()-2; i<dpm.wallHandler.getNumberOfObjects(); ++i) {
    //    dpm.wallHandler.getObject(i)->setVTKVisibility(false);
    //}
    // add particles in a limited region
    CubeInsertionBoundary b;
    SphericalParticle p(s);
    b.set(p,1000,dpm.getMin()+Vec3D(0,2*halfWidth.Y,0),dpm.getMax()+Vec3D(0,.5*halfWidth.Y,0),{0,0,0},{0,0,0},radius,radius);
    b.insertParticles(&dpm);
    logger(INFO,"#particles %",dpm.particleHandler.getNumberOfObjects());
    // turn on vtk output
    dpm.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    dpm.setParticlesWriteVTK(true);
    // start the simulation
    dpm.solve();
    return 0;
}
