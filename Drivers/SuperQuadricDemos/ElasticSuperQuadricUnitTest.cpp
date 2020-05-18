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

#include <Species/HertzianViscoelasticMindlinSpecies.h>
#include <Walls/InfiniteWall.h>
#include "Mercury3D.h"
#include "Particles/SuperQuadricParticle.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Math/ExtendedMath.h"

/**
 * Takes an oblique collision of elastic superellipsoids and checks whether the collision
 *  - conserves momentum
 *  - conserves kinetic+potential energy
 */

int main(int argc, char* argv[])
{
    Mercury3D problem;
    problem.setName("ElasticSuperQuadricUnitTest");
    problem.setDomain({-1,-1,-1},{1,1,1});
    problem.setSaveCount(NEVER);
    problem.dataFile.setSaveCount(10);
    
    auto s = problem.speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
    s->setStiffness(2e5);
    s->setDensity(6.0/constants::pi);
    
    SuperQuadricParticle p;
    p.setSpecies(s);
    p.setAxes(.5,.5,.5);
    p.setExponents(.5,.9);
    p.setPosition({.1,0,.5});
    p.setVelocity({0,0,-.5});
    auto p0 = problem.particleHandler.copyAndAddObject(p);
    
    p.setPosition(-p.getPosition());
    p.setVelocity(-p.getVelocity());
    //p.rotate({0,constants::pi,0});
    auto p1 = problem.particleHandler.copyAndAddObject(p);
    
    problem.setTimeMax(1e-0);
    problem.setTimeStep(1e-4);
    
    const Vec3D momentum0 = problem.particleHandler.getMomentum();
    const Vec3D angularMomentum0 = problem.particleHandler.getAngularMomentum();
    const Mdouble kineticEnergy0 = problem.particleHandler.getKineticEnergy()
            +problem.particleHandler.getRotationalEnergy();

    // comment the  next to lines to turn on file output
    problem.setSuperquadricParticlesWriteVTK(true);
    problem.setFileType(FileType::NO_FILE);
    problem.solve();

    const Vec3D momentum1 = problem.particleHandler.getMomentum();
    const Vec3D angularMomentum1 = problem.particleHandler.getAngularMomentum();
    const Mdouble kineticEnergy1 = problem.particleHandler.getKineticEnergy()
                                   +problem.particleHandler.getRotationalEnergy();

    // Check conservation properties
    helpers::check(momentum0-momentum1,{0,0,0},1e-10*momentum0.getLength(),
            "Conservation of momentum");
    helpers::check(angularMomentum0-angularMomentum1,{0,0,0},1e-10*angularMomentum0.getLength(),
    "Con. of angular momentum");
    helpers::check(kineticEnergy0-kineticEnergy1,0,1e-4*kineticEnergy0,
            "Con. of kinetic energy  ");
    // Note, energy conservation is not very good; it does not decay quickly with time step

    return 0;
}
