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
#include "Walls/TriangleWall.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Logger.h"

int main () {
    Mercury3D dpm;
    dpm.setTimeMax(1);
    dpm.setName("TriangleWallSelfTest");
    dpm.setDomain({0,0,0},{1,1,1});
    dpm.setGravity({0,0,0});

    LinearViscoelasticSpecies species;
    species.setDensity(constants::pi/6.0);
    Mdouble d = .1;
    species.setCollisionTimeAndRestitutionCoefficient(0.5,1.0,d*d*d);
    auto s = dpm.speciesHandler.copyAndAddObject(species);
    dpm.setTimeStep(0.01);

    SphericalParticle particle;
    particle.setSpecies(s);
    particle.setRadius(0.5*d);
    particle.setVelocity({0,0,-.1});
//    particle.setPosition({0,0,.5*d});
//    dpm.particleHandler.copyAndAddObject(particle);
    for (Mdouble x=-d/8; x<1+d; x+=d+d/22) {
        for (Mdouble y=-d/8; y<1+d-x; y+=d+d/22) {
            particle.setPosition(Vec3D(x, y, .5*d));
            dpm.particleHandler.copyAndAddObject(particle);
        }
    }

    TriangleWall wall;
    Vec3D A = {0,0,0};
    Vec3D B = {1,0,0};
    Vec3D C = {0,1,0};
    wall.setSpecies(s);
    wall.setVertices(A,B,C);
    dpm.wallHandler.copyAndAddObject(wall);
    dpm.setParticlesWriteVTK(true);
    dpm.setWallsWriteVTK(FileType::ONE_FILE);
    dpm.solve();
}
