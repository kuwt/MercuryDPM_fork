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

#include <Species/LinearViscoelasticSpecies.h>
#include <Walls/Screw.h>
#include "Mercury3D.h"
/**
 * Testing if screw geometry is implemented well
 */
int main()
{
    Mercury3D dpm;
    
    LinearViscoelasticSpecies species;
    species.setDensity(1);
    species.setStiffness(1);
    auto s = dpm.speciesHandler.copyAndAddObject(species);
    
    Screw screw({0,0,0}, 1, 1, 1, 0, .1, ScrewType::singleHelix);
    screw.setSpecies(s);
    auto w = dpm.wallHandler.copyAndAddObject(screw);
    
    dpm.setDomain({-0.2,-1,-1},{1.2,1,1});
    dpm.setTimeStep(1e-12);
    dpm.setTimeMax(dpm.getTimeStep());
    dpm.setName("ScrewUnitTest");
    dpm.setParticlesWriteVTK(true);
    dpm.setWallsWriteVTK(true);
    
    Mdouble h = 0.05;
    Mdouble distance;
    Vec3D normal;
    SphericalParticle particle(s);
    particle.setRadius(0.5*h);
    for (Mdouble x = dpm.getXMin(); x<dpm.getXMax(); x+=h) {
        for (Mdouble y = dpm.getYMin(); y<dpm.getYMax(); y+=h) {
            for (Mdouble z = dpm.getZMin(); z<dpm.getZMax(); z+=h)
            {
                particle.setPosition({x, y, z});
                if (w->getDistanceAndNormal(particle,distance,normal))
                {
                    dpm.particleHandler.copyAndAddObject(particle);
                    //logger(INFO,"p %",particle.getPosition());
                }
            }
        }
    }
    logger(INFO,"Inserted % particles",dpm.particleHandler.getSize());

    //dpm.solve();
    helpers::check(dpm.particleHandler.getSize(),8213,0,"Screw surface was wrongly detected");

    return 0;
}
