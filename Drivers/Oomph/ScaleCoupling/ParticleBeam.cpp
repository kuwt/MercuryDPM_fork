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

#include "Mercury3D.h"
#include "Species/LinearViscoelasticReversibleAdhesiveSpecies.h"

/// A 1D "beam" structure composed of particles
class ParticleBeam : public Mercury3D {
    double length = 20.0; // length of beam
    double distance = 0.4; // distance between particles
    double bulkDensity = 1309; // bulk density (assuming a cubic packing)
    double elasticModulus = 1e8; // =stiffness/distance
    std::function<Vec3D(double)> velocityAtBoundary =
            [] (double) {return Vec3D(0.1,0,0);}; // velocity of particle at left boundary

public:
    ParticleBeam()
    {
        double overlap = 0.05*distance; //overlap between particles
        double radius = 0.5*(distance+overlap); // particle radius
        double particleDensity = bulkDensity*mathsFunc::cubic(distance) / (constants::pi/6.0*mathsFunc::cubic(2.0*radius));
        double stiffness = elasticModulus*distance; // particle stiffness

        setName("ParticleBeam");
        setDomain(Vec3D(-radius,-radius,0),Vec3D(radius,radius,length));
        setParticlesWriteVTK(true);

        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticReversibleAdhesiveSpecies());
        species->setDensity(particleDensity);
        species->setStiffness(stiffness);
        species->setAdhesionStiffness(stiffness);
        species->setAdhesionForceMax(stiffness*overlap);

        SphericalParticle p(species);
        p.setRadius(radius);
        auto n = (unsigned)(length/distance);
        for (int i = 0; i<n; ++i) {
            for (int j = 0; j<2; ++j) {
                for (int k = 0; k<2; ++k) {
                    p.setPosition(Vec3D(1+2*i, j==0?-1:1, k==0?-5:-3)*0.5*distance);
                    particleHandler.copyAndAddObject(p);
                    if (i==0) {
                        particleHandler.getLastObject()->fixParticle();
                        particleHandler.getLastObject()->setPrescribedVelocity(velocityAtBoundary);
                    }
                }
            }
        }

        setTimeStep(species->computeTimeStep(p.getMass()));

        double waveSpeed = sqrt(elasticModulus/bulkDensity);
        double propagationTime = length/waveSpeed;
        setTimeMax(1.25*propagationTime);
        setSaveCount(getTimeMax()/getTimeStep()/100.0);

        logger(INFO,"propagationTime % timeMax %", propagationTime, getTimeMax());
    }

    void printTime() const override
    {
        logger(INFO, "t=%3.6, tMax=%3.6, mom=%3.6",
               getTime(), getTimeMax(), particleHandler.getMomentum());
    }
};

int main() {
    ParticleBeam dpm;
    dpm.write(std::cout);
    dpm.solve();
    helpers::writeToFile(dpm.getName()+".gnu","p 'ParticleBeam.data' u 1 ev 101::99");
    return 0;
}
