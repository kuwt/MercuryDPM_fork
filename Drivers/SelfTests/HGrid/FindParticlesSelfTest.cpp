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

#include "Mercury2D.h"
#include "Walls/IntersectionOfWalls.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Boundaries/CubeInsertionBoundary.h"

using mathsFunc::cubic;
using mathsFunc::square;
using constants::pi;

/**
 * \brief This tests the hGridFindParticleContacts routine.
 * \details Particle are inserted into a square domain using an insertion boundary.
 * Then we find all particles close to the center and give them a velocity.
 */
int main() {
    // create dpm setup, set necessary variables (time step, domain, name, species)
    Mercury2D dpm;
    dpm.setTimeStep(1e-4);
    dpm.setDomain({0,0,0},{10,10,0});
    dpm.setName("FindParticles");

    // define and add species to dpm
    LinearViscoelasticSpecies species;
    species.setDensity(6.0 / pi);
    species.setCollisionTimeAndRestitutionCoefficient(1.0, 1.0, 1.0);
    dpm.speciesHandler.copyAndAddObject(species);

    // define an insertion boundary and use it to insert particles;
    // note, the boundary is not added to the dpm, it is only used to insert particles initially
    SphericalParticle particle;
    particle.setSpecies(dpm.speciesHandler.getLastObject());
    particle.setRadius(0.5);
    int maxFailed = 100;
    Vec3D posMin = dpm.getMin();
    Vec3D posMax = dpm.getMax();
    Vec3D velMin = {0, 0, 0};
    Vec3D velMax = {0, 0, 0};
    double radMin = 0.1;
    double radMax = 0.1;
    CubeInsertionBoundary insertionBoundary;
    insertionBoundary.set(&particle, maxFailed, posMin, posMax, velMin, velMax, radMin, radMax);
    //insertionBoundary.setVolumeFlowRate(10);
    insertionBoundary.checkBoundaryBeforeTimeStep(&dpm);
    logger(INFO,"number of particles % ",dpm.particleHandler.getNumberOfObjects());

    // use a search particle to find all particles close to the center of the domain, and give those particles a velocity
    SphericalParticle searchParticle;
    searchParticle.setSpecies(dpm.speciesHandler.getLastObject());
    searchParticle.setRadius(3);
    searchParticle.setPosition({5,5,0});
    auto particlesInContact = dpm.hGridFindParticleContacts(&searchParticle);
    for (auto p : particlesInContact) {
        p->setVelocity({1,0,0});
    }
    logger(INFO,"number of particles in contact with the search particle %",particlesInContact.size());

    //run the dpm to create the output files; it can be seen in the output that the particles in the center have velocity
    dpm.solve();
    logger(INFO,"All particles in contact to the search particle have a different velocity, see findCloseParticles.xballs");

    return 0;
}
