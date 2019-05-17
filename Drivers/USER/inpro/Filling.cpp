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
#include "Particles/ThermalParticle.h"
#include "Walls/InfiniteWall.h"
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
    Mercury3D dpm;
    dpm.setDomain({0,0,0},{2e-2,2e-2,2e-2});
    dpm.setName("Filling");
    dpm.setTimeMax(0.2);
    dpm.setGravity({0,0,-9.8});
    dpm.setSaveCount(500);

    // define and add species to dpm
    LinearViscoelasticSpecies species;
	species.setDensity(100);
    species.setHandler(&dpm.speciesHandler);
    //species.setCollisionTimeAndRestitutionCoefficient(5e-3, 0.5, 0.01);
    double r = 1.2e-3;
    double E = 250e6;
    double vmax = 10;
    double m = species.getMassFromRadius(r);
    //double k = E*r*0.2;
    double k = m*square(vmax)/square(0.1*r);
    species.setStiffnessAndRestitutionCoefficient(k,0.5,m);
    double tc = species.getCollisionTime(m);
    dpm.setTimeStep(0.1*tc);
    logger(INFO,"Timestep %",dpm.getTimeStep());
    dpm.speciesHandler.copyAndAddObject(species);

    // define an insertion boundary and use it to insert particles;
    // note, the boundary is not added to the dpm, it is only used to insert particles initially
    ThermalParticle particle;
	particle.setTemperature(23);
    particle.setSpecies(dpm.speciesHandler.getLastObject());
	//particle.setTimeDependentTemperature(
	//	[] (double time) {return 296+20*time;}	
	//);
    int maxFailed = 100;
    Vec3D posMin = dpm.getMin();
    Vec3D posMax = dpm.getMax();
    Vec3D velMin = {0, 0, 0};
    Vec3D velMax = {0, 0, 0};
    double radMin = r;
    double radMax = r;
    //dpm.particleHandler.copyAndAddObject(particle);
	CubeInsertionBoundary insertionBoundary;
    insertionBoundary.set(&particle, maxFailed, posMin, posMax, velMin, velMax, radMin, radMax);
    //insertionBoundary.setVolumeFlowRate(10);
    //insertionBoundary.checkBoundaryBeforeTimeStep(&dpm);
	dpm.boundaryHandler.copyAndAddObject(insertionBoundary);
    logger(INFO,"number of particles % ",dpm.particleHandler.getNumberOfObjects());

    InfiniteWall w;
    w.setSpecies(dpm.speciesHandler.getLastObject());
    w.set({0,0,-1},dpm.getMin());
    dpm.wallHandler.copyAndAddObject(w);
    w.set({0,0,1},dpm.getMax());
    dpm.wallHandler.copyAndAddObject(w);
    w.set({0,-1,0},dpm.getMin());
    dpm.wallHandler.copyAndAddObject(w);
    w.set({0,1,0},dpm.getMax());
    dpm.wallHandler.copyAndAddObject(w);
    w.set({-1,0,0},dpm.getMin());
    dpm.wallHandler.copyAndAddObject(w);
    w.set({1,0,0},dpm.getMax());
    dpm.wallHandler.copyAndAddObject(w);

    dpm.solve();

    return 0;
}
