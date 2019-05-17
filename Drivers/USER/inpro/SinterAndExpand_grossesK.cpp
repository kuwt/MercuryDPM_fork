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
//#include "Species/LinearViscoelasticSpecies.h"
#include "Species/SinterFrictionSpecies.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "/usr/lib/openmpi/include/mpi.h"
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
	dpm.setDomain({0,0,0},{0.030,0.040,0.010});
    dpm.setName("SinterAndExpand_grossesK");
    dpm.setTimeMax(0.2);
    dpm.setGravity({0,-9.80665,0});
    dpm.setSaveCount(1000);
	//dpm.readRestartFile("Filling.restart");
	std::cout << "Hello" << std::endl; 
	SinterFrictionSpecies species;
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
	double stiffness=species.getLoadingStiffness();
	species.setPlasticParameters(stiffness, stiffness, 0, 0.16);
	species.setSinterAdhesion(0.0013*stiffness);
    species.setSlidingStiffness(2.0/7.0*stiffness);
	species.setSlidingDissipation(2.0/7.0*species.getDissipation());
	species.setSlidingFrictionCoefficient(0.5);
    species.setRollingStiffness(2.0/5.0*stiffness);
	species.setRollingDissipation(2.0/5.0*species.getDissipation());
	species.setRollingFrictionCoefficient(0.1);
	
	
	
	dpm.setTimeStep(0.05*tc);
    logger(INFO,"Timestep %",dpm.getTimeStep());
	//ParticleSpecies* species = dpm.speciesHandler.getLastObject();
	species.setTemperatureDependentDensity(
		[] (double temperature) {return 100*296/temperature;}
	);
	species.setSinterRate(0.01);
    dpm.speciesHandler.copyAndAddObject(species);
	
	dpm.setNumberOfDomains({2,1,1});
	
	ThermalParticle particle;
	
	particle.setSpecies(dpm.speciesHandler.getLastObject());
	dpm.particleHandler.copyAndAddObject(particle);
	dpm.readDataFile("Filling_Stl_Test.data",14); //neue species erzeugen, wenn datafile eingelesen wird
	std::cout << "Hello" << std::endl; 
	//dpm.setRestarted(false);
	
	//dpm.setSaveCount(1);

	dpm.write(std::cout,false);
	dpm.boundaryHandler.clear();

	
	for (BaseParticle* particle : dpm.particleHandler) {
		ThermalParticle* thermalParticle = dynamic_cast<ThermalParticle*>(particle);
		logger.assert_always(thermalParticle,"Particles need to be of type ThermalParticle");
		thermalParticle->setTimeDependentTemperature(
			[] (double time) {return 296+200*time;}	
		);
	}
	
	
	dpm.setParticlesWriteVTK(true);
	dpm.wallHandler.readTriangleWall("grossesK_binaer.stl",&species,0.001);
    dpm.solve();
	
	//logger(INFO,"Hello");
    return 0;
}
