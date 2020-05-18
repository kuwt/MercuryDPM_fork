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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include "Funnel.h"

int main(int argc UNUSED, char *argv[] UNUSED)
{
	Funnel problem;
	
	// Problem parameters
	problem.setName("dec13_A254_Hi0075_RC06_MU05");

    auto species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
    species->setSlidingFrictionCoefficient(0.5);
	species->setDensity(1442.0);

	//problem.set_HGRID_max_levels(2);
	//problem.set_HGRID_num_buckets(1e6);
	// Particle properties

	problem.setFixedParticleRadius(300e-6);
	problem.setInflowParticleRadius(300e-6);

	species->setCollisionTimeAndRestitutionCoefficient(4e-4, 0.6, species->getMassFromRadius(problem.getFixedParticleRadius()));


	//problem.setStiffness(1e5*4/3*pi*problem.getInflowParticleRadius()*9.81*problem.getDensity());
	//problem.set_disp(50*sqrt(9.81/(2*problem.getInflowParticleRadius());
	std::cout << "Setting k to " << species->getStiffness() << " and disp to " << species->getDissipation() << std::endl;

 	species->setSlidingStiffness(species->getStiffness()*2.0/7.0);
 	species->setSlidingDissipation(species->getDissipation()*2.0/7.0);

    double mass = species->getMassFromRadius(0.5 * (problem.getMinInflowParticleRadius() + problem.getMaxInflowParticleRadius()));
    problem.setTimeStep(0.02 * species->getCollisionTime(mass));
	problem.setTimeMax(3);
	problem.setSaveCount(100);
	
	problem.setChuteLength(0.6);
	problem.setChuteWidth(0.25);
	problem.setChuteAngle(25.4);
	problem.setRoughBottomType(MONOLAYER_DISORDERED);
	


	double funx = problem.getXMin()+0.5*(problem.getXMax()-problem.getXMin());
	

	problem.set_funa(60.0);
	problem.set_funD(0.015);
	problem.set_funHf(0.075+(problem.getXMax()-funx)*sin(problem.getChuteAngle()));
	problem.set_funnz(50.0);
	problem.set_funO(-funx, 0.5*(problem.getYMax()-problem.getYMin()));
	problem.set_funfr(0.3);
	
	problem.setInflowVelocity(0);
	problem.setInflowVelocityVariance(0.01);
	problem.setMaxFailed(1);
	

	
	//solve
	//cout << "Maximum allowed speed of particles: " << problem.getMaximumVelocity() << endl; // speed allowed before particles move through each other!
	

	std::cout << "dt=" << problem.getTimeStep() << std::endl;
	problem.solve();
	problem.writeRestartFile();
}
