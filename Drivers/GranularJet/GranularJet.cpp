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

	// Particle properties
	problem.setFixedParticleRadius (.3e-3);
	problem.setInflowParticleRadius(.3e-3);

	//Contact properties
    auto species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
 	species->setSlidingFrictionCoefficient(0.5);
	species->setDensity(1442.0); //sand
	species->setCollisionTimeAndRestitutionCoefficient(4e-4, 1.0-(0.12*(50.0/15.0)), species->getMassFromRadius(problem.getInflowParticleRadius()));//eps=0.6
 	species->setSlidingStiffness(species->getStiffness()*2.0/7.0);
 	species->setSlidingDissipation(species->getDissipation()*2.0/7.0);

	std::cout << "Setting k to " << species->getStiffness() << " and disp to " << species->getDissipation() << " with radius: " <<  problem.getInflowParticleRadius()<< std::endl;
	
	// Default problem parameters
	problem.setChuteLength(0.25);
	problem.setChuteWidth(0.15);
	problem.setChuteAngle(26.7);
	problem.setRoughBottomType(MONOLAYER_DISORDERED);

	//funr; // Funnel radius.
	//funO[2]; // (x,y) location of the center of the top of the funnel.
	//funa; // Angle of the funnel
	//funH; // Heigth of the funnel
	//funHf; // Falling heigth
	//funD; // Funnel diameter at the downside of the funnel.
	//funnz; //Number of particles along the heigth of the funnel
	//funfr; // Filling ratio of the funnel
	//fundiag; //The diagonal of the filling region
	//funrmax; //The maximum range for r
	problem.set_funa(60.0); // Angle of the funnel
	problem.set_funD(0.015); // Funnel diameter at the downside of the funnel.
	problem.set_funHf(0.05); // Falling heigth
	problem.set_funnz(25.0); //Number of particles along the heigth of the funnel
	//problem.set_funO(0, 0.5*(problem.getYMax()+problem.getYMin())); // (x,y) location of the center of the top of the funnel.
	problem.set_funfr(0.3); // Filling ratio of the funnel
	
	problem.setInflowVelocity(0);
	problem.setInflowVelocityVariance(0.01);
	problem.setMaxFailed(1);
	
	std::cout << "Chute inflow height: " << problem.getInflowHeight() << "Chute inflow velocity: " << problem.getInflowVelocity() << "Chute inflowvelocityvariance: " << problem.getInflowVelocityVariance() << std::endl;

	//Discretization parameters
	//problem.setHGridMaxLevels(1);
	//problem.setHGridNumberOfBucketsToPower(1e6); //automated
	double mass = species->getMassFromRadius(0.5 * (problem.getMinInflowParticleRadius() + problem.getMaxInflowParticleRadius()));
    problem.setTimeStep(0.02 * species->getCollisionTime(mass));
	problem.setTimeMax(.5e-2);
	problem.setSaveCount(10);

	std::cout << " dt = " << problem.getTimeStep() << " tmax = " << problem.getTimeMax() << std::endl;
	//	std::cout << "Maximum allowed speed of particles: " << problem.getMaximumVelocity() << endl; // speed allowed before particles move through each other!


	problem.readArguments(argc, argv);
	problem.set_funO(-0.9*problem.get_funHf()*sin(problem.getChuteAngle()), 0.5*(problem.getYMax()+problem.getYMin())); // (x,y) location of the center of the top of the funnel.
	problem.setName_();
	problem.solve();
}
