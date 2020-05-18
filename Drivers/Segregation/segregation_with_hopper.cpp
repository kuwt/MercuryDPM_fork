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

#include "ChuteWithHopper.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <Species/LinearViscoelasticFrictionSpecies.h>

class SegregationWithHopper : public ChuteWithHopper {
public:
	
	void create_inflow_particle()
	{
		
		int N= getNCreated();


		//This is 33% chance the particle is large
		int temp=(N%3)/2;

		//Difference in radius is sqrt 2 so factor of 2 in size.
		inflowParticle_.setRadius(0.3e-3*(1.0+(sqrt(2)-1.0)*temp));
		
		//inflowParticle_.computeMass();
		
		static double s = sin(getChuteAngle());
		static double c = cos(getChuteAngle());
		static double Ht = tan(getHopperAngle());
		static double Hc = cos(getHopperAngle());
		double gamma = random.getRandomNumber(0.5,1.0);
		// setting up particle position
		static Vec3D AB = Vec3D(c,0.0,s);
		AB = AB*(random.getRandomNumber(-1.0,1.0)*(0.5*gamma*getHopperLength() - inflowParticle_.getRadius()/Hc));
		static Vec3D AC = Vec3D(-s,0.0,c);
		AC = AC * (gamma*0.5*getHopperLength()/Ht);
		static Vec3D A = Vec3D(0.0, 0.0, s*(-0.5*(getHopperLength()-getHopperExitLength())) + c*getHopperHeight()) + AB*0.5*getHopperLength() + AC*(-0.5*getHopperLength()/Ht);
		static Vec3D AY = Vec3D(0.0,random.getRandomNumber(getYMin()+ inflowParticle_.getRadius(),getYMax()- inflowParticle_.getRadius()),0.0);
		//double gamma = random.get_RN(0.5,1.0);
		//P0.Position = A + AB * (random(-1.0,1.0)*(0.5*gamma*hopperLength_ - P0.Radius/Hc)) + AC * (gamma*0.5*hopperLength_/Ht);
		//P0.Position.Y = random(ymin+P0.Radius, ymax-P0.Radius);
		static Vec3D position = A + AB + AC + AY;
		inflowParticle_.setPosition(position);
	}
	
    int getNCreated() const
    {
        return nCreated_;
    }

    void increaseNCreated()
    {
        nCreated_++;
    }

    int nCreated_;
	SphericalParticle inflowParticle_;

};
int main()
{
	SegregationWithHopper problem;
	
	// Problem parameters
	problem.setName("segregation");
	//Should be 10 for full length problem, but here I keep it low for a test case
	problem.setTimeMax(10);
	
	
	// Particle properties
    auto species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    
    species->setDensity(2400.0);
	
	problem.setInflowParticleRadius(0.3e-3,0.60e-3);
    double mass = 0.5*species->getMassFromRadius(0.5*(problem.getMinInflowParticleRadius() + problem.getMaxInflowParticleRadius()));
    species->setCollisionTimeAndRestitutionCoefficient(4e-4,0.6,mass);
	species->setSlidingDissipation(species->getDissipation());
	species->setSlidingFrictionCoefficient(0.8);
	problem.setFixedParticleRadius(0.3e-3);
	problem.setRoughBottomType(MONOLAYER_DISORDERED);
	
	// Chute properties
	problem.setChuteAngle(25.0);
	problem.setChuteLength(600.0e-3);
	problem.setChuteWidth(3e-3);
	problem.setMaxFailed(6);
	problem.makeChutePeriodic();
	double ExitHeight = 12.0e-3, ExitLength = 1.0 * ExitHeight, hopperAngle_ = 60.0, hopperLength_ = 6.0 * ExitLength;
	Mdouble hopperLowestPoint_ = ExitHeight - ExitLength * tan(problem.getChuteAngle());
	Mdouble hopperHeight_=hopperLowestPoint_ + 1.1 * 0.5*(hopperLength_+ExitLength) / tan(hopperAngle_*constants::pi/180.0);
	Mdouble HopperCornerHeight = hopperHeight_ - 0.5*(hopperLength_-ExitLength) / tan(hopperAngle_*constants::pi/180.0);
	if (HopperCornerHeight<=0.0) 
	  { 
	    hopperHeight_ += -HopperCornerHeight + problem.getMaxInflowParticleRadius();
	    HopperCornerHeight = problem.getMaxInflowParticleRadius();
	  }
    problem.setHopper(ExitLength, ExitHeight, hopperAngle_, hopperLength_, hopperHeight_);
	
	//solve
	//std::cout << "Maximum allowed speed of particles: " << problem.particleHandler.getSmallestParticle()->calculateMaximumVelocity() << std::endl; // speed allowed before particles move through each other!
    problem.setTimeStep(0.02 * species->getCollisionTime(mass));
    problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(25, problem.getTimeMax(),problem.getTimeStep()));

    std::cout << "dt=" << problem.getTimeStep() << std::endl;
	
	problem.autoNumber();
	problem.write(std::cout,false);
	problem.writeRestartFile();
	
	//This set to colouring based of size and small vectors
	problem.setXBallsColourMode(7);
	problem.setXBallsVectorScale(1);
	problem.setXBallsScale(20.0);
	problem.setXBallsAdditionalArguments("-v0 -solidf");
	
	
	problem.solve();
	problem.write(std::cout);
	//problem.HGRID_base::write(std::cout);
	
	
	//std::cout << problem << std::endl;
	problem.writeRestartFile();
}
