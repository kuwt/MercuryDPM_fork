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

#include <sstream>
#include <iostream>


//This is only used to change the file permission of the xball script create, at some point this code may be moved from this file to a different file.
#include <sys/types.h>
#include <sys/stat.h>

#include "DPMBase.h"

/// This problem consists of a ball sitting on the bottom of a gravity-bound box with horizontal velocity. Tangential dissipation is turned on, so the ball gradually looses velocity and gains angular velocity due to bottom friction, until the relative velocity reaches zero.
class PROBLEM : public DPMBase { 
						
	// initialise particle position, velocity, radius
	void setupInitialConditions()
	{
		set_N(1);
		particleHandler.getObject(0)->Position=Vec3D(0.005,0.0004,0.0);
		particleHandler.getObject(0)->setVelocity(Vec3D(0.0,0.0,0.0));
		particleHandler.getObject(0)->Angle=Vec3D(0.0,0.0,0.0);
		particleHandler.getObject(0)->setAngularVelocity(Vec3D(0.0,0.0,0.0));
		particleHandler.getObject(0)->getRadius()=0.0004;
		
		// Walls
		set_NWall(1);
		Walls[0].set(Vec3D(0,-1,0), -getYMin());

		// Periodic Walls
		set_NWallPeriodic(1);
		WallsPeriodic[0].set(Vec3D(1,0,0), getXMin(), getXMax());
	}
	
};


int main(int argc UNUSED, char *argv[] UNUSED)
{
	///Start off my solving the default problem
	PROBLEM problem;
	problem.setName("rotation");
	problem.setSystemDimensions(2);
	problem.setTimeMax(1.0);
	
	problem.getSpecies(0)->setDissipation(0.0);
	problem.setParticleDimensions(3);
	problem.speciesHandler.getObject(0)->setSlidingFrictionCoefficient(tan(20.0*constants::pi/180.0));
	problem.speciesHandler.getObject(0)->setSlidingDissipation(1e-2);
	problem.speciesHandler.getObject(0)->setStiffness(1e2);
	problem.setTimeStep(0.02*helpers::getCollisionTimeFromKAndDispAndEffectiveMass(problem.getSpecies(0)->getStiffness(), problem.getSpecies(0)->getDissipation(), 0.5*problem.getSpecies(0)->getMassFromRadius(0.0004)));
	problem.setGravity(Vec3D(0.01,-1.0,0.0));
	problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimestep(100,problem.getTimeMax(),problem.getTimeStep()));
	problem.solve();
	std::cout << problem;
	// std::cout << problem;
}
