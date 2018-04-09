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
#include <iomanip>
#include <cmath>

#include "scr/ChuteWithHopper.h"
using namespace std;

///This code does the MD of a normal shock into a wall, it has 
///autorestart protection, which is currently under testing.

class AirySavageHutter : public ChuteWithHopper
{
	
	public:
	
	void setupInitialConditions()
		{
			cout << "Entering the solve now what happens" << endl;
			cout << "Problem name " <<getName() <<endl;
			
			
			
			if (load_restart_data()==1)
				{
			
				
			
				write(std::cout,false);
	
				if (getTimeMax()>getTime())
					{
					cout << "Problem is not complete, will restart shortly" << endl;
					}
			
				else
					{
					cout << "Problem has been run and is complete: About to quit" << endl;
			
					exit(0);
					}
				}
				//end if load restart date true
			else 
			//If you are in here it is a new problem and about to be started
				{
	
				//Setup the base i.e. the chute particles
				ChuteWithHopper::setupInitialConditions();
		
		
				//Set up the walls
				int Nwall=get_NWall();
				set_NWall(Nwall+1);
				//Put a wall one particle dimater from the end of the chute, this stops gaps
				Walls[Nwall].set(Vec3D(1.0, 0.0, 0.0), getXMax()-2.0*getFixedParticleRadius ());
		
		
		
				}
		}
	
	
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	AirySavageHutter problem;
	

	
	// Problem parameters
	problem.setName("AirySavageHutter");
	problem.setTimeStep(1e-4);
	problem.setTimeMax(500.0);
	problem.set_HGRID_max_levels(2);
	
	// Particle properties
	problem.setInflowParticleRadius(0.5);
	problem.setFixedParticleRadius(0.25);
	problem.setRoughBottomType(MULTILAYER);
	problem.setDensity(6/pi);
	problem.speciesHandler.getObject(0)->setStiffness(2e5);
	problem.speciesHandler.getObject(0)->setSlidingStiffness(2.0/7.0*problem.speciesHandler.getObject(0)->getStiffness());
	problem.setDissipation(25.0);
	problem.speciesHandler.getObject(0)->setSlidingDissipation(0);
	problem.speciesHandler.getObject(0)->setSlidingFrictionCoefficient(0.5);
 	
	// Chute properties
	problem.setChuteAngle(27.0,1.0);
	problem.setChuteLength(250);
	problem.setChuteWidth(10);
	problem.setMaxFailed(6);
	//problem.makeChutePeriodic();
	
	//Hopper properties
	double ExitHeight = 10.0, ExitLength = 10.0, hopperAngle_ = 45.0, hopperLength_ = 2.0 * ExitLength;
	problem.set_Hopper(ExitLength,ExitHeight,hopperAngle_,hopperLength_);
	problem.setHopperFillPercentage(50.0);
	
	// Inflow values
	//problem.setInflowHeight(10);
	//problem.setInflowVelocity(0.7);
	//problem.setInflowVelocityVariance(0.02);
	
	// Xballs tunning
	problem.setXBallsScale(0.012);
	problem.setXBallsAdditionalArguments("-sort -v0 -oh 500 -cmode 4");
	
	//solve
	cout << "Maximum allowed speed of particles: " << problem.getMaximumVelocity() << endl; // speed allowed before particles move through each other!
	//You really do not need data more often that once a second do not change this number away from this without really good reason.
 	problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimestep(500,problem.getTimeMax(),problem.getTimeStep()));
	cout << "dt=" << problem.getTimeStep() << endl;
	
	problem.auto_number();
	problem.solve();
	//problem.writeRestartFile();
}
