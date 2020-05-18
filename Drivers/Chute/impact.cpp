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

#include "scr/ChuteWithHopper.h"

#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

class GranularJet : public ChuteWithHopper {
public:

	//void actionsBeforeTimeLoop(){write(std::cout,false); save_info_to_disk();}
	void computeExternalForces(int CI)
	{
		
		/// Now add on gravity
		if (Particles[CI].Position.Z>getHopperLift()+hopperExitHeight_)
			Particles[CI].Force += MassFlowFactor * gravity * Particles[CI].get_mass();
		else
			Particles[CI].Force += gravity * Particles[CI].get_mass();
		
		///Finally walls
		if (!Particles[CI].is_fixed()) computeForceDueToWalls(CI);
		
	}
	
	double InitialVelocity;
	double MassFlowFactor;

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	GranularJet problem;
	
	// Problem parameters
	problem.setName("impact");
	//Should be 10 for full length problem, but here I keep it low for a test case
	problem.setTimeMax(25);
	
	
	// Particle properties
	problem.setDensity(2400.0);
	
	problem.setInflowParticleRadius(0.5,0.5);
	problem.speciesHandler.getObject(0)->setCollisionTimeAndRestitutionCoefficient(4e-3, 0.6);
	problem.speciesHandler.getObject(0)->setSlidingDissipation(problem.get_dissipation());
	problem.speciesHandler.getObject(0)->setSlidingFrictionCoefficient(0.8);
	problem.setFixedParticleRadius(0.0);
	problem.setRoughBottomType(MONOLAYER_ORDERED);
	problem.setAlignBase(false);
	
	// Chute properties
	problem.setChuteAngle(0.0);
	problem.setChuteLength(100.0);
	problem.setChuteWidth(100.0);
	problem.setMaxFailed(6);
	//problem.makeChutePeriodic();
	problem.setHopperDimension(2);
	problem.setHopperLift(20);
	double ExitHeight = 8.0, ExitLength = 1.0 * ExitHeight, hopperAngle_ = 30.0, hopperLength_ = 6.0 * ExitLength;
	problem.setIsHopperCentred(true);
	problem.setHopperLowerFillingHeight(.75);
	problem.setAlignBase(false);
	problem.set_Hopper(ExitLength,ExitHeight,hopperAngle_,hopperLength_);
	problem.MassFlowFactor=20;
	//solve
	cout << "Maximum allowed speed of particles: " << problem.getMaximumVelocity() << endl; // speed allowed before particles move through each other!
	problem.setTimeStep();
	//problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(100,problem.getTimeMax(),problem.getTimeStep()));
	problem.setSaveCount(500);
	cout << "dt=" << problem.getTimeStep() << endl;
	
	problem.setXBallsAdditionalArguments("-sort -v0 -solidf -v0 -oh -200 -p 20");
	
	problem.auto_number();
	problem.readArguments(argc, argv);
	problem.solve();
	problem.write(cout);
	//problem.HGRID_base::write(cout);
	
	
	//cout << problem << endl;
	problem.writeRestartFile();
}
