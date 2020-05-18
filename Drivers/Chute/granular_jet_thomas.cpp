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

	void set_silbert() {
		//time stepping
		setTimeStep(1e-4);
		setTimeMax(1e20);
		setSaveCount(1e4); //save every unit time (\hat{t}=sqrt(d/g)=~0.008s)
	 
		//particle radii
		setInflowParticleRadius(.5*8 ); //.5
		setFixedParticleRadius(.5*16); //.5
		setRoughBottomType(MULTILAYER);

		//particle properties
		setDensity(6/pi);
		setStiffness(2e5);
		setSlidingStiffness(2.0/7.0*getStiffness());
		setDissipation(25.0);
		setSlidingDissipation(getDissipation());
		setSlidingFrictionCoefficient(0.5);
		
		//chute properties
		setChuteAngle(0.0, 1.0);
	}

	bool readNextArgument(unsigned int& i, char *argv[]){
		if (!strcmp(argv[i],"-ExitLength")) {
			double ExitHeight = atof(argv[i+1]), 
			       ExitLength = 1.0 * ExitHeight, 
			       hopperAngle_ = 45.0, 
			       hopperLength_ = 4.0 * ExitLength;
			set_Hopper(ExitLength,ExitHeight,hopperAngle_,hopperLength_);
		} else return Chute::readNextArgument(i, argv);
		return true; //returns true if argv is found
	}

	void actionsBeforeTimeLoop(){write(std::cout,false); save_info_to_disk();}
	

};

int main(int argc, char *argv[]) 
{
	GranularJet problem;

	// Particle properties
	problem.set_silbert();	
	// Problem parameters
	problem.setName("granular_jet");

	//Corrections
	problem.setRoughBottomType(MONOLAYER_DISORDERED); //to save base particles
	problem.setAlignBase(false);
	
	// Chute properties
	problem.setChuteAngle(20);
	problem.setChuteLength(700); //700=40cm
	problem.setChuteWidth(400); //400=24cm
	problem.setMaxFailed(6);
	problem.makeChutePeriodic();
	problem.setHopperDimension(2);
	problem.setHopperLift(250); //250=15cm
	double ExitHeight = 25.0, //25=15mm
	       ExitLength = 3.0 * ExitHeight, 
	       hopperAngle_ = 45.0, 
	       hopperLength_ = 4.0 * ExitLength;
	problem.set_Hopper(ExitLength,ExitHeight,hopperAngle_,hopperLength_);
	problem.setZMax(200);
	
	//solve
	problem.setTimeStep(); 
	problem.setSaveCount(5000);
	problem.dataFile.setFileType(FileType::ONE_FILE);
	problem.fStatFile.setFileType(FileType::NO_FILE);
	problem.restartFile.setFileType(FileType::ONE_FILE);
	problem.setXBallsAdditionalArguments("-v0 ");
	cout << "Maximum allowed speed of particles: " << problem.getMaximumVelocity() 
	<< endl; // speed allowed before particles move through each other!
	cout << "dt=" << problem.getTimeStep() << endl;
	
	problem.auto_number();

	//solve
	problem.readArguments(argc, argv);
	problem.solve();
}
