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
using namespace std;

class VariableBottom : public ChuteWithHopper {
public:
	
	void setupInitialConditions()
	{
		cout << "restarted " << restarted << endl;
		if (!restarted) {
			set_NWallPeriodic(1);
			set_NWall(0);
			WallsPeriodic[0].set(Vec3D( 0.0, 1.0, 0.0), getYMin(), getYMax());
			add_hopper();	
		}
		set_HGRID_num_buckets_to_power();
		write(cout);
	}

	void cleanChute() 
	{
		//clean outflow every 100 time steps
		static int count = 0, maxcount = 100;
		if (count>maxcount)
		{
			count = 0;
			// delete all outflowing particles
			for (unsigned int i=0;i<Particles.size();)
			{
				if (Particles[i].Position.Z<-10*Particles[i].Radius)
				{
					#ifdef DEBUG_OUTPUT_FULL
						cout << "erased:" << Particles[i] << endl;
					#endif
					removeParticle(i);
				}	
				else i++;
			}
		} else count++;
	}


	///sets parameters of particles and time stepping to the L3 type used in Silbert's papers
	void set_silbert_parameters() {
		//time stepping
		if (true) {
			setTimeStep(1e-4*10);
			setTimeMax(2e6*getTimeStep());
			setStiffness(2e5/10);
			setSlidingStiffness(2.0/7.0*getStiffness());
		} else {
			setTimeStep(1e-4);
			setTimeMax(2e7*getTimeStep());
			setStiffness(2e5);
			setSlidingStiffness(2.0/7.0*getStiffness());
		}
		//particle properties
		setInflowParticleRadius(.5);
		setFixedParticleRadius(.5);
		setDensity(6/pi);
		setDissipation(25.0);
		setSlidingDissipation(0);
		setSlidingFrictionCoefficient(0.5);
	}

	///sets bottom of the chute
	void set_chute_parameters() {
		setZMax(30);

		unsigned int NFirstBottom = 40;
		vector<CParticle> FirstBottom;
		if (readDataFile("../ini_files/bottom0.25.ini")) {
			FirstBottom.resize(NFirstBottom*Particles.size());
			for (unsigned int j=0; j<NFirstBottom; j++)
			for (unsigned int i=0; i<Particles.size(); i++) {
				FirstBottom[j*Particles.size()+i]=Particles[i];
				FirstBottom[j*Particles.size()+i].Position.X += j*20;
				FirstBottom[j*Particles.size()+i].Position /= 2.0;
			}
			setChuteLength(20*NFirstBottom);
		} else {
			cerr << "Input data not found exiting " << endl;
			exit(-1);
		}
		unsigned int NSecondBottom = 20;
		vector<CParticle> SecondBottom;
		if (readDataFile("../ini_files/bottom1.0.ini")) {
			SecondBottom.resize(NSecondBottom*Particles.size());
			for (unsigned int j=0; j<NSecondBottom; j++)
			for (unsigned int i=0; i<Particles.size(); i++) {
				SecondBottom[j*Particles.size()+i]=Particles[i];
				SecondBottom[j*Particles.size()+i].Position.X += (NFirstBottom+j)*20;
				SecondBottom[j*Particles.size()+i].Position /= 2.0;
			}
			setChuteLength(2*getChuteWidth()*(NFirstBottom+NSecondBottom));
		} else {
			cerr << "Input data not found exiting " << endl;
			exit(-1);
		}
		set_Nmax(FirstBottom.size()+SecondBottom.size()+getChuteLength()*getChuteWidth()*getZMax());
		Particles.resize(FirstBottom.size()+SecondBottom.size());
		for (unsigned int i=0; i<FirstBottom.size(); i++) {
			Particles[i]=FirstBottom[i];
			Particles[i].fixParticle();
		}
		for (unsigned int i=0; i<SecondBottom.size(); i++) {
			Particles[FirstBottom.size()+i]=SecondBottom[i];
			Particles[FirstBottom.size()+i].fixParticle();
		}

		//chute properties
		setChuteWidth(5);
		setChuteLength(2*getChuteWidth()*(NFirstBottom+NSecondBottom));
		setChuteAngle(32.0, 1.0);
	}
		
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	VariableBottom problem;
	problem.setName("VariableBottomQuick");
	//problem.makeChutePeriodic();
	problem.set_chute_parameters();
	//problem.load_restart_data("../VariableBottomQuick.3/VariableBottomQuick.restart");
	problem.set_silbert_parameters();
	problem.setSaveCount(1000);

	//Hopper properties
	double ExitHeight = 20.0, ExitLength = 20.0, hopperAngle_ = 45.0, hopperLength_ = 4.0 * ExitLength;
	problem.set_Hopper(ExitLength,ExitHeight,hopperAngle_,hopperLength_);
	problem.setHopperFillPercentage(50.0);
	problem.setMaxFailed(40);
	problem.readArguments(argc, argv);
	problem.solve();
}
