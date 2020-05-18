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
#include <sys/types.h>
#include <sys/stat.h>
#include "scr/Chute.h"
#include "scr/Time.h"

using namespace std;
bool quick = true;

class VariableBottom : public Chute
{
public:
	void setupInitialConditions(){
		createBottom();
	}	

	///sets parameters of particles and time stepping to the L3 type used in Silbert's papers
	void set_silbert_parameters() {
		setInflowParticleRadius(.5);
		setFixedParticleRadius(.5);
		setDensity(6/pi);
		//time stepping
		if (quick) {
			setTimeStep(1e-3);
			setCollisionTimeAndRestitutionCoefficient(0.05, .88);
		} else {
			setTimeStep(1e-4);
			setCollisionTimeAndRestitutionCoefficient(0.005, .88);
		}
		//particle properties
		setSlidingStiffness(2.0/7.0*getStiffness());
		//setSlidingDissipation(0.0);
		setSlidingDissipation(getDissipation());
		setSlidingFrictionCoefficient(0.5);
	}

	///sets parameters of the chute
	void set_chute_parameters() {
		//chute properties
		setChuteAngle(24, 1.0);
		setChuteWidth(10);
		setChuteLength(80);
		setZMax(30);
		setTimeMax(1e20);
		setSaveCount(1e4);
		//create a wall below the chute
		set_NWall(0);
		//~ Walls[0].set(Vec3D( 0.0, 0.0, -1.0), 3.4*MaxInflowParticleRadius);
		//create periodic walls
		set_NWallPeriodic(2);
		setXMax(20*max(1.,round(getXMax()/20)));
		WallsPeriodic[0].set(Vec3D( 0.0, 1.0, 0.0), getYMin(), getYMax());
		WallsPeriodic[1].set(Vec3D( 1.0, 0.0, 0.0), getXMin(), getXMax());
	}

	void createBottom() {
		unsigned int NFirstBottom = round(getXMax()/20);
		DPMBase FirstBottom;
		if (!FirstBottom.readDataFile(quick?"../ini_files/H20A24L0.5M0.5Quick.ini":"../ini_files/H20A22L0.5M0.5.ini",14)) {
			cerr << "1st input data not found exiting " << endl;
			exit(-1);
		}
		set_Nmax(FirstBottom.get_N()*NFirstBottom
			+getChuteLength()*getChuteWidth()*getZMax());
		set_N(0);
		for (int j=0; j<FirstBottom.get_N(); j++) {
			if (FirstBottom.getObjects()[j].Velocity.GetLength2()==0.0) 
				FirstBottom.getObjects()[j].fixParticle();
		}
		for (unsigned int i=0; i<NFirstBottom; i++)
		for (int j=0; j<FirstBottom.get_N(); j++) {
			Particles.push_back(FirstBottom.getObjects()[j]);
			Particles.back().Position.X += i*20;
		}
	}

	void actionsBeforeTimeStep(){}

	void printTime() const {
		static int Nold = get_N();
		static double told = getTime();
		cout << "t=" << setprecision(3) << left << setw(6) << getTime() 
			<< ", tmax=" << setprecision(3) << left << setw(6) << getTimeMax()
			<< ", N=" << setprecision(3) << left << setw(6) << Particles.size()
			//~ << ", dN/dt=" << setprecision(3) << left << setw(6) << (get_N()-Nold)/(getTime()-told)
			<< endl;
		Nold=get_N();
		told=getTime();
		//~ static unsigned int counter=0;
		//~ if (++counter>10) {counter=0; cout.flush();}
	}

};

int main(int argc, char *argv[])
{
	VariableBottom problem;
	problem.setName("LongChute");

	problem.set_silbert_parameters();
	problem.set_chute_parameters();
	problem.readArguments(argc, argv);
	problem.solve();
}
