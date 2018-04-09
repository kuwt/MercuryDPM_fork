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

#include<iostream>
#include "scr/Mercury3D.h"


/// In this file, 5 Particles (4 fixed) are loaded from files "../demos/5Particles.ini" and "../demos/5Particles.restart". The particles are aligned such that the single nonfixed particle rotates sinusoidally without moving. This is to test the behaviour of the tangential spring and the file loading routines.
class FiveParticles_restart : public Mercury3D{
public:

	void setupInitialConditions() {
 		//loads data file "../demos/5Particles.ini"
		if (!readDataFile("../demos/5Particles.ini")) {
			std::cerr << "Input data not found exiting " << std::endl;
			exit(-1);
		}
		//write(std::cout,true);
		getObjects().rbegin()->getAngularVelocity().Y /= 2;
		setTimeMax(0.00829999);
	}

};

int main(int argc UNUSED, char *argv[] UNUSED)
{

 	FiveParticles_restart problem;
 	//loads restart file "../demos/5Particles.restart"
 	problem.set_name("../demos/5Particles");
	problem.load_restart_data();
	
		
	const int imax=10;
	double a[imax], w[imax], dt[imax];
	for (int i=0; i<imax; i++) {
		stringstream ss("");
		ss << "5Particles_convergence" << i;
		problem.set_name(ss.str().c_str());
		problem.setTimeStep(0.5*problem.getTimeStep());
		problem.setSaveCount(2*problem.get_savecount());
		problem.solve();
		dt[i] = problem.getTimeStep();
		a[i] = problem.getObjects().rbegin()->Angle.Y;
		a[i] = problem.getObjects().rbegin()->TangentialSprings[0].delta.GetLength();
		w[i]= problem.getObjects().rbegin()->getAngularVelocity().Y;		
		std::cout << " t: " << setw(10) << setprecision(10) << problem.getTime() 
			 << ", a:" << setw(10) << setprecision(10) << problem.getObjects().rbegin()->Angle.Y
			 << ", w:" << setw(10) << setprecision(10) << problem.getObjects().rbegin()->getAngularVelocity().Y << std::endl;
		problem.writeRestartFile();
	}
	for (int i=1; i<imax-2; i++) {
		std::cout 
		 << " dt=" << setw(10) << setprecision(4) << dt[i]
		 << ", a=" << setw(12) << setprecision(6) << a[i]
		 << ", w=" << setw(12) << setprecision(6) << w[i]
		 << ", c_a=" << setw(8) << setprecision(2) << log((a[i]-a[imax-1])/(a[i-1]-a[imax-1])) / log(dt[i-1]/dt[i])
		 << ", c_w=" << setw(8) << setprecision(2) << log((w[i]-w[imax-1])/(w[i-1]-w[imax-1])) / log(dt[i-1]/dt[i])
		 << std::endl;
	}
}
