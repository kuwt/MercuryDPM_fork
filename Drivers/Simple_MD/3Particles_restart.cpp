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
#include "Mercury3D.h"


/** In this file, 3 Particles (2 fixed) are loaded from files 
 * "../demos/3Particles.ini" and "../demos/3Particles.restart". 
 * The particles are aligned such that the single nonfixed particle 
 * moves sinusoidally in x-direction. This is to test the behaviour 
 * of the normal spring and the file loading routines. **/
class ThreeParticles_restart : public Mercury3D{
public:

	void setupInitialConditions() {
 		//loads data file "../demos/5Particles.ini"
		if (readDataFile("../demos/3Particles.ini")) {
			//write(std::cout,false);
			write(std::cout,true);
		} else {
			std::cerr << "Input data not found exiting " << std::endl;
			exit(-1);
		}
	}

};

int main(int argc UNUSED, char *argv[] UNUSED)
{

 	ThreeParticles_restart problem;
 	//loads restart file
 	problem.set_name("../demos/3Particles");
	problem.load_restart_data();
	 	
	problem.solve();
	
	problem.write(std::cout,false);
	problem.writeRestartFile();
}

/* 3 particles of radius .5 in a domain of volume 2.5
 * >> r=.5; Vp=4/3*constants::pi*r^3; Nu=3*Vp/2.5
 * Nu = 0.6283 */
