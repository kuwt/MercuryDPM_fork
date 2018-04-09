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
#include "scr/DPMBase.h"
#include "scr/Mercury3D.h"
#include "scr/Chute.h"
using namespace std;

///This code examines the flow front of rough-bottom chute flow. The 
///flow is initialised on \f$x\in[0,FlowLength]\f$. Then the chute gradually 
///expands to the right as the flow develops, and is cut on the left 
///to minimize computation.
class FlowFrontChute : public Chute{
public:

	void actionsBeforeTimeStep(){};
		
	void setupInitialConditions() {};
	
	void stretch() {
		///prolong the chute 10-fold
		int stretchFactorBottom = 4;
		int stretchFactorFlow = 2;
		
		//count fixed particles
		int FixedParticles = 0;
		for (vector<CParticle>::iterator it = Particles.begin(); it!=Particles.end(); ++it) if (it->is_fixed()) {
			FixedParticles++;
		}
		
		set_Nmax(get_N()+(stretchFactorBottom-1)*FixedParticles+(stretchFactorFlow-1)*(get_N()-FixedParticles));
		//add new flow particles
		for (vector<CParticle>::reverse_iterator it = Particles.rbegin(); it!=Particles.rend(); ++it) if (!it->is_fixed()) {
			for (int i=1; i<stretchFactorFlow; i++) {
				Particles.push_back(*it);
				Particles.back().Position.X += i*(xmax-getXMin());
			}
		}
		//add new fixed particles
		for (vector<CParticle>::reverse_iterator it = Particles.rbegin(); it!=Particles.rend(); ++it) if (it->is_fixed()) {
			for (int i=1; i<stretchFactorBottom; i++) {
				Particles.push_back(*it);
				Particles.back().Position.X += i*(xmax-getXMin());
			}
		}
		
		//stretch domain
		setXMax(getXMin()+stretchFactorBottom*(xmax-getXMin()));
		for (vector<CWallPeriodic>::iterator it = WallsPeriodic.begin(); it!=WallsPeriodic.end(); ++it) {
			if (it->getNormal().X==1.0) it->set(Vec3D( 1.0, 0.0, 0.0), getXMin(), xmax);
		}
		
		//other settings
		set_HGRID_num_buckets_to_power();
		//setSaveCount(.2/dt);
		//setTimeMax(2000);
		
	}

};

int main(int argc UNUSED, char *argv[] UNUSED)
{

 	FlowFrontChute problem;
 	problem.setName("ini/silbert.theta.26.z.40");
	problem.load_restart_data();
	problem.stretch();
	 	
	problem.write(std::cout,false);
	problem.solve();
	
	problem.write(std::cout,false);
	problem.writeRestartFile();
}
