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

#include "SilbertPeriodic.h"
#include "scr/ChuteBottom.h"

using namespace std;

class FlowRule : public SilbertPeriodic {
public:

	void createBottom() {
		//create random bottom
		ChuteBottom bottom(*this);
		bottom.setThickness(1.0);
		bottom.set_periodicbottom(false);
		bottom.setInflowParticleRadius(getFixedParticleRadius());
		bottom.make_rough_bottom(Particles);

		//add particles
		hGridActionsBeforeTimeLoop();
		hGridActionsBeforeTimeStep();
		//HGridActionsBeforeTimeStep();
		int failed = 0, max_failed = 1000000;
		//try max_failed times to find new insertable particle
		while (failed<=max_failed){
			P0.Radius = FixedParticleRadius;
			P0.computeMass(Species);
			P0.Position.X = random.getRandomNumber(P0.Radius, xmax-P0.Radius);
			P0.Position.Y = random.getRandomNumber(P0.Radius, getYMax()-P0.Radius);
			P0.Position.Z = random.getRandomNumber(-bottom.getThickness()*getFixedParticleRadius(), -0.0);
			P0.Velocity.set_zero();
			if (IsInsertable(P0)) 
			{
				cout << "Particle added" << endl;
				failed = 0; 
				num_created++;
			} 
			else failed++;
		}

		//finally, fix particles to the floor
		for (vector<CParticle>::iterator it= this->Particles.begin(); it!=this->Particles.end(); ++it) 
			it->fixParticle();

		//now run it quickly
		setStiffness(2e5);
		setSlidingStiffness(2.0/7.0*getStiffness());
		setDissipation(25.0);
		setSlidingDissipation(getDissipation());

		Mdouble dist2=1e20;
		for (unsigned int i=0; i<particleHandler.getNumberOfObjects(); i++) 
		{
			for (int j=i+1; j<particleHandler.getNumberOfObjects(); j++) 
			{
				dist2 = min(dist2, GetDistance2(particleHandler.getObject(i)->getPosition(), Particles[j].Position));
			}
		}
		cout << "Dist" << setprecision(10) << dist2 << endl;
	};
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	FlowRule problem;
	//problem.randomize();
	problem.setName("bottom");
	problem.setTimeMax(0.0);
	problem.setSaveCount(5000);
	problem.setXMax(60);
	problem.setYMax(80);
	problem.setFixedParticleRadius(.502);
	problem.set_H(2e-4);
	
	problem.solve();
}
