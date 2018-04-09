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


#include "DPMBase.h"
#include <iostream>

class finite_walls : public DPMBase{
public:
	
	void computeInternalForces(int PI UNUSED, int PJ UNUSED) {}
	
	void computeExternalForces(int CI) {computeWalls(CI);}
	
	void setupInitialConditions()
	{
/*		set_N(1);
		particleHandler.getObject(0)->getRadius()=0.005;
		particleHandler.getObject(0)->Position=Vec3D(0.1001,0.12,0.0);
		particleHandler.getObject(0)->setVelocity(Vec3D(0.0,-0.1,0.0));*/
		double radius = 0.05;
		double dist = 0.05* 2.0 * radius;
		int N = floor(getXMax()/dist);
		set_N(N*N);
		int id = 0;
		for (int i=1; i<=N; i++) {
			for (int j=1; j<=N; j++) {
				Particles[id].Radius=radius;
				Particles[id].Position=Vec3D(i*dist,j*dist,0.0);
				Particles[id].Velocity=Vec3D(0.0,0.0,0.0);
				id++;
			}
		}

		set_NWall(1);
		Walls[0].addObject(Vec3D( 1.0,-1.0, 0.0), 0.0);
		Walls[0].addObject(Vec3D(-1.0,-2.0, 1.0),-0.2);
		Walls[0].addObject(Vec3D(-1.0,-2.0,-1.0),-0.2);
	}
	
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	
	///Start off my solving the default problem
	finite_walls problem;
	problem.set_name("finite_walls");
	problem.setTimeMax(2e-4);
	problem.setTimeStep(1e-4);
	problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimestep(1000,problem.getTimeMax(),problem.getTimeStep()));
	//problem.setSaveCount(1);
	problem.speciesHandler.getObject(0)->setCollisionTimeAndRestitutionCoefficient(5e-3,1.0,0.157);
	problem.setXMax(0.2);
	problem.setYMax(0.2);
	problem.setZMax(0.2);
	
	problem.solve();
	
	//now run ./finite_walls/finite_walls.disp -rmult .05
}
