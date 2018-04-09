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

class spring : public DPMBase {

	void setupInitialConditions()
	{
		//dimensions
		setSystemDimensions(2);
		setParticleDimensions(3);

		//box size
		setXMax(2);
		setYMax(2);
		setZMax(1);
		
		set_NWall(4);
		Walls[0].set(Vec3D(1.0,0.0,0.0),getXMax());
		Walls[1].set(Vec3D(-1.0,0.0,0.0),-getXMin());
		Walls[2].set(Vec3D(0.0,1.0,0.0),getYMax());
		Walls[3].set(Vec3D(0.0,-1.0,0.0),-getYMin());

		//time stepping
		setTimeStep(1e-4);
		setTimeMax(10000*getTimeStep());
		setSaveCount(10);
	 
		//particle properties
		setDensity(6/constants::pi);
		setCollisionTimeAndRestitutionCoefficient(5e-3,0.88,1.0);
		//setStiffness(2e5);
		//setDissipation(25.0);
		setSlidingStiffness(2.0/7.0*getStiffness());
		setSlidingDissipation(0);
		setSlidingFrictionCoefficient(0.5);
		
		//Particles
		set_N(2);

		particleHandler.getObject(0)->Position=Vec3D(0.50, 1.0, 0.0);
		particleHandler.getObject(1)->Position=Vec3D(1.50, 1.0, 0.0);
	
		particleHandler.getObject(0)->setVelocity(Vec3D( 0.4, 0.1, 0.0));
		particleHandler.getObject(1)->setVelocity(Vec3D(-0.4,-0.1, 0.0));
	
		particleHandler.getObject(0)->getRadius()=0.5;
		particleHandler.getObject(1)->getRadius()=0.5;

		//~ set_N(3);
//~ 
		//~ particleHandler.getObject(0)->Position=Vec3D(0.50, 1.0, 0.0);
		//~ particleHandler.getObject(1)->Position=Vec3D(1.20, 0.5, 0.0);
		//~ particleHandler.getObject(2)->Position=Vec3D(1.50, 1.5, 0.0);
	//~ 
		//~ particleHandler.getObject(0)->setVelocity(Vec3D( 0.4, 0.1, 0.0));
		//~ particleHandler.getObject(1)->setVelocity(Vec3D(-0.4,-0.1, 0.0));
		//~ particleHandler.getObject(2)->Velocity=Vec3D(0.4,-0.1, 0.0);
	//~ 
		//~ particleHandler.getObject(0)->getRadius()=0.5;
		//~ particleHandler.getObject(1)->getRadius()=0.5;
		//~ particleHandler.getObject(2)->Radius=0.5;

		
	}

};

int main(int argc UNUSED, char *argv[] UNUSED)
{

	///Start off by solving the default problem
	spring problem;
	problem.set_name("springs");
	problem.setGravity(Vec3D(0.0,0.0,0.0));
	
	problem.solve();
	problem.write(std::cout,true);
}
