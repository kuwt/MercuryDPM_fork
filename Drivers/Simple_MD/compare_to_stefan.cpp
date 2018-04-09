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

///This code can be tested against Stefan's code
///\todo{It tests reasonably well against Stefans code with position verlet, but we have to test against velocity verlet.}

#include "DPMBase.h"
#include <iostream>
enum Spec {normal_spring, normal_spring_damper, normal_spring_damper_tangential_friction_damper, normal_spring_damper_tangential_friction_damper_spring};

class myMD : public DPMBase 
{
public:
	Spec species;
	
	void setupInitialConditions()
	{
		// set species comparable to stefan's species.ini
		if (species==normal_spring) {
			std::cout << "normal spring: " << std::endl;
			setStiffness(5.0);
			setSlidingStiffness(0.0);
			setSlidingFrictionCoefficient(0.0);
			set_dissipation(0.0);
			setSlidingDissipation(0.0);
		} else if (species==normal_spring_damper) {
			std::cout << "normal spring and damper: " << std::endl;
			setStiffness(5.0);
			setSlidingStiffness(0.0);
			setSlidingFrictionCoefficient(0.0);
			set_dissipation(1e-5);
			setSlidingDissipation(0.0);
		} else if (species==normal_spring_damper_tangential_friction_damper) {
			std::cout << "normal spring and damper, tangential friction with damper: " << std::endl;
			setStiffness(5.0);
			setSlidingStiffness(0.0);
			setSlidingFrictionCoefficient(0.5);
			set_dissipation(1e-5);
			setSlidingDissipation(1e-5);
		} else if (species==normal_spring_damper_tangential_friction_damper_spring) {
			std::cout << "normal spring and damper, tangential friction with damper and spring: " << std::endl;
			setStiffness(5.0);
			setSlidingStiffness(5.0);
			setSlidingFrictionCoefficient(0.5);
			set_dissipation(1e-5);
			setSlidingDissipation(1e-5);
		} else { std::cerr << "error: species not known" << std::endl;}

		// set variables comparable to stefans par.ini
		setGravity(Vec3D(0.0,0.0,-9.8));
		setSystemDimensions(3);
		setParticleDimensions(3);
		setTimeStep(1e-7);
		setSaveCount(100);
		setTimeMax(0.004+1.5*getTimeStep());
		setXMax(0.05);
		setYMax(0.05);
		setZMax(0.05);
		set_NWall(6);
		Walls[0].set(Vec3D(-1.0, 0.0, 0.0), -getXMin());
		Walls[1].set(Vec3D( 1.0, 0.0, 0.0),  getXMax());
		Walls[2].set(Vec3D( 0.0,-1.0, 0.0), -getYMin());
		Walls[3].set(Vec3D( 0.0, 1.0, 0.0),  getYMax());
		Walls[4].set(Vec3D( 0.0, 0.0,-1.0), -getZMin());
		Walls[5].set(Vec3D( 0.0, 0.0, 1.0),  getZMax());
		setDensity(2000);
		
		// set variables comparable to stefans c3d.ini
		set_N(2);
		// .004 .01 .01  0 0 1  .005 0
		particleHandler.getObject(0)->Position=Vec3D(.01, .01, .01);
		particleHandler.getObject(0)->setVelocity(Vec3D(0, 0, 0));
		particleHandler.getObject(0)->getRadius()=.005;
		particleHandler.getObject(1)->Position=Vec3D(.018, .01, .01);
		particleHandler.getObject(1)->setVelocity(Vec3D(0, 0, 1));
		particleHandler.getObject(1)->getRadius()=.005;
		
		fstat_file.precision(12);
		data_file.precision(10);
	}

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	myMD problem;
	problem.set_name("1");
	problem.species = normal_spring;
	problem.solve();

	myMD problem2; 
	problem.set_name("2");
	problem.species = normal_spring_damper;
	problem.solve();

	myMD problem3; 
	problem3.set_name("3");
	problem3.species = normal_spring_damper_tangential_friction_damper;
	problem3.solve();
	//std::cout << problem3 << std::endl;	

	///\todo{for some reason the results change if I reuse problem}
	myMD problem4; 
	problem2.set_name("4");
	problem2.species = normal_spring_damper_tangential_friction_damper_spring;
	problem2.solve();
	//std::cout << problem2 << std::endl;	
}
