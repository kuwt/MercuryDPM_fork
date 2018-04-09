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

#include "scr/Mercury3D.h"
#include "scr/StatisticsVector.h"
#include <cmath>
#include <iostream>
#include <iomanip>


//creates fcc-layers
class HexagonalPacking : public Mercury3D {

public:

	void setupInitialConditions()
	{
		// make sure the number of particles is right
		set_N(4);
		setXMax(1);
		setYMax(1);
		setZMax(5);
		get_ParticleHandler().get_Particle(0)->getRadius() = 0.25;
		get_ParticleHandler().get_Particle(1)->getRadius() = 0.5;
		particleHandler.getObject(2)->Radius = 0.5;
		Particles[3].Radius = 0.5;
		
		get_ParticleHandler().get_Particle(0)->Position = Vec3D(.5,.5,0.75);
		get_ParticleHandler().get_Particle(1)->Position = Vec3D(.5,.5,1.5);
		particleHandler.getObject(2)->Position = Vec3D(.5,.5,3.5);
		Particles[3].Position = Vec3D(.5,.5,4.5);
		
		get_ParticleHandler().get_Particle(0)->fixParticle();
		get_ParticleHandler().get_Particle(1)->getVelocity() = Vec3D(0,0,0);
		particleHandler.getObject(2)->fixParticle();
		Particles[3].Velocity = Vec3D(0,0,0);
	}
	
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main2(int argc UNUSED, char *argv[] UNUSED)
{
	//~ HexagonalPacking problem;
	//~ problem.set_name("HexagonalPacking");
	//~ 
	//~ //set the number of particles
	//~ problem.set_N(mathsFunc::cubic(5));
	//~ problem.setSystemDimensions(3);
	//~ problem.setDensity(1.9098593);
	//~ problem.setGravity(Vec3D(0.,0.,-.1));
	//~ problem.speciesHandler.getObject(0)->setCollisionTimeAndRestitutionCoefficient(1.0,0.2,1.);
	//~ problem.setTimeStep(.02);
	//~ problem.setTimeMax(200.0);
	//~ problem.setSaveCount(250);
	//~ std::cout << "Relax the packing" << std::endl;
	//~ problem.solve();
	//~ 
	//~ std::cout << "Get statistics" << std::endl;
	//~ StatisticsVector<Z> stats3("HexagonalPacking");
	//~ double n =300;
	//~ stats3.setN(n);
	//~ stats3.setCGWidth(problem.getXMax()/n/2.);
	//~ stats3.statistics_from_fstat_and_data();
}

int main(int argc UNUSED, char *argv[] UNUSED)
{

	
	double dim = 3;
	double numcoeff = 5;
	double lucycoeff[] = {-3,8,-6,0,1};
	NORMALIZED_POLYNOMIAL poly(lucycoeff, numcoeff, dim);
	// (16 * pi) / 105 = 0.478718881
	
	numcoeff = 1;
	double heavisidecoeff[] = {1};
	NORMALIZED_POLYNOMIAL poly2(heavisidecoeff, numcoeff, dim);
	// 	(8 * pi) / 6 = 4.1887902
	return 0;
}
