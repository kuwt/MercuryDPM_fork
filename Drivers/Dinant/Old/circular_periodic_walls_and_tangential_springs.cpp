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

//#include "../../src/Mercury3D.h"
#include "DPMBase.h"
#include <sstream>
#include <iostream>
#include <cstdlib>

//This test is a multi puprose test for:
//Periodic Particles in combination with HGRID (early a bug arose when a particle had a collision with a periodic particle, but where in the same cell)
//Tangential Springs over periodic contacts and during transitions from periodic to normal and reverse
//Removal of Particles while tangential spring information has te be conserved

//class circular_periodic_walls_and_tangential_springs : public Mercury3D {
class circular_periodic_walls_and_tangential_springs : public DPMBase
{
  
	void setupInitialConditions()
	{
		particleHandler.clear();
		boundaryHandler.clear();
		
		//dimensions
		setSystemDimensions(2);
		setParticleDimensions(2);

		//box size
		setXMin(-2);
		setYMin(-2);
		setXMax(2);
		setYMax(2);
		
		CircularPeriodicBoundary b0;
		get_BoundaryHandler().copyAndAddBoundary(b0);
		
		///Test case 2, Particle moving inward and outward agian
		/*Particle P0;
		setTimeMax(16.0);
		P0.setPosition(Vec3D(8.5 ,0.5  , 0.0));
		P0.setVelocity(Vec3D(-1.0, 0.0, 0.0));
		P0.setRadius(0.1);
		particleHandler.copyAndAddObject(P0);*/
		
		///Test case 2, Particles rotating (note external forces have to be enabled)
		/*Particle P0;
		setTimeMax(16.0); 
		int N=20;
		for(int i=0;i<N;i++)
		{
			P0.setPosition(Vec3D(0.5*i,0.0,0.0));
			P0.setVelocity(Vec3D(0,-0.5*i, 0.0));
			P0.setRadius(0.1);
			particleHandler.copyAndAddObject(P0);
		}*/
		
		///Test case 3, Periodic collisions
		Particle P0;
		P0.setRadius(0.1);

		double alpha=constants::pi*6/5;
		double R=0.8;
		double V=1.0;
		P0.setPosition(Vec3D(cos(alpha)*R,sin(alpha)*R, 0.0));
		P0.setVelocity(Vec3D(cos(alpha)*V,sin(alpha)*V, 0.0));
		particleHandler.copyAndAddObject(P0);		

		alpha-=constants::pi;
		P0.setPosition(Vec3D(cos(alpha)*R,sin(alpha)*R, 0.0));
		P0.setVelocity(Vec3D(cos(alpha)*V,sin(alpha)*V, 0.0));
		particleHandler.copyAndAddObject(P0);		

		R=1.2;
		V=-1.0;
		P0.setPosition(Vec3D(cos(alpha)*R,sin(alpha)*R, 0.0));
		P0.setVelocity(Vec3D(cos(alpha)*V,sin(alpha)*V, 0.0));
		particleHandler.copyAndAddObject(P0);		


	}
	
	protected:
	
	///Couts time
	virtual void printTime() const {
		cout << "t=" << setprecision(3) <<fixed<< left << setw(6) << getTime() 
			<< ", tmax=" << setprecision(3) << left << setw(6) << getTimeMax() << endl;
		cout.flush();
	}
	
	/*void computeExternalForces(BaseParticle* CI)
	{
		/// Now add on gravity
		CI->addForce(getGravity() * CI->getMass());
		///Finally walls
		if (!CI->isFixed()) computeWalls(CI);
		CI->addForce(-CI->getPosition()*CI->getMass());
	}*/
	
	
};

int main(int /*argc*/, char **/*argv[]*/)
{
	///Start off by solving the default problem
	circular_periodic_walls_and_tangential_springs problem;
	problem.setName("circular_periodic_walls_and_tangential_springs");
	problem.setDensity(20);
	problem.species->setStiffness(1000);
	problem.species->setSlidingStiffness(10);
	problem.species->setSlidingFrictionCoefficient(1.0);
	problem.setGravity(Vec3D(0.0,0.0,0.0));
	problem.setTimeMax(0.14);
	problem.setTimeStep();
	problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimestep(200,getTimeMax(),getTimeStep()));
	problem.write(std::cout,false);
	problem.solve();	
}
