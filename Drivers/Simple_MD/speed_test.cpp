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
#include "scr/Time.h"


//#define DEBUG_OUTPUT

class my_problem_HGRID : public Mercury3D
	{
	public:

//	my_problem_HGRID(my_problem& other) : Mercury3D(other){}
	
	void setupInitialConditions()
		{

		set_NWall(6);
		Walls[0].set(Vec3D(-1,0,0), -getXMin());
		Walls[1].set(Vec3D( 1,0,0),  getXMax());
		Walls[2].set(Vec3D(0,-1,0), -getYMin());
		Walls[3].set(Vec3D(0, 1,0),  getYMax());
		Walls[4].set(Vec3D(0,0,-1), -getZMin());
		Walls[5].set(Vec3D(0,0,1), getZMax());

	
	
		int N=get_N();
		int N1=static_cast<int>(pow(N,0.33))+1;
		
		
		for (int i=0;i<N;i++)
			{
		
			int ix=static_cast<int>(i%N1);
			int iz=static_cast<int>(i/N1/N1);
			
			int iy=(i-ix-N1*N1*iz)/N1;
		
			double x=(getXMax()-getXMin())*(ix+1)/(N1+1);
			double y=(getYMax()-getYMin())*(iy+1)/(N1+1);
			double z=(getZMax()-getZMin())*(iz+1)/(N1+1);
			
		
			particleHandler.getObject(i)->Position=Vec3D(x,y,z);
			particleHandler.getObject(i)->setVelocity(Vec3D(random(-1.0,1.0),random(-1.0,1.0),random(-1.0,1.0));
			particleHandler.getObject(i)->getRadius()=1.0/200.0;
	
			}	
			
			
	
		}
	
	private:
		void computeExternalForces(int P){computeWalls(P);}
		
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc UNUSED, char *argv[] UNUSED)
{
	
	

	///Start off my solving the default problem
 	my_problem_HGRID problem;
 	problem.auto_number();
 	
	int N=10;
	
	for (int i=1;i<problem.get_counter();i++)
		{
			N=N*10;
		}
 	
 	problem.set_N(N);
 	problem.set_name("speed_test");
 	problem.set_dissipation(0.005);
	problem.speciesHandler.getObject(0)->setStiffness(1e3);
	problem.setTimeStep(5e-5);
	problem.setSaveCount(1000);
 	problem.setTimeMax(1.0);
	problem.setSystemDimensions(3);
	problem.setDensity(2500);
	problem.setParticleDimensions(3);
	problem.setXMax(2);
	problem.setYMax(2);
	problem.setZMax(2);
	problem.save_info_to_disk();
	problem.set_HGRID_num_buckets(N);

	
	Time clock;
	
	clock.tic();
	
	problem.solve();
	problem.save_info_to_disk();
	
	std::cout << clock.toc();
	
	if (problem.get_counter() > 6)
		{
			exit(0);
		}
	
	problem.launch_new("speed_test",true);
		
	
}

