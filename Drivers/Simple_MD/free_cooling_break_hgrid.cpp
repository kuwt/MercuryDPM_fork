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

#include "DPMBase.h"
#include "scr/Mercury2D.h"


//#define DEBUG_OUTPUT

class my_problem : public DPMBase{
public:
	
	void setupInitialConditions()
	{
	
		int N=get_N();
		int N1=static_cast<int>(sqrt(N))+1;
		
		for (int i=0;i<N;i++)
		{
		
			int ix=static_cast<int>(i%N1);
			int iy=static_cast<int>(i/N1);
		
			double x=(getXMax()-getXMin())*(ix+1)/(N1+1);
			double y=(getYMax()-getYMin())*(iy+1)/(N1+1);
			
		
			particleHandler.getObject(i)->Position=Vec3D(x,y,0.0);
			particleHandler.getObject(i)->setVelocity(Vec3D(0.1,0.1,0.0));
			particleHandler.getObject(i)->getRadius()=0.0001;
	
		}	
	
	}

private:
	void computeExternalForces(int P){computeWalls(P);}
		
};

class my_problem_HGRID : public Mercury2D{
public:

	my_problem_HGRID(my_problem& other) : Mercury2D(other){}
	
	void setupInitialConditions()
	{
	
		int N=get_N();
		int N1=static_cast<int>(sqrt(N))+1;
		
		for (int i=0;i<N;i++)
		{
		
			int ix=static_cast<int>(i%N1);
			int iy=static_cast<int>(i/N1);

			
		
			double x=(getXMax()-getXMin())*(ix+1)/(N1+1);
			double y=(getYMax()-getYMin())*(iy+1)/(N1+1);
			
		
			particleHandler.getObject(i)->Position=Vec3D(x,y,0.0);
			particleHandler.getObject(i)->setVelocity(Vec3D(0.1,0.1,0.0));
			particleHandler.getObject(i)->getRadius()=0.0001;
	
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
	//Note should be 1000 particles

	

	///Start off my solving the default problem
 	my_problem problem;
 	problem.set_N(1024);
 	problem.set_name("free_cooling");
 	problem.set_dissipation(0.005);
	problem.speciesHandler.getObject(0)->setStiffness(1e3);
	problem.setTimeStep(5e-5);
	problem.setSaveCount(10);
 	problem.setTimeMax(0.5);
 	std::cout << "Solving by checking all collisions " << std::endl;
 	problem.solve();

	///Now do the same problem with the HGrid turned on, simple done by create a new problem from the old problem.
	my_problem_HGRID problem2(problem);
	problem2.set_name("free_cooling_hgrid");
	//problem2.set_HGRID_min_cell_size(0.0001);
	std::cout << "... now using the Hgrid" << std::endl;
	//The idea is with only one bucket the Hgrid is also testing all collision so now the different is really only in the order. This also some the small driff between the two solutions.
	problem2.set_HGRID_num_buckets(1);
	problem2.solve();
		
/*	std::cout << std::endl << "v=" << problem2.getMaximumVelocity_of_smallest_particle() << std::endl;
	std::cout << "tc=" << problem2.getCollisionTime_for_smallest_particle() << std::endl;
	std::cout << "eps=" << problem2.getRestitutionCoefficient_for_smallest_particle() << std::endl;*/
	
}
