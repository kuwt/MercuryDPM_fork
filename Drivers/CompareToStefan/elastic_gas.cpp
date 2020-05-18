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

#include<iostream>
#include <Species/LinearViscoelasticFrictionSpecies.h>

#include "Mercury2D.h"
#include "StatisticsVector.h"

//#define DEBUG_OUTPUT

class my_problem: public Mercury2D{
public:

	void setupInitialConditions()
	{
		
		if (readDataFile("c3d.ini",7))
			{
			  PeriodicBoundary b0;
			  b0.set(Vec3D(1,0,0), getXMin(), getXMax());
              boundaryHandler.copyAndAddObject(b0);
              b0.set(Vec3D(0,1,0), getYMin(), getYMax());
			  boundaryHandler.copyAndAddObject(b0);
			  std::cout << " Problem setup, about to solve" << std::endl;
			
			}
			
			else
			
			{
			  std::cerr << "Input data not found exiting " << std::endl;
				exit(-1);
			}
		
	}

private:
	void computeExternalForces(BaseParticle* P){}
		
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//* This code does the free cooling problem and compares it to the results of Stefans code all the particle data is loaded from an ini file
int main(int argc, char *argv[])
{

	

	///Start off my solving the default problem
 	my_problem problem;
 	problem.setName("free_cooling");
 	
 	//Stefan problem
 	//problem.set_dissipation(0.0);
//	species->setStiffness(1e5);
//	problem.setTimeStep(4e-5);
//	problem.setSaveCount(250000);
// 	problem.setTimeMax(10);

	//Vits problem
    auto species = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
	species->setDissipation(0.0);
	problem.setTimeStep(4e-5);
    species->setStiffness(1e6);
	species->setDensity(2000);
	problem.setTimeMax(12e-5);
	//problem.set_HGRID_max_levels(1);
	//problem.set_HGRID_num_buckets(1e5);
	problem.setParticleDimensions(2);
	problem.setGravity(Vec3D(0.0,0.0,0.0));
	problem.setSaveCount(1e5);
	
	

 	problem.write(std::cout,false);
	problem.solve();
	problem.write(std::cout,false);
	problem.writeRestartFile();
	
	//StatisticsVector problemstats(problem.getName(), 100, 100, 0, 0.05, Gaussian);

	StatisticsVector<XY> problemstats;
	problemstats.setName(problem.getName());
	problemstats.setNX(100);
        problemstats.setNY(100);
        problemstats.setNZ(0);
        problemstats.setCGWidth(0.05);
        problemstats.setCGShape(Gaussian);
	// 0, 0.05, Gaussian);
	problemstats.statistics_from_fstat_and_data();

	
		
/*	cout << endl << "v=" << problem2.getMaximumVelocity_of_smallest_particle() << endl;
	cout << "tc=" << problem2.getCollisionTime_for_smallest_particle() << endl;
	cout << "eps=" << problem2.getRestitutionCoefficient_for_smallest_particle() << endl;*/
}
