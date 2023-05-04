//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "Chute.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"

// New class that is the same as chute, but can have more walls
class Chute_tessa : public Chute{
public:

    Chute_tessa()
    {
        species  = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
    }
// this method overwrites the setupInitialConditions method in Chute
void setupInitialConditions() override
{
	//set side walls - solid if not a periodic
	if (getIsPeriodic())
		{
		  PeriodicBoundary b0;//set_NWallPeriodic(1);
		  wallHandler.clear();//set_NWall(0);
	          b0.set(Vec3D( 0.0, 1.0, 0.0), getYMin(), getYMax());
                  boundaryHandler.copyAndAddObject(b0);
		}
	else
		{
		  boundaryHandler.clear();//set_NWallPeriodic(0);
		  InfiniteWall w0;//set_NWall(3);
               	  w0.set(Vec3D( 0.0,-1.0, 0.0), getMin());
                  wallHandler.copyAndAddObject(w0);
		  w0.set(Vec3D( 0.0, 1.0, 0.0),  getMax());
                  wallHandler.copyAndAddObject(w0);
		
		// ------------------------------
		// extra wall:
		// muur op xmin, en de create_inflow_particle functie houdt er al rekening mee dat de deeltjes niet in de muur geplaatst worden
		  w0.set(Vec3D(-1.0, 0.0, 0.0), getMin());
                  wallHandler.copyAndAddObject(w0);
		// ------------------------------
		}

	//creates the bottom of the chute
	createBottom();
}
    LinearViscoelasticSlidingFrictionSpecies* species;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	Chute_tessa problem;
	
	// Problem parameters
	problem.setName("chute_demo_111025");
	problem.setTimeMax(10);
	//problem.set_HGRID_max_levels(4); // 2^(... -1) moet minimaal >= FixedParticleRadius/InflowParticleRadius (of andersom)
	
	// Particle properties
	problem.setFixedParticleRadius(0.0015);
	problem.setInflowParticleRadius(0.001,0.0015);
	//problem.species->setCollisionTimeAndRestitutionCoefficient(2.5e-3, 0.8);
	problem.species->setCollisionTimeAndRestitutionCoefficient(5e-4, 0.5,problem.species->getMassFromRadius(problem.getInflowParticleRadius()));
 	problem.species->setSlidingDissipation(0.0);
	
	// Chute properties
	//problem.setChuteAngle(30.0);
	problem.setChuteAngle(5);
	problem.setXMax(0.5);
	//problem.setZMax(5*problem.getInflowParticleRadius()); // dit doet niks
	problem.setYMax(2*problem.getInflowParticleRadius()); // de "dikte"
	
	
	// Inflow values
	problem.setInflowHeight(0.02);
	problem.setInflowVelocity(0.02);
	problem.setInflowVelocityVariance(0.02);
	
	//solve
	//std::cout << "Maximum allowed speed of particles: " << problem.getMaximumVelocity() << std::endl; // speed allowed before particles move through each other!
	problem.setTimeStep(5e-4/50.0); 
 	problem.setSaveCount(500);
	std::cout << "dt=" << problem.getTimeStep() << std::endl;
	problem.readArguments(argc, argv);
	problem.solve();
	problem.writeRestartFile();
}
