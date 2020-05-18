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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "scr/Chute.h"

#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

///This class does segregation problems in a periodic chute 
class SegregationPeriodic : public Chute{
public:


//This code requires you do not nothing special after each time step
void actionsBeforeTimeStep(){};

///This is the info call
void write(  	std::ostream &   	 os, bool  	print_all = false) 	
{
	os << "This is a segregation chute code problem " << endl;
	os << "\n \n \n"<< endl;
	
		
	 MercuryBase::write(os, print_all); 
	 
	 os << "Large particle size : " << radius_l << endl;
	 os << "Small particle size : " << radius_s << endl;

	
	
}


/// This setup the intial conditions, generates small volume fraction of particles. 
/// Sets the program to be periodic in x.
/// \bug This code is not non-dimensionalised at the moment, should do this shortly, but at the moment
void setupInitialConditions()
{
	
	//Check if the run has been done before. If yes, skip and start next run
		if (FileExists(data_filename.str()))
		{
			//If it has move on to teh next run immedently
			cout << "This run has been done " << endl;
			launchNewRun("segregation",true);
			exit(0);
		}
		
		
	//Set up a 10 by 10 study
	vector<int> study_num=get2DParametersFromRunNumber(10,1);
	
	
	//If study 0 is complete quit
		if (study_num[0] > 0)
			{
				cout << "Study is complete " << endl;
				exit(0);
			}
		else
		//If the study is not complete save the data to disk and move on
			{
				writeRestartFile();
				launchNewRun("segregation");
			}
			
	//Now setup the particles
	
	//Setup the base i.e. the chute particles
	Chute::setupInitialConditions();
	
	
	
	
	
	//Set up the walls
	set_NWallPeriodic(2);
	WallsPeriodic[1].set(Vec3D( 1.0, 0.0, 0.0), getXMin(), xmax);

	
	//Number of small particles
	int Ns=5000;
	//Small partilce radius
	radius_s=0.3e-3;
	//Radius of large particles, changes from study to study.
	radius_l=radius_s*(1.0+study_num[1]/10.0);
	//Number of large partices, fixed to the keep the volume fraction of large and small paritlces equal.
	int Nl=pow(radius_s/radius_l,3)*Ns;
	
	
	
	
	while ((Ns>0) && (Nl>0))
		{


		//random to see if want to generate a large or small particles, helps makes the initial conditions homogenious 
		if (random(1.0,Nl+Ns) > Nl)
			{
			//Generate a small particle: set radius to small radius subtract one off the list of small particles to be generated
			P0.Radius=radius_s;
			Ns--;
			}
		else
			{
			//Generate a large particle: set radius to large radius subtract one of the list of large particles to be generated
			P0.Radius=radius_l;
			Nl--;
			}
		
		
		//P0.Radius = 0.3e-3*(1.0+(sqrt(2)-1.0)*temp);	
		P0.Angle.set_zero();
		P0.AngularVelocity.set_zero();
		P0.computeMass(Species);
	
		//randomize particle position, zero intial velocity
		P0.Position.X = random(getXMin(),getXMax());
		P0.Position.Y = random(getYMin(),getYMax());
		P0.Position.Z = random(getZMin(),getZMax());
		P0.Velocity = Vec3D(0.0,0.0,0.0);
		
		
		//Add the new particle to the list of current particles
		Particles.push_back (P0); 

		
		}
		
		//Write the info to the screen and save a copy to the disk
		write(std::cout,false);
		writeRestartFile();
	
}

private:
double radius_s;
double radius_l;

};


int main(int argc UNUSED, char *argv[] UNUSED)
{
	SegregationPeriodic problem;
	
	// Problem parameters
	problem.setName("segregation");
	problem.setTimeMax(100);
	
	
	// Particle properties
	problem.setDensity(2400.0);
	
	problem.setInflowParticleRadius(0.3e-3,0.60e-3);
	problem.speciesHandler.getObject(0)->setCollisionTimeAndRestitutionCoefficient(4e-4, 0.6);
	problem.speciesHandler.getObject(0)->setSlidingDissipation(problem.get_dissipation());
	problem.speciesHandler.getObject(0)->setSlidingFrictionCoefficient(0.8);
	problem.setFixedParticleRadius(0.3e-3);
	problem.setRoughBottomType(MONOLAYER_DISORDERED);
	
	
	// Chute properties
	problem.setChuteAngle(25.0);
	problem.setChuteLength(50.0e-3);
	problem.setChuteWidth(3e-3);
	problem.setMaxFailed(6);
	problem.makeChutePeriodic();
	
	//solve
	cout << "Maximum allowed speed of particles: " << problem.getMaximumVelocity() << endl; // speed allowed before particles move through each other!
	problem.setTimeStep(); 
	//This is based on the fact in general you get too much data, so prob at worst you want to turn it into a 20 at 60fps (which is its self overkill)
	problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(20*60,getTimeMax(),getTimeStep()));
	//problem.setSaveCount(1);
	cout << "dt=" << problem.getTimeStep() << endl;
	
	problem.auto_number();

	
	//This set to colouring based of size and small vectors
	problem.setXBallsColourMode(7);
	problem.setXBallsVectorScale(1);
	problem.setXBallsAdditionalArguments("-sort -v0 -solidf");
	
	//problem.setSaveCount(1000);
	
	problem.solve();
	//problem.HGRID_base::write(cout);
	
	
	//cout << problem << endl;
	problem.writeRestartFile();
}





