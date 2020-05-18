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

#include "scr/ChuteWithHopper.h"

#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

class SegregationWithHopper : public ChuteWithHopper {
public:
	
	void create_inflow_particle()
	{
		
		int N=num_created;


		//This is 33% chance the particle is large
		int temp=(N%3)/2;

		//Difference in radius is sqrt 2 so factor of 2 in size.
		P0.Radius = 0.3e-3*(1.0+(sqrt(2)-1.0)*temp);		
		
		P0.computeMass(Species);
		
		static double s = sin(getChuteAngle());
		static double c = cos(getChuteAngle());
		static double Ht = tan(hopperAngle_);
		static double Hc = cos(hopperAngle_);
		static Vec3D AB = Vec3D(c,0.0,s);
		static Vec3D AC = Vec3D(-s,0.0,c);
		static Vec3D A = Vec3D(0.0, 0.0, s*(-0.5*(hopperLength_-hopperExitLength_)) + c*hopperHeight_) + AB*0.5*hopperLength_ + AC*(-0.5*hopperLength_/Ht);
		double gamma = random(0.5,1.0);
		P0.Position = A + AB * (random(-1.0,1.0)*(0.5*gamma*hopperLength_ - P0.Radius/Hc)) + AC * (gamma*0.5*hopperLength_/Ht);
		P0.Position.Y = random(getYMin()+P0.Radius, getYMax()-P0.Radius);
		
		
		
	}
	

};
int main(int argc UNUSED, char *argv[] UNUSED)
{
	SegregationWithHopper problem;
	
	// Problem parameters
	problem.setName("segregation");
	//Should be 10 for full length problem, but here I keep it low for a test case
	problem.setTimeMax(10);
	
	
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
	problem.setChuteLength(600.0e-3);
	problem.setChuteWidth(3e-3);
	problem.setMaxFailed(6);
	problem.makeChutePeriodic();
	double ExitHeight = 12.0e-3, ExitLength = 1.0 * ExitHeight, hopperAngle_ = 60.0, hopperLength_ = 6.0 * ExitLength;
	problem.set_Hopper(ExitLength,ExitHeight,hopperAngle_,hopperLength_);
	
	//solve
	cout << "Maximum allowed speed of particles: " << problem.getMaximumVelocity() << endl; // speed allowed before particles move through each other!
	problem.setTimeStep(); 
	problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(2000,problem.getTimeMax(),problem.getTimeStep()));
	cout << "dt=" << problem.getTimeStep() << endl;
	
	problem.auto_number();
	problem.write(std::cout,false);
	problem.save_info_to_disk();
	
	//This set to colouring based of size and small vectors
	problem.setXBallsColourMode(7);
	problem.setXBallsVectorScale(1);
	problem.setXBallsScale(20.0);
	
	
	problem.solve();
	problem.write(cout);
	//problem.HGRID_base::write(cout);
	
	
	//cout << problem << endl;
	problem.writeRestartFile();
}
