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
#include <sys/types.h>
#include <sys/stat.h>
#include "scr/Chute.h"
//#include "scr/fstatistics.h"
using namespace std;

class ChutePeriodic : public Chute{
public:

	void actionsBeforeTimeStep(){
		if  (getTime()>1&&getTime()-getTimeStep()<1) {
			Walls[Walls.size()-1].set(Vec3D(1,0,0),getXMax());
		}
	};
		
	void create_inflow_particle()
	{
		P0.Radius = random(MinInflowParticleRadius,MaxInflowParticleRadius);
		P0.computeMass(Species);
		
		P0.Position.X = random(getXMin()+P0.Radius,get_hopperLength_()-P0.Radius);
		P0.Position.Y = random(getYMin()+P0.Radius,getYMax()-P0.Radius);
		P0.Position.Z = random(getZMin()+2.0*P0.Radius,getZMax());
		//the 1.8 is there because the initial packing compacts
		P0.Velocity = Vec3D(0.0,0.0,0.0);
	}
		
	void setupInitialConditions()
	{
		
		Chute::setupInitialConditions();
		
		//~ set_Nmax(10000);
		HGridActionsBeforeTimeLoop();
		HGridActionsBeforeTimeStep();

		add_particles();
		
		int N=Walls.size();
		Walls.resize(N+2);
		Walls[N].set(Vec3D(-1,0,0),0);
		Walls[N+1].set(Vec3D(1,0,0),get_hopperLength_());
		fix_hgrid();
		write(std::cout,false);
	}
	
	void set_hopperLength_(double new_) {hopperLength_=new_;}
	double get_hopperLength_() {return hopperLength_;}

private:
	double hopperLength_;

};


int main(int argc, char *argv[])
{
	ChutePeriodic md;
	
	// Problem parameters
	md.setName("Lawine");
	md.fStatFile.setFileType(FileType::ONE_FILE);
	md.dataFile.setFileType(FileType::ONE_FILE);
	md.autoNumber();
	md.setTimeMax(100);
 
	// Particle properties
	md.setDensity(2400.0);
	md.setInflowParticleRadius(2.5e-3);//d=5mm
	md.speciesHandler.getObject(0)->setCollisionTimeAndRestitutionCoefficient(5e-4, 0.8);
	md.setTimeStep(1e-5); 
	md.setTimeMax(25); 
	md.speciesHandler.getObject(0)->setSlidingDissipation(md.get_dissipation());
	md.speciesHandler.getObject(0)->setSlidingStiffness(2./7.*md.speciesHandler.getObject(0)->getStiffness());
	md.speciesHandler.getObject(0)->setSlidingFrictionCoefficient(0.5);
	cout << "friction angle " << atan(md.getSlidingFrictionCoefficient())*180./pi << endl;
	md.setFixedParticleRadius(0);
	md.setFixedParticleRadius(md.getInflowParticleRadius());
	md.setRoughBottomType(MONOLAYER_ORDERED);
	
	// Chute properties
	md.setChuteAngle(0.0);
	md.setChuteWidth(450e-3);
	md.setZMax(150e-3);//initial height 150-200mm
	md.setChuteLength(1);
	md.setMaxFailed(1000);
	
	if (argc>1&&argv[1][0]!='-') {
		double ratio = pow(atof(argv[1]),1./3.);
		md.set_hopperLength_(md.getZMax()/mathsFunc::square(ratio));
		md.setChuteWidth(md.getChuteWidth()*ratio);
		md.setZMax(md.getZMax()*ratio);
		argc--;	argv++;
	} else md.set_hopperLength_(md.getZMax());
		
	
	cout << "hopperLength_=" << md.get_hopperLength_() << endl; 
	cout << "ChuteWidth=" << md.getChuteWidth() << endl; 
	cout << "zmax=" << md.getZMax() << endl; 
	cout << "Maximum allowed speed of particles: " << md.getMaximumVelocity() << endl; 
	md.setSaveCount(.5e4);
	cout << "dt=" << md.getTimeStep() << endl;
	
	//This set to colouring based of size and small vectors
	md.setXBallsColourMode(7);
	md.setXBallsVectorScale(1);
		
	//solve
	md.readArguments(argc, argv);
	md.solve();
	
	md.write(std::cout,false);
	md.writeRestartFile();
	md.writeRestartFile();
	
}
