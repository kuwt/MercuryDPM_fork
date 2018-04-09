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

void actionsBeforeTimeStep(){};
	
void create_inflow_particle()
{
	P0.Radius = random(MinInflowParticleRadius,MaxInflowParticleRadius);
	P0.computeMass(Species);
	
	P0.Position.X = random(getXMin()+2.0*P0.Radius,getXMax());
	P0.Position.Y = random(getYMin()+2.0*P0.Radius,getYMax());
	P0.Position.Z = random(getZMin()+2.0*P0.Radius,1.8*getZMax());
	P0.Velocity = Vec3D(0.0,0.0,0.0);
}
	
void setupInitialConditions()
{
	
	Chute::setupInitialConditions();
	
	set_NWallPeriodic(2);
	WallsPeriodic[0].set(Vec3D( 1.0, 0.0, 0.0), getXMin(), xmax);
	WallsPeriodic[1].set(Vec3D( 0.0, 1.0, 0.0), getYMin(), getYMax());
	
	set_Nmax(10000);
	HGridActionsBeforeTimeLoop();
	add_particles();
}

void setup() {
	// Problem parameters
	setName("chute_periodic");
	fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
	dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
	autoNumber();
	setTimeMax(10);
 
	// Particle properties
	setDensity(2400.0);
	setInflowParticleRadius(.5e-3);
	setCollisionTimeAndRestitutionCoefficient(4e-4, 0.8);
	setSlidingDissipation(get_dissipation());
	setSlidingStiffness(getStiffness());
	setSlidingFrictionCoefficient(0.5);
	setFixedParticleRadius(getInflowParticleRadius());
	setRoughBottomType(MULTILAYER);
	
	// Chute properties
	setChuteAngle(25.0);
	setChuteLength(20e-3);
	setChuteWidth(10e-3);
	setZMax(20e-3);
	setMaxFailed(1000);
	makeChutePeriodic();
	
	cout << "Maximum allowed speed of particles: " << getMaximumVelocity() << endl; 
	setTimeStep(); 
	setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimestep(150,getTimeMax(),getTimeStep()));
	cout << "dt=" << getTimeStep() << endl;
	
	//This set to colouring based of size and small vectors
	setXBallsColourMode(7);
	setXBallsVectorScale(1);
}

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	ChutePeriodic problem;
	problem.setup();
	if (argc>1) problem.setChuteAngle(atof(argv[1]));
	else {cerr << "Argument needed!" << endl; exit(-1);}
	cout << "Chute Angle: " << problem.getChuteAngle() << endl;
	problem.setZMax(5e-3);
	problem.solve();
	problem.write(std::cout,false);
	problem.save_info_to_disk();
	problem.writeRestartFile();
	
	ChutePeriodic problem2;
	problem2.setup();
	if (argc>1) problem2.setChuteAngle(atof(argv[1]));
	else {cerr << "Argument needed!" << endl; exit(-1);}
	cout << "Chute Angle: " << problem2.getChuteAngle() << endl;
	problem2.setZMax(10e-3);
	problem2.solve();
	problem2.write(std::cout,false);
	problem2.save_info_to_disk();
	problem2.writeRestartFile();

	ChutePeriodic problem3;
	problem3.setup();
	if (argc>1) problem3.setChuteAngle(atof(argv[1]));
	else {cerr << "Argument needed!" << endl; exit(-1);}
	cout << "Chute Angle: " << problem3.getChuteAngle() << endl;
	problem3.setZMax(15e-3);
	problem3.solve();
	problem3.write(std::cout,false);
	problem3.save_info_to_disk();
	problem3.writeRestartFile();
	
	ChutePeriodic problem4;
	problem4.setup();
	if (argc>1) problem4.setChuteAngle(atof(argv[1]));
	else {cerr << "Argument needed!" << endl; exit(-1);}
	cout << "Chute Angle: " << problem4.getChuteAngle() << endl;
	problem4.setZMax(20e-3);
	problem4.solve();
	problem4.write(std::cout,false);
	problem4.save_info_to_disk();
	problem4.writeRestartFile();
	
	ChutePeriodic problem5;
	problem5.setup();
	if (argc>1) problem5.setChuteAngle(atof(argv[1]));
	else {cerr << "Argument needed!" << endl; exit(-1);}
	cout << "Chute Angle: " << problem5.getChuteAngle() << endl;
	problem5.setZMax(25e-3);
	problem5.solve();
	problem5.write(std::cout,false);
	problem5.save_info_to_disk();
	problem5.writeRestartFile();
}
