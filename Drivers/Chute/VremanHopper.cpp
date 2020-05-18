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
#include "scr/ChuteWithHopper.h"
#include "scr/Time.h"
//#include "scr/fstatistics.h"
using namespace std;

class Vreman : public ChuteWithHopper
{
public:
	void setupInitialConditions()
	{
		ChuteWithHopper::setupInitialConditions();
		set_symmetric_contraction(.3,.5,.03);
		write(cout);
		cout<< "tc=" << getCollisionTime() 
			<< ", eps="	<< getRestitutionCoefficient()
			<< ", vmax=" << getMaximumVelocity()
			<< ", InflowHeight/getZMax()=" << getInflowHeight()/getZMax()
			<< endl << endl;
	}
	
	void set_symmetric_contraction(double x_min, double x_max, double delta_y) {
		Walls.resize(Walls.size()+1);
		//back wall
		Vec3D normalIntoWall = Vec3D(-1,0,0);
		Vec3D point = Vec3D(x_max,0,0);
		Walls.back().addObject(normalIntoWall, Vec3D::dot(normalIntoWall,point));
		//slanted wall
		double delta_x = x_max-x_min;
		normalIntoWall = Vec3D(delta_y,-delta_x,0)/sqrt(mathsFunc::square(delta_x)+mathsFunc::square(delta_y));
		point = Vec3D(x_min,0,0);
		Walls.back().addObject(normalIntoWall, Vec3D::dot(normalIntoWall,point));

		Walls.resize(Walls.size()+1);
		//back wall
		normalIntoWall = Vec3D(-1,0,0);
		point = Vec3D(x_max,getChuteWidth(),0);
		Walls.back().addObject(normalIntoWall, Vec3D::dot(normalIntoWall,point));
		//slanted wall
		delta_x = x_max-x_min;
		normalIntoWall = Vec3D(delta_y,delta_x,0)/sqrt(mathsFunc::square(delta_x)+mathsFunc::square(delta_y));
		point = Vec3D(x_min,getChuteWidth(),0);
		Walls.back().addObject(normalIntoWall, Vec3D::dot(normalIntoWall,point));
	}
	
	void actionsAfterTimeStep() {
		cout << "t=" << setprecision(3) << left << setw(6) << getTime() << "N=" << setprecision(3) << left << setw(6) << get_N() << endl;
	}
	
	void printTime() const {
		//cout << "t=" << setprecision(3) << left << setw(6) << getTime() << "N=" << setprecision(3) << left << setw(6) << get_N() << endl;
	}

};

int main(int argc, char *argv[])
{
	Vreman problem;
	problem.setName("Vreman");
	problem.setChuteLength(.7);
	problem.setChuteWidth(.13);
	problem.setChuteAngle(19);
	problem.setDensity(2470);
	problem.setInflowParticleRadius(1e-3/2.);
	
	//Hopper properties
	double ExitHeight = 15e-3, ExitLength = 15e-3, hopperAngle_ = 45.0, hopperLength_ = 2.0 * ExitLength;
	problem.set_Hopper(ExitLength,ExitHeight,hopperAngle_,hopperLength_);
	problem.setMaxFailed(100);
	
	problem.setZMax(10e-3);
	problem.setFixedParticleRadius(0);
	problem.setStiffnessAndRestitutionCoefficient(100, 0.97, problem.get_mass_from_Radius(problem.getInflowParticleRadius()));
	problem.set_HGRID_num_buckets_to_power(4e5);
	problem.set_HGRID_max_levels(1);
	problem.setXBallsAdditionalArguments("-v0 -solidf -sort");
	problem.setTimeStep(problem.getCollisionTime()/50);
	problem.setTimeMax(20);
	problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(2000,problem.getTimeMax(),problem.getTimeStep()));
	problem.readArguments(argc, argv);
	problem.solve();
}
