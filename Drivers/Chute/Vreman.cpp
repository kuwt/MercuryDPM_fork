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
#include "scr/Time.h"
//#include "scr/fstatistics.h"
using namespace std;

class Vreman : public Chute
{
public:
	Vreman() {
		InflowNu=0;
	}

	void setupInitialConditions()
	{
		Chute::setupInitialConditions();
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
	
	void set_InflowNu(double Nu) {InflowNu=Nu; cout << "Nu=" <<  Nu << endl;};
	double InflowNu;

	void add_particles() 
	{
		static double Radius = 1.04*(1.+2.*InflowVelocityVariance)*MaxInflowParticleRadius;
		static double frequency = InflowVelocity/(4./3.*sqrt(6)*Radius);
		static int N = round(InflowNu*InflowVelocity/frequency*getChuteWidth()*InflowHeight/(8.*pi/6.*mathsFunc::cubic(getInflowParticleRadius())));
		static double fy = 2.*Radius;
		static double fz = 2.*sqrt(3)*Radius;
		static double fx = 2./3.*sqrt(6)*Radius;
		if (ceil(getTime()*frequency)!=ceil((getTime()-getTimeStep())*frequency)) {
			P0.Radius = random(MinInflowParticleRadius,MaxInflowParticleRadius);
			P0.computeMass(Species);
			int Ny = (getYMax()-Radius-(getYMin() + 2.*Radius))/fy;
			int Nz = ceil(InflowHeight/fz);
			if (N>4*Ny*Nz) {cout << "InflowNu too big: " << InflowNu << ", could only place" << (4.*Ny*Nz)/N*100. << "%" << endl; exit(-1); }
			for (int j=0;j<Nz;j++)
			for (int i=0;i<Ny;i++) 
			{
				P0.Position.Y = getYMin() + 2.*Radius + fy*i;
				P0.Position.Z = getZMin() + Radius + fz*j;
				P0.Position.X = getXMin() + Radius;
				P0.Velocity.X = InflowVelocity * random(-InflowVelocityVariance,InflowVelocityVariance) + InflowVelocity;
				P0.Velocity.Y = 0;
				P0.Velocity.Z = 0;
				if (random(0,4*Ny*Nz)<N) add_particle(P0);
				P0.Position.Y = getYMin() + Radius + fy*i;
				P0.Position.Z = getZMin() + (1+sqrt(3))*Radius + fz*j;
				P0.Position.X = getXMin() + Radius;
				P0.Velocity.X = InflowVelocity * random(-InflowVelocityVariance,InflowVelocityVariance) + InflowVelocity;
				P0.Velocity.Y = 0;
				P0.Velocity.Z = 0;
				if (random(0,4*Ny*Nz)<N) add_particle(P0);
				P0.Position.Y = getYMin() + 2.*Radius + fy*i;
				P0.Position.Z = getZMin() + 4./sqrt(3)*Radius + fz*j;
				P0.Position.X = getXMin() + Radius+fx;
				P0.Velocity.X = InflowVelocity * random(-InflowVelocityVariance,InflowVelocityVariance) + InflowVelocity;
				P0.Velocity.Y = 0;
				P0.Velocity.Z = 0;
				if (random(0,4*Ny*Nz)<N) add_particle(P0);
				P0.Position.Y = getYMin() + Radius + fy*i;
				P0.Position.Z = getZMin() + 1./sqrt(3)*Radius +fz + fz*j;
				P0.Position.X = getXMin() + Radius + fx;
				P0.Velocity.X = InflowVelocity * random(-InflowVelocityVariance,InflowVelocityVariance) + InflowVelocity;
				P0.Velocity.Y = 0;
				P0.Velocity.Z = 0;
				if (random(0,4*Ny*Nz)<N) add_particle(P0);
			}
			cout << getTime() <<" " << get_N() <<endl;
		}
	}
	
	void actionsAfterTimeStep() {
		static double InflowSize = 2.*4./3.*sqrt(6)*MaxInflowParticleRadius;
		for (vector<CParticle>::iterator it = Particles.begin(); it!=Particles.end(); ++it) {
			if (it->Position.X<getXMin()+InflowSize) {
				it->Force.set_zero();
				it->Torque.set_zero();
			}
		}
	}
	
	void printTime() const {
		//cout << "t=" << setprecision(3) << left << setw(6) 
		//<< getTime() << "N=" << setprecision(3) << left << setw(6) << get_N() << endl;
	}

};

int main(int argc, char *argv[])
{
	Vreman problem;
	problem.setName("Vreman");
	problem.setChuteLength(.7);
	problem.setChuteWidth(.13);
	problem.setChuteAngle(19);
	//problem.makeChutePeriodic();
	problem.setDensity(2470);
	problem.setInflowParticleRadius(1e-3/2);
	problem.setInflowVelocity(.17);
	problem.setInflowVelocityVariance(.02);
	problem.setInflowHeight(8e-3);
	problem.setZMax(10e-3);
	problem.set_InflowNu(0.29/(problem.getInflowHeight()*problem.getChuteWidth()*problem.getInflowVelocity()*problem.getDensity()));
	problem.setFixedParticleRadius(0);
	problem.setStiffnessAndRestitutionCoefficient(100, 0.97, problem.get_mass_from_Radius(problem.getInflowParticleRadius()));
	problem.set_HGRID_num_buckets_to_power(4e5);
	problem.set_HGRID_max_levels(1);
	problem.setXBallsAdditionalArguments("-v0 -solidf -sort");
	problem.setTimeStep(problem.getCollisionTime()/50);
	problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(2000,problem.getTimeMax(),problem.getTimeStep()));
	problem.setTimeMax(20);
	problem.readArguments(argc, argv);
	problem.solve();
}
