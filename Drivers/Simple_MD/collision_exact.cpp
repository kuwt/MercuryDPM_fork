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

#include "DPMBase.h"
#include <sstream>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>

class collision : public DPMBase {
	void computeExternalForces(int CI)
	{
		//Only walls, no gravity
		computeWalls(CI);
	}
	
	void setupInitialConditions()
	{
		// Two particles that immediately collide
		particleHandler.getObject(0)->Position=Vec3D(0.5*getXMax()-0.0005,0.005,0.0);
		particleHandler.getObject(1)->Position=Vec3D(0.5*getXMax()+0.0005,0.005,0.0);
		
		particleHandler.getObject(0)->setVelocity(Vec3D(0.1,0.0,0.0));
		particleHandler.getObject(1)->setVelocity(Vec3D(-0.1,0.0,0.0));
		
		particleHandler.getObject(0)->getRadius()=0.0005;
		particleHandler.getObject(1)->getRadius()=0.0005;
	}
	
};

void collision_numerical(double tmax, double dt, double dissipation, double k, int savecount);

void collision_exact(double tmax, double dt, double dissipation, double k, int savecount);

int main(int argc UNUSED, char *argv[] UNUSED)
{
	//set common parameters
	double tmax = 1.7e-5;
	double dt = 1e-7;
	double dissipation = 0.01;
	double k = 1e4;
	int savecount = 1;
	
	// get exact and numerical solution
	collision_numerical(tmax, dt, dissipation, k, savecount);
	collision_exact(tmax, dt, dissipation, k, savecount);
}

void collision_numerical(double tmax, double dt, double dissipation, double k, int savecount)
{
	//solving it numerically
	collision problem;
	problem.set_name("Collision");
	problem.setParticleDimensions(3);
	problem.setTimeMax(tmax);
	problem.setTimeStep(dt);
	problem.set_dissipation(dissipation);
	problem.speciesHandler.getObject(0)->setStiffness(k);
	problem.setSaveCount(savecount);
	problem.solve();
	//std::cout << problem;
	//std::cout << "tc " << problem.getCollisionTime_for_smallest_particle() << std::endl;
}

void collision_exact(double tmax, double dt, double dissipation, double k, int savecount)
{
	fstream data_file;
	data_file.open( "Collision_exact.data" , fstream::out);
	
	double t=0;
	dt *= savecount;
	double xmin = 0.0, xmax = 0.01, ymin = 0.0, ymax = 0.01;	
	double PositionX = 0.5*xmax-0.0005, PositionY = .005;//, PositionZ= 0.0;
	double VelocityX = 0.1, VelocityY = 0.0;//, VelocityZ= 0.0;
	double Radius = 0.0005;
	
	double rho = 2000;
	double mass = constants::pi/1.5 * Radius * Radius * Radius * rho;
	std::cout << "mass " << mass << std::endl;
	double mu = dissipation/mass;
	double mu_wall = dissipation/(2*mass);
	double w = sqrt(k/(mass/2) - mu*mu);
	double w_wall = sqrt(k/mass - mu*mu);
	
	double T0 = 0, T1;
	double Position0 = PositionX;
	double Velocity0 = VelocityX;
	while (t<tmax)
	{
		// move freely to the right
		T1 = min(tmax, max(T0, T0 + (xmax/2 - Radius - Position0) / Velocity0));
		std::cout << T1 << std::endl;
		while(t<T1) 
		{
			PositionX = Position0 + Velocity0 * (t-T0);
			VelocityX =             Velocity0;
			data_file << 2 << " " << t << " " <<xmin<<" "<<ymin<<" "<<xmax<<" "<<ymax<< std::endl;
			data_file << PositionX                 << " " << PositionY <<" " <<  VelocityX <<" " << VelocityY<< " " << Radius <<" -0 -0 0"<<std::endl;
			data_file << xmax - (PositionX - xmin) << " " << PositionY <<" " << -VelocityX <<" " << VelocityY<< " " << Radius <<" -0 -0 0"<<std::endl;
			t += dt;
		}
		Position0 = Position0 + Velocity0 * (T1-T0);
		T0 = T1;
		// collision with particle
		T1 = min(tmax, max(T0, T0 + constants::pi/w));
		std::cout << T1 << std::endl;
		while(t<T1) 
		{
			PositionX = Position0 + Velocity0/w * exp(-mu*(t-T0)) * sin(w*(t-T0));
			VelocityX =             Velocity0/w * exp(-mu*(t-T0)) * (-mu * sin(w*(t-T0)) + w * cos(w*(t-T0))) ;
			data_file << 2 << " " <<t << " " <<xmin<<" "<<ymin<<" "<<xmax<<" "<<ymax<< std::endl;
			data_file << PositionX                 << " " << PositionY <<" " <<  VelocityX <<" " << VelocityY<< " " << Radius <<" -0 -0 0"<<std::endl;
			data_file << xmax - (PositionX - xmin) << " " << PositionY <<" " << -VelocityX <<" " << VelocityY<< " " << Radius <<" -0 -0 0"<<std::endl;
			t += dt;
		}
		Position0 = Position0 + Velocity0/w * exp(-mu*(T1-T0)) * sin(w*(T1-T0));
		Velocity0 =             Velocity0/w * exp(-mu*(T1-T0)) * (-mu * sin(w*(T1-T0)) + w * cos(w*(T1-T0))) ;
		T0 = T1;
		// move freely to the left
		T1 = min(tmax, max(T0, T0 + (xmin + Radius - Position0) / Velocity0));
		std::cout << T1 << std::endl;
		while(t<T1) 
		{
			PositionX = Position0 + Velocity0 * (t-T0);
			VelocityX =             Velocity0;
			data_file << 2 << " " <<t << " " <<xmin<<" "<<ymin<<" "<<xmax<<" "<<ymax<< std::endl;
			data_file << PositionX                 << " " << PositionY <<" " <<  VelocityX <<" " << VelocityY<< " " << Radius <<" -0 -0 0"<<std::endl;
			data_file << xmax - (PositionX - xmin) << " " << PositionY <<" " << -VelocityX <<" " << VelocityY<< " " << Radius <<" -0 -0 0"<<std::endl;
			t += dt;
		}
		Position0 = Position0 + Velocity0 * (T1-T0);
		T0 = T1;
		// collision with wall
		T1 = min(tmax, max(T0, T0 + constants::pi/w_wall));
		std::cout << T1 << std::endl;
		while(t<T1) 
		{
			PositionX = Position0 + Velocity0/w_wall * exp(-mu_wall*(t-T0)) * sin(w_wall*(t-T0));
			VelocityX =             Velocity0/w_wall * exp(-mu_wall*(t-T0)) * (-mu_wall * sin(w_wall*(t-T0)) + w_wall * cos(w_wall*(t-T0))) ;
			data_file << 2 << " " <<t << " " <<xmin<<" "<<ymin<<" "<<xmax<<" "<<ymax<< std::endl;
			data_file << PositionX                 << " " << PositionY <<" " <<  VelocityX <<" " << VelocityY<< " " << Radius <<" -0 -0 0"<<std::endl;
			data_file << xmax - (PositionX - xmin) << " " << PositionY <<" " << -VelocityX <<" " << VelocityY<< " " << Radius <<" -0 -0 0"<<std::endl;
			t += dt;
		}
		Position0 = Position0 + Velocity0/w_wall * exp(-mu_wall*(T1-T0)) * sin(w_wall*(T1-T0));
		Velocity0 =             Velocity0/w_wall * exp(-mu_wall*(T1-T0)) * (-mu_wall * sin(w_wall*(T1-T0)) + w_wall * cos(w_wall*(T1-T0))) ;
		T0 = T1;
	} 
	data_file.close();
	
	fstream script_file;
	script_file.open("Collision_exact.disp" , fstream::out);
	double scale = 1 / max( ymax-ymin, xmax-xmin );
	script_file << "../xballs -format 8 -f " << "Collision_exact.data"
		<<" -s "<<scale<<" -cmode 0 -cmax -scala 4 $*";
	script_file.close();
	chmod("Collision_exact.disp",S_IRWXU);
}
