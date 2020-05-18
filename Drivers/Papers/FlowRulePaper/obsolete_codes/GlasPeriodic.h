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

#include "scr/Chute.h"
//~ #include "scr/Time.h"
using namespace std;

class SilbertPeriodic : public Chute
{
public:

	SilbertPeriodic() {
		// Problem parameters
		setName("silbert");
		
		//time stepping
		setTimeStep(1e-4);
		setTimeMax(2000);

		//output parameters
		setSaveCount(50e4);
	 
		//particle radii
		setInflowParticleRadius(.5);
		setFixedParticleRadius(.5);//getInflowParticleRadius());
		setRoughBottomType(MULTILAYER);

		//particle properties
		setDensity(6/pi);
		//~ setStiffnessAndRestitutionCoefficient(2e5,0.97,1);
		double tc=5e-3, r=0.97, beta=0.44, mu=0.092, mur=0.042;
		setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc,r,beta,1.0,1.0);
		setSlidingFrictionCoefficient(mu);
		setSlidingFrictionCoefficientr();
		
		//chute properties
		setChuteAngle(24.0, 1.0);
		setChuteLength(20);
		setChuteWidth(10);
		set_H(20);
			
	}
	
	void fix_hgrid() {
		//assume 1-2 levels are optimal (which is the case for mono and bidispersed) and set the cell size to min and max 
		// !this is not optimal for polydispersed
		double minCell = 2.*min(getFixedParticleRadius(),getMinInflowParticleRadius());
		double maxCell = 2.*max(getFixedParticleRadius(),getMaxInflowParticleRadius());
		if ((minCell==maxCell)|(minCell==0.)) set_HGRID_max_levels(1);
		else set_HGRID_max_levels(2);
		set_HGRID_cell_to_cell_ratio (1.0000001*maxCell/minCell);
	}

	double getSlidingFrictionCoefficientBottom() { 
		if (speciesHandler.getNumberOfObjects()>1) return speciesHandler.getMixedObject(1,0)->getSlidingFrictionCoefficient(); 
		else return getSlidingFrictionCoefficient(); 
	}
	
	void setSlidingFrictionCoefficientBottom(double new_) { 
		createBaseSpecies(); 
		speciesHandler.getMixedObject(1, 0)->setSlidingFrictionCoefficient(new_); 
	}
	
	void createBaseSpecies() {
		//only create once
		static bool created=false;
		if (!created) {
			speciesHandler.copyAndAddObject(speciesHandler.getObject(0));
			created = true;
			for (int i=0; i<get_N(); i++) {
				if (Particles[i].is_fixed()) Particles[i].indSpecies=1;
			}
		}
	}

	void set_study() {
		stringstream name;
		name << "H" << getInflowHeight() 
			 << "A" << getChuteAngleDegrees() 
			 << "L" << round(100.*getFixedParticleRadius()*2.)/100.
			 << "M" << getSlidingFrictionCoefficient()
			 << "B" << getSlidingFrictionCoefficientBottom();
		setName(name.str().c_str());
		set_data_filename();
	}

	void set_study(int study_num) {
		//S=0-5: lambda = 0, 3./6., 4./6., 5./6., 1, 2
		//S=6-8: mu = 0, 1, inf
		//S=9-13: mub = 0,1,inf,1/4,1/8
		//S=14-15: mu = 1/4, 1/8
		//S=16-19: lambda = 1./6., 2./6., 1.5, 4
		//S=21-25: mub=1/16,1/32,1/64,1/128,1/1024
		//S=26-28: lambda=1/2, mub=1/16,1/128,1/1024
		//S=29-32: lambda=0, mub=1/16,1/128,1/1024,0
		
		if (study_num < 6) {
			// set mu_all = 0.5, vary lambda
			double Lambdas[] = {0, 3./6., 4./6., 5./6., 1, 2};
			setFixedParticleRadius(Lambdas[study_num]/2.);
		} else {
			//If study_num is complete quit
			cout << "Study is complete " << endl;
			exit(0);
		}
		//Note make sure h and a is defined
		set_study(); 
	}
	
	void set_study(vector<int> study_num) {
		double Heights[] = {10, 20, 30, 40};
		double Angles[] = {20, 22, 24, 26, 28, 30, 40, 50, 60};
		setInflowHeight(Heights[study_num[1]-1]);
		setChuteAngle(Angles[study_num[2]-1]);
		set_study(study_num[0]);
	}

	//Do not add or remove particles
	void actionsBeforeTimeStep(){ };
		
	//Set up periodic walls, rough bottom, add flow particles
	void setupInitialConditions()
	{
		fix_hgrid();
		set_Nmax(get_N()+getChuteLength()*getChuteWidth()*getZMax());//why is this line needed?
		
		createBottom();
		//~ write(std::cout,false);
		//cout << "correct fixed" << endl;
		if (Species.size()>1) {
			for (int i=0; i<get_N(); i++)
				if (Particles[i].is_fixed()) 
					Particles[i].indSpecies=1;
		}

		set_NWall(1);
		if (getFixedParticleRadius()) {
			wallHandler.getObject(0)->set(Vec3D(0,0,-1), 3.4*MaxInflowParticleRadius);
		} else {
			wallHandler.getObject(0)->set(Vec3D(0,0,-1), 0.);
		}

		set_NWallPeriodic(2);
		WallsPeriodic[0].set(Vec3D( 1.0, 0.0, 0.0), getXMin(), getXMax());
		WallsPeriodic[1].set(Vec3D( 0.0, 1.0, 0.0), getYMin(), getYMax());
		
		add_flow_particles();

		cout << endl << "Status before solve:" << endl;
		cout 
			<< "tc=" << getCollisionTime() 
			<< ", eps="	<< getRestitutionCoefficient()
			<< ", vmax=" << getMaximumVelocity()
			<< ", InflowHeight/zmax=" << getInflowHeight()/getZMax()
			<< endl << endl;
		//~ timer.set(t,tmax);

		//optimize number of buckets
		cout << "Nmax" << get_Nmax() << endl;
		set_HGRID_num_buckets_to_power(get_N()*1.5);
	}

	//add flow particles
	void add_flow_particles() 
	{
                set_HGRID_num_buckets_to_power(get_Nmax());
                hGridActionsBeforeTimeLoop();
                hGridActionsBeforeTimeStep();
                unsigned int N=get_N()+getChuteLength()*getChuteWidth()*InflowHeight;
		set_Nmax(N);
		double H = InflowHeight;
		setZMax(1.2*InflowHeight);
		
		writeRestartFile();
		//try to find new insertable particles
		while (Particles.size()<N){
			create_inflow_particle();
			if (IsInsertable(P0)) {
				num_created++;
			} else InflowHeight += .0001* MaxInflowParticleRadius;
		}
		InflowHeight = H;
		set_HGRID_num_buckets_to_power();
		write(std::cout,false);
	}
	
	//defines type of flow particles	
	void create_inflow_particle()
	{
		P0.Radius = MaxInflowParticleRadius;
		P0.computeMass(Species);
		
		P0.Position.X = random(getXMin()+2.0*P0.Radius,getXMax());
		P0.Position.Y = random(getYMin()+2.0*P0.Radius,getYMax());
		P0.Position.Z = random(getZMin()+2.0*P0.Radius,getInflowHeight());
		P0.Velocity = Vec3D(0.0,0.0,0.0);
	}

	//set approximate height of flow
	void set_H(double new_) {InflowHeight=new_; setZMax(InflowHeight);}
	double get_H() {return InflowHeight;}

	void printTime() {
		cout << "t=" << setprecision(3) << left << setw(6) << getTime() 
			<< ", tmax=" << setprecision(3) << left << setw(6) << getTimeMax()
			<< ", N=" << setprecision(3) << left << setw(6) << Particles.size()
			//<< ", time left=" << setprecision(3) << left << setw(6) << timer.getTime2Finish(t)
			//~ << ", finish by " << setprecision(3) << left << setw(6) << timer.getFinishTime(t)
			<< ". " << endl;
	}
	
	int readNextArgument(unsigned int& i, unsigned int& argc, char *argv[]) {
		if (!strcmp(argv[i],"-muBottom")) {
			setSlidingFrictionCoefficientBottom(atof(argv[i+1]));
			cout << "muB=" << getSlidingFrictionCoefficientBottom() << endl;
		} else return Chute::readNextArgument(i, argc, argv); //if argv[i] is not found, check the commands in Chute
		return true; //returns true if argv[i] is found
	}
};
