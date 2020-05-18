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
#include "scr/Statistics.h"

using namespace std;

class AngleOfRepose : public ChuteWithHopper {
public:
	void setupInitialConditions()
	{
		createBottom();

		set_NWallPeriodic(0);
		set_NWall(1);
		wallHandler.getObject(0)->set(Vec3D(0.0, 0.0, -1.0), -(getZMin()-getFixedParticleRadius()));
		
		//clean up chute
		for (unsigned int i=0;i<particleHandler.getNumberOfObjects();)
		{
			if (mathsFunc::square(particleHandler.getObject(i)->getPosition().X-getChuteLength()*0.5) + mathsFunc::square(particleHandler.getObject(i)->getPosition().Y-getChuteLength()*0.5) > mathsFunc::square(getChuteLength()*0.5)) 
				particleHandler.removeObject(i);
			else i++;
		}
		
		
		particleHandler.set_StorageCapacity((int)min(getXMax()*getYMax()*getZMax()/mathsFunc::cubic(2.0*getInflowParticleRadius()),1e6));
		set_HGRID_num_buckets_to_power(particleHandler.getStorageCapacity());
		write(std::cout,false);
	}	
	
	void create_inflow_particle()
	{
		P0.setRadius(random.getRandomNumber(0.0,1.0)<0.1?MinInflowParticleRadius:MaxInflowParticleRadius);
		P0.computeMass(Species);
		
		Mdouble gamma = random.getRandomNumber(-InflowRadius,InflowRadius);
		P0.setPosition(Vec3D(
			0.5*getChuteLength() + gamma,
			0.5*getChuteWidth() + random.getRandomNumber(-1.0,1.0)*sqrt(mathsFunc::square(InflowRadius)-mathsFunc::square(gamma)),
			FixedParticleRadius + P0.getRadius() + random.getRandomNumber(0.0,1.0)*Height));
		//~ cout << P0.getPosition() << endl;
	}

	void actionsBeforeTimeStep() 
	{
		if (getTime()<getTimeMax()*0.7) {
			//~ InflowRadius = .6*getChuteLength()*std::max(0.0,0.5-t/tmax);
			InflowRadius = 5.*getInflowParticleRadius();

			int failed = 0;
			
			//try max_failed times to find new insertable particle
			while (failed<=max_failed){
				create_inflow_particle();
				//~ cout << P0.getPosition() << endl;
				if (IsInsertable(P0)) {
					failed = 0; 
					num_created++;
				} else failed++;
			};
		}
		
		//clean up chute
		static int count = 0;
		if (count>10) 
		{
			count = 0;
			for (unsigned int i=0;i<particleHandler.getNumberOfObjects();)
			{
				if (mathsFunc::square(particleHandler.getObject(i)->getPosition().X-getChuteLength()*0.5) + mathsFunc::square(particleHandler.getObject(i)->getPosition().Y-getChuteLength()*0.5) > mathsFunc::square(getChuteLength()*0.5))
				{
					#ifdef DEBUG_OUTPUT_FULL
						cout << "erased:" << i << endl;
					#endif
					particleHandler.removeObject(i);
				} //end if particle can be erased
				else i++;
			}
		} else count++;

		
	}
		
	Mdouble Height;
	Mdouble InflowRadius;
	
	//////////////////following copied from SilbertPeriodic.h
	
	AngleOfRepose() {
		// Problem parameters
		setName("silbert");
		
		//time stepping
		setTimeStep(1e-4);
		setTimeMax(2000);

		//output parameters
		setSaveCount(50e1);
	 
		//particle radii
		setInflowParticleRadius(.5);
		setFixedParticleRadius(.5);//getInflowParticleRadius());
		setRoughBottomType(MULTILAYER);

		//particle properties
		setDensity(6/constants::pi);
		setStiffness(2e5);
		setSlidingStiffness(2.0/7.0*getStiffness());
		setDissipation(25.0);
		setSlidingDissipation(getDissipation());
		setSlidingFrictionCoefficient(0.5);
		
		//chute properties
		setChuteAngle(24.0, 1.0);
		setChuteLength(20);
		setChuteWidth(10);
		set_H(20);
			
	}

	void set_H(Mdouble new_) {InflowHeight=new_; setZMax(InflowHeight);}
	Mdouble get_H() {return InflowHeight;}

	void run(int study_num)
	{
		//Set up a parameter study
		set_study(study_num);
		
		stringstream name;
		name << "AngleOfRepose"
			 << "L" << round(100.*getFixedParticleRadius()*2.)/100.
			 << "M" << getSlidingFrictionCoefficient()
			 << "B" << getSlidingFrictionCoefficientBottom();
		setName(name.str().c_str());
		
		//Save the info to disk
		writeRestartFile();

		//Check if the run has been done before
		if (FileExists(data_filename.str())) {
			//If it has move on to teh next run immedently
			cout << "Run " << getName() << " has already been done " << endl;
		} else {
			//launch this code
			stringstream com("");
			com << "echo started \tstudy_num \t" 
				 << study_num << " \tname \t" 
				 << getName() << " &>>Report_AngleOfRepose";
			int sysret;
			sysret = system(com.str().c_str());
			cout << "started study_num " 
				 << study_num << ", name " 
				 << getName() << endl;

			restartFile.setFileType(FileType::ONE_FILE);
			dataFile.setFileType(FileType::ONE_FILE);
			fStatFile.setFileType(FileType::NO_FILE);
			eneFile.setFileType(FileType::ONE_FILE);
			setSaveCount(1e4);
			solve();

			com.str("");
			com << "echo finished \tstudy_num \t" 
				 << study_num << " \tname \t" 
				 << getName() << " &>>Report_AngleOfRepose";
			sysret = system(com.str().c_str());
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
		if (study_num < 6) {
			// set mu_all = 0.5, vary lambda
			Mdouble Lambdas[] = {0, 3./6., 4./6., 5./6., 1, 2};
			setFixedParticleRadius(Lambdas[study_num]/2.);
			setSlidingFrictionCoefficient(0.5);
		} else if (study_num < 9) {
			// set lambda = 1, vary mu_all
			Mdouble MuAll[] = {0, 1., 1e20};
			setSlidingFrictionCoefficient(MuAll[study_num-6]);
			setFixedParticleRadius(0.5);
		} else if (study_num < 12) {
			// set lambda = 1, mu_all = 0.5, vary mu_bottom
			Mdouble MuBottom[] = {0, 1., 1e20};
			setSlidingFrictionCoefficient(0.5);
			setSlidingFrictionCoefficientBottom(MuBottom[study_num-9]);
			setFixedParticleRadius(0.5);
		} else {
			//If study_num is complete quit
			cout << "Study is complete " << endl;
			exit(0);
		}
		//Note make sure h and a is defined
		set_study(); 
	}
	
	Mdouble getSlidingFrictionCoefficientBottom() { 
		if (speciesHandler.getNumberOfObjects()>1) return speciesHandler.getMixedObject(1,0)->getSlidingFrictionCoefficient(); 
		else return getSlidingFrictionCoefficient(); 
	}
	void setSlidingFrictionCoefficientBottom(Mdouble new_) { createBaseSpecies(); speciesHandler.getMixedObject(1, 0)->setSlidingFrictionCoefficient(new_); }
	void createBaseSpecies() {
		//only create once
		static bool created=false;
		if (!created) {
			speciesHandler.copyAndAddObject(speciesHandler.getObject(0));
			for (unsigned int i=0; i<particleHandler.getNumberOfObjects(); i++) {
				if (particleHandler.getObject(i)->isFixed()) particleHandler.getObject(i)->setSpecies(speciesHandler.getObject(1));
			}
		}
	}

};

//this function has been changed to become a demo code
int main(int argc, char *argv[])
{
	AngleOfRepose problem;

	//Problem parameters
	//~ problem.setTimeMax(1000.0);
	problem.setTimeMax(1.0);
	
	//Chute properties
	problem.setChuteAngle(0.0);
	//~ problem.setChuteLength(200.0);
	problem.setChuteLength(20.0);
	problem.setChuteWidth(problem.getChuteLength());
	problem.Height = max(5.,problem.getChuteLength() * 0.5 * tan(  25.0  *constants::pi/180.0));
	problem.setMaxFailed(1);

	int study_num;
	if (argc>1) {
		study_num=atoi(argv[1]);
	} else {
		study_num=0;
		//~ cout << "Please enter study number" << endl;
		//exit(-1);
	}
	problem.run(study_num);
}
