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

#include "ChuteWithHopper.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "Walls/InfiniteWall.h"

class AngleOfRepose : public ChuteWithHopper {
public:
	void setupInitialConditions() override {
		createBottom();

		//set_NWallPeriodic(0);
		//set_NWall(1);
		///todo{I(Dinant) had to clear the WallHandler to prevent it from inserting the same wall twice, why?}
		wallHandler.clear();
		InfiniteWall w0;
		w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0,0,getZMin()-getFixedParticleRadius()));
		wallHandler.copyAndAddObject(w0);
		
		//clean up chute
		for (unsigned int i=0;i<particleHandler.getNumberOfObjects();)
		{
			if (mathsFunc::square(particleHandler.getObject(i)->getPosition().X-getChuteLength()*0.5) + mathsFunc::square(particleHandler.getObject(i)->getPosition().Y-getChuteLength()*0.5) > mathsFunc::square(getChuteLength()*0.5)) 
				particleHandler.removeObject(i);
			else i++;
		}
		
		
		particleHandler.setStorageCapacity(static_cast<unsigned>(std::min(getXMax()*getYMax()*getZMax()/mathsFunc::cubic(2.0*getInflowParticleRadius()),1e6)));
		//setHGridNumberOfBucketsToPower(particleHandler.getStorageCapacity());
		write(std::cout,false);
		nCreated_=0;
	}	
	
	void create_inflow_particle()
	{
		inflowParticle_.setRadius(random.getRandomNumber(0.0,1.0)<0.1? getMinInflowParticleRadius() : getMaxInflowParticleRadius());
		///\bug This line should be longer be required but if this code still works should be tested.
        //inflowParticle_.computeMass();
		
		Mdouble gamma = random.getRandomNumber(-InflowRadius,InflowRadius);
		inflowParticle_.setPosition(Vec3D(
			0.5*getChuteLength() + gamma,
			0.5*getChuteWidth() + random.getRandomNumber(-1.0,1.0)*sqrt(mathsFunc::square(InflowRadius)-mathsFunc::square(gamma)),
                getFixedParticleRadius() + inflowParticle_.getRadius() + random.getRandomNumber(0.0,1.0)*Height));
		//~ cout << P0.getPosition() << endl;
	}

	void actionsBeforeTimeStep() override {
		if (getTime()<getTimeMax()*0.7) {
			//~ InflowRadius = .6*getChuteLength()*std::max(0.0,0.5-t/tmax);
			InflowRadius = 5.*getInflowParticleRadius();

			int failed = 0;
			
			//try max_failed times to find new insertable particle
			while (failed<= getMaxFailed()){
				create_inflow_particle();
				//~ cout << P0.getPosition() << endl;
				if (checkParticleForInteraction(inflowParticle_)) {
					failed = 0;
                    increaseNCreated();
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
		setSaveCount(500);
	 
		//particle radii
		setInflowParticleRadius(.5);
		setFixedParticleRadius(.5);//getInflowParticleRadius());
		setRoughBottomType(MULTILAYER);

		//particle properties
        baseSpecies = nullptr;
        species = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        species->setDensity(6/constants::pi);
        species->setStiffness(2e5);
		species->setSlidingStiffness(2.0/7.0* species->getStiffness());
		species->setDissipation(25.0);
		species->setSlidingDissipation(species->getDissipation());
		species->setSlidingFrictionCoefficient(0.5);
		
		//chute properties
        setChuteAngleAndMagnitudeOfGravity(24.0, 1.0);
		setChuteLength(20);
		setChuteWidth(10);
		set_H(20);
			
	}

	void set_H(Mdouble new_) {
        setInflowHeight(new_); setZMax(new_);
    }

	Mdouble get_H() {return getInflowHeight();}

	void run(int study_num)
	{
		//Set up a parameter study
		set_study(study_num);
		
		std::stringstream name;
		name << "AngleOfRepose"
			 << "L" << round(100.*getFixedParticleRadius()*2.)/100.
		     << "M" << species->getSlidingFrictionCoefficient()
			 << "B" << getSlidingFrictionCoefficientBottom();
		setName(name.str().c_str());
		
		//Save the info to disk
		writeRestartFile();

		//Check if the run has been done before
		if (helpers::fileExists(dataFile.getName())) {
			//If it has move on to teh next run immedently
		  std::cout << "Run " << getName() << " has already been done " << std::endl;
		} else {
			//launch this code
		  std::stringstream com("");
			com << "echo started \tstudy_num \t" 
				 << study_num << " \tname \t" 
				 << getName() << " &>>Report_AngleOfRepose";
			///todo{Change this line to have some meaningfull behaviour}
			if(system(com.str().c_str())){}
			std::cout << "started study_num " 
				 << study_num << ", name " 
			     << getName() << std::endl;

			restartFile.setFileType(FileType::ONE_FILE);
			dataFile.setFileType(FileType::ONE_FILE);
			fStatFile.setFileType(FileType::ONE_FILE);
			eneFile.setFileType(FileType::ONE_FILE);
			setSaveCount(10000);
			solve();

			com.str("");
			com << "echo finished \tstudy_num \t" 
				 << study_num << " \tname \t" 
				 << getName() << " &>>Report_AngleOfRepose";
			///todo{Change this line to have some meaningfull behaviour}
			if(system(com.str().c_str())){}
		}
	}
	
	void set_study() {
	  std::stringstream name;
		name << "H" << getInflowHeight() 
			 << "A" << getChuteAngleDegrees() 
			 << "L" << round(100.*getFixedParticleRadius()*2.)/100.
		     << "M" << species->getSlidingFrictionCoefficient()
			 << "B" << getSlidingFrictionCoefficientBottom();
		dataFile.setName(name.str());
		//set_data_filename();
	}

	void set_study(int study_num) {
		if (study_num < 6) {
			// set mu_all = 0.5, vary lambda
			Mdouble Lambdas[] = {0, 3./6., 4./6., 5./6., 1, 2};
			setFixedParticleRadius(Lambdas[study_num]/2.);
			species->setSlidingFrictionCoefficient(0.5);
		} else if (study_num < 9) {
			// set lambda = 1, vary mu_all
			Mdouble MuAll[] = {0, 1., 1e20};
			species->setSlidingFrictionCoefficient(MuAll[study_num-6]);
			setFixedParticleRadius(0.5);
		} else if (study_num < 12) {
			// set lambda = 1, mu_all = 0.5, vary mu_bottom
			Mdouble MuBottom[] = {0, 1., 1e20};
			species->setSlidingFrictionCoefficient(0.5);
			setSlidingFrictionCoefficientBottom(MuBottom[study_num-9]);
			setFixedParticleRadius(0.5);
		} else {
			//If study_num is complete quit
		  std::cout << "Study is complete " << std::endl;
			exit(0);
		}
		//Note make sure h and a is defined
		set_study(); 
	}
	
	Mdouble getSlidingFrictionCoefficientBottom() { 
		if (baseSpecies!= nullptr)
            return baseSpecies->getSlidingFrictionCoefficient();
		else return species->getSlidingFrictionCoefficient(); 
	}
	void setSlidingFrictionCoefficientBottom(Mdouble new_) {
        createBaseSpecies();
        baseSpecies->setSlidingFrictionCoefficient(new_);
    }
	void createBaseSpecies() {
		//only create once
		static bool created=false;
		if (!created) {
			auto species1 = speciesHandler.copyAndAddObject(species);
            baseSpecies = speciesHandler.getMixedObject(species, species1);
			for (unsigned int i=0; i<particleHandler.getNumberOfObjects(); i++) {
				if (particleHandler.getObject(i)->isFixed())
					particleHandler.getObject(i)->setSpecies(speciesHandler.getObject(1));
			}
		}
	}

    int getNCreated() const
    {
        return nCreated_;
    }

    void increaseNCreated()
    {
        nCreated_++;
    }

    int nCreated_;
	SphericalParticle inflowParticle_;
    LinearViscoelasticSlidingFrictionSpecies* species;
    LinearViscoelasticSlidingFrictionMixedSpecies* baseSpecies;
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
	problem.Height = std::max(5.,problem.getChuteLength() * 0.5 * tan(  25.0  *constants::pi/180.0));
	problem.setMaxFailed(1);

	int study_num;
	if (argc>1) {
		study_num=atoi(argv[1]);
	} else {
		study_num=0;
		std::cout << "Please enter study number" << std::endl;
		exit(-1);
	}
	problem.run(study_num);
}
