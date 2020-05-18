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

#include <string.h>
#include "Chute.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"

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
        baseSpecies = nullptr;
        species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
		species->setDensity(6/constants::pi);
		//~ setStiffnessAndRestitutionCoefficient(2e5,0.97,1);
		double tc=5e-3, r=0.97, beta=0.44, mu=0.092, mur=0.042;
		species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc,r,beta,1.0);// need to consider effective mass
		species->setSlidingFrictionCoefficient(mu);
        species->setRollingFrictionCoefficient(mur);
		
		//chute properties
		setChuteAngleAndMagnitudeOfGravity(24.0, 1.0);
		setChuteLength(20);
		setChuteWidth(10);
		set_H(20);
		nCreated_=0;
			
	}
	
	//void fix_hgrid() {
		//assume 1-2 levels are optimal (which is the case for mono and bidispersed) and set the cell size to min and max 
		// !this is not optimal for polydispersed
	//	double minCell = 2.*min(getFixedParticleRadius(),getMinInflowParticleRadius());
	//	double maxCell = 2.*max(getFixedParticleRadius(),getMaxInflowParticleRadius());
	//	if ((minCell==maxCell)|(minCell==0.)) set_HGRID_max_levels(1);
	//	else set_HGRID_max_levels(2);
	//	set_HGRID_cell_to_cell_ratio (1.0000001*maxCell/minCell);
	//}

    Mdouble getSlidingFrictionCoefficientBottom() {
        if (baseSpecies!= nullptr)
            return baseSpecies->getSlidingFrictionCoefficient();
        else return species->getSlidingFrictionCoefficient();
    }
    void setSlidingFrictionCoefficientBottom(Mdouble new_) {
        createBaseSpecies();
        baseSpecies->setSlidingFrictionCoefficient(new_);
    }
    virtual void createBaseSpecies() {
        //only create once
        static bool created=false;
        if (!created) {
            auto species1 = speciesHandler.copyAndAddObject(species);
            baseSpecies = speciesHandler.getMixedObject(species, species1);
            for (unsigned int i=0; i<particleHandler.getNumberOfObjects(); i++) {
                if (particleHandler.getObject(i)->isFixed()) particleHandler.getObject(i)->setSpecies(speciesHandler.getObject(1));
            }
        }
    }

	void set_study() {
	  std::stringstream name;
		name << "H" << getInflowHeight() 
			 << "A" << getChuteAngleDegrees() 
			 << "L" << round(100.*getFixedParticleRadius()*2.)/100.
		     << "M" << species->getSlidingFrictionCoefficient()
			 << "B" << getSlidingFrictionCoefficientBottom();
		dataFile.setName(name.str().c_str());
		//set_data_filename();
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
		  std::cout << "Study is complete " << std::endl;
			exit(0);
		}
		//Note make sure h and a is defined
		set_study(); 
	}
	
	void set_study(std::vector<int> study_num) {
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
	  //fix_hgrid();
	  //set_Nmax(particleHandler.getNumberOfObjects()+getChuteLength()*getChuteWidth()*getZMax());//why is this line needed?
		
		createBottom();
		//~ write(std::cout,false);
		//cout << "correct fixed" << endl;
		if (speciesHandler.getNumberOfObjects()>1) {
			for (int i=0; i<particleHandler.getNumberOfObjects(); i++)
			  if (particleHandler.getObject(i)->isFixed())
			    particleHandler.getObject(i)->setSpecies(speciesHandler.getObject(1));
		}

		//set_NWall(1);
		InfiniteWall w0;
		if (getFixedParticleRadius()) {
			w0.set(Vec3D(0,0,-1), Vec3D(0,0,-3.4* getMaxInflowParticleRadius()));
		} else {
			w0.set(Vec3D(0,0,-1), Vec3D(0,0,0));
		}
		wallHandler.copyAndAddObject(w0);

		PeriodicBoundary b0;//set_NWallPeriodic(2);
		b0.set(Vec3D( 1.0, 0.0, 0.0), getXMin(), getXMax());
         	boundaryHandler.copyAndAddObject(b0);
		b0.set(Vec3D( 0.0, 1.0, 0.0), getYMin(), getYMax());
		boundaryHandler.copyAndAddObject(b0);
		add_flow_particles();

//		std::cout << std::endl << "Status before solve:" << std::endl;
//		std::cout
//			<< "tc=" << getCollisionTime()
//			<< ", eps=" << getRestitutionCoefficient()
//		  //<< ", vmax=" << getMaximumVelocity()
//			<< ", getInflowHeight()/zmax=" << getInflowHeight()/getZMax()
//			<< std::endl << std::endl;
		//~ timer.set(t,tmax);

		//optimize number of buckets
		//std::cout << "Nmax" << get_Nmax() << std::endl;
		//setHGridNumberOfBucketsToPower(particleHandler.getNumberOfObjects()*1.5);
	}

	//add flow particles
	void add_flow_particles() 
	{
	        //setHGridNumberBucketsToPower(get_Nmax());
                hGridActionsBeforeTimeLoop();
                hGridActionsBeforeTimeStep();
                unsigned int N=particleHandler.getNumberOfObjects()+getChuteLength()*getChuteWidth()*getInflowHeight();
		//set_Nmax(N); // automated in the new version
		double H = getInflowHeight();
		setZMax(1.2*getInflowHeight());
		
		//writeRestartFile();
		//try to find new insertable particles
		while (particleHandler.getNumberOfObjects()<N){
			create_inflow_particle();
			if (checkParticleForInteraction(inflowParticle_)) {
                increaseNCreated();
			} else setInflowHeight(getInflowHeight() + .0001* getMaxInflowParticleRadius());
		}
		setInflowHeight(H);
		//setHGridNumberOfBucketsToPower();
		write(std::cout,false);
	}
	
	//defines type of flow particles	
	void create_inflow_particle()
	{
	    inflowParticle_.setRadius(getMaxInflowParticleRadius());
		//inflowParticle_.computeMass();
		Vec3D position;
		position.X = random.getRandomNumber(getXMin()+2.0*inflowParticle_.getRadius(),getXMax());
		position.Y = random.getRandomNumber(getYMin()+2.0*inflowParticle_.getRadius(),getYMax());
		position.Z = random.getRandomNumber(getZMin()+2.0*inflowParticle_.getRadius(),getInflowHeight());
		inflowParticle_.setPosition(position);
		inflowParticle_.setVelocity(Vec3D(0.0,0.0,0.0));
	}

	//set approximate height of flow
	void set_H(double new_) {setInflowHeight(new_); setZMax(getInflowHeight());}
	double get_H() {return getInflowHeight();}

	void printTime() const {
	  std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime() 
		     << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
		     << ", N=" << std::setprecision(3) << std::left << std::setw(6) << particleHandler.getNumberOfObjects()
			//<< ", time left=" << setprecision(3) << left << setw(6) << timer.getTime2Finish(t)
			//~ << ", finish by " << setprecision(3) << left << setw(6) << timer.getFinishTime(t)
		     << ". " << std::endl;
	}
	
	bool readNextArgument(int& i, int argc, char *argv[]) {
		if (!strcmp(argv[i],"-muBottom")) {
			setSlidingFrictionCoefficientBottom(atof(argv[i+1]));
			std::cout << "muB=" << getSlidingFrictionCoefficientBottom() << std::endl;
		} else return Chute::readNextArgument(i, argc, argv); //if argv[i] is not found, check the commands in Chute
		return true; //returns true if argv[i] is found
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
public:
    LinearViscoelasticFrictionSpecies* species;
    LinearViscoelasticFrictionMixedSpecies* baseSpecies;
};
