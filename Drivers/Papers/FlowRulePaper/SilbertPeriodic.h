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

#include <iomanip>
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
		species->setStiffness(2e5);
		species->setSlidingStiffness(2.0/7.0* species->getStiffness());
		species->setDissipation(25.0);
		//setSlidingDissipation(2.0/7.0*getDissipation());
		species->setSlidingDissipation(species->getDissipation());
		species->setSlidingFrictionCoefficient(0.5);
        inflowParticle_.setSpecies(species);
		
		//chute properties
		setChuteAngleAndMagnitudeOfGravity(24.0, 1.0);
		setChuteLength(20);
		setChuteWidth(10);
		setInflowVelocity(0);
		set_H(20);
		
		randomiseSpecies=false;
		nCreated_=0;
        
        restartFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
		dataFile.setFileType(FileType::NO_FILE);
		fStatFile.setFileType(FileType::NO_FILE);
        eneFile.setFileType(FileType::ONE_FILE);
	}
	
	//void fix_hgrid() {
		//assume 1-2 levels are optimal (which is the case for mono and bidispersed) and set the cell size to min and max 
		// !this is not optimal for polydispersed
	//	Mdouble minCell = 2.*min(getFixedParticleRadius(),getMinInflowParticleRadius());
	//	Mdouble maxCell = 2.*max(getFixedParticleRadius(),getMaxInflowParticleRadius());
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
                if (particleHandler.getObject(i)->isFixed())
                	particleHandler.getObject(i)->setSpecies(speciesHandler.getObject(1));
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
		setName(name.str().c_str());
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
		//Case 37 set getSlidingDissipation()=0
		//Case 38 set eps=0.97
		//Case 39 set hertzian = true
		//Case 40, 41, 42: set Foerster glass, Lorenz steel, Lorenz glassv
		//Case 43, 44, 45, 46: set Silbert, Foerster glass, Lorenz steel, Lorenz glass with rolling friction
	  std::cout << "Study " << study_num << std::endl;
		
		if (study_num < 6) {
			// set mu_all = 0.5, vary lambda
			Mdouble Lambdas[] = {0, 3./6., 4./6., 5./6., 1, 2};
			setFixedParticleRadius(Lambdas[study_num]/2.);
			species->setSlidingFrictionCoefficient(0.5);
		} else if (study_num < 9) { //Case 6,7,8
			// set lambda = 1, vary mu_all
			Mdouble MuAll[] = {0, 1., 1e20};
			species->setSlidingFrictionCoefficient(MuAll[study_num-6]);
			setFixedParticleRadius(0.5);
		} else if (study_num < 12) { //Case 9,10,11
			// set lambda = 1, mu_all = 0.5, vary mu_bottom
			Mdouble MuBottom[] = {0, 1., 1e20};
			species->setSlidingFrictionCoefficient(0.5);
			setSlidingFrictionCoefficientBottom(MuBottom[study_num-9]);
			setFixedParticleRadius(0.5);
		} else if (study_num < 14) { //Case 12,13
			// set lambda = 1, mu_all = 0.5, vary mu_bottom
			Mdouble MuBottom[] = {0.25, 0.125};
			species->setSlidingFrictionCoefficient(0.5);
			setSlidingFrictionCoefficientBottom(MuBottom[study_num-12]);
			setFixedParticleRadius(0.5);
		} else if (study_num < 16) { //Case 14,15
			// set lambda = 1, vary mu_all
			Mdouble MuAll[] = {0.25, 0.125};
			species->setSlidingFrictionCoefficient(MuAll[study_num-14]);
			setFixedParticleRadius(0.5);
		} else if (study_num < 21) { //Case 16,17,18,19,20
			// set mu_all = 0.5, vary lambda
			Mdouble Lambdas[] = {1./6., 2./6., 1.5, 4, 1./12};
			setFixedParticleRadius(Lambdas[study_num-16]/2.);
			species->setSlidingFrictionCoefficient(0.5);
		} else if (study_num < 26) { //Case 21 22 23 24 25
			// set lambda = 1, mu_all = 0.5, vary mu_bottom
			Mdouble MuBottom[] = {1./16.,1./32.,1./64.,1./128.,1./1024.};
			species->setSlidingFrictionCoefficient(0.5);
			setSlidingFrictionCoefficientBottom(MuBottom[study_num-21]);
			setFixedParticleRadius(0.5);
		} else if (study_num < 29) { //Case 26 27 28
			// set lambda = 1/2, mu_all = 0.5, vary mu_bottom
			Mdouble MuBottom[] = {1./16.,1./128.,1./1024.};
			species->setSlidingFrictionCoefficient(0.5);
			setSlidingFrictionCoefficientBottom(MuBottom[study_num-26]);
			setFixedParticleRadius(0.25);
		} else if (study_num < 33) { //Case 29 30 31 32
			// set lambda = 0, mu_all = 0.5, vary mu_bottom
			Mdouble MuBottom[] = {1./16.,1./128.,1./1024.,0};
			species->setSlidingFrictionCoefficient(0.5);
			setSlidingFrictionCoefficientBottom(MuBottom[study_num-29]);
			setFixedParticleRadius(0);
		} else if (study_num < 37) { //Case 33-36
		  std::cout << "S" << study_num << std::endl;
			// set lambda = 1, mu_b = 0.5, vary mu
			Mdouble Mu[] = {1e20,1,1./64,0};
			species->setSlidingFrictionCoefficient(Mu[study_num-33]);
			setSlidingFrictionCoefficientBottom(0.5);
			setFixedParticleRadius(1);
		} else if (study_num < 38) { //Case 37
			// set getSlidingDissipation()=0
		  species->setSlidingDissipation(0);			
		} else if (study_num < 39) { //Case 38
			// set eps=0.97
			Mdouble eps = 0.97;
			species->setStiffnessAndRestitutionCoefficient(species->getStiffness(), eps, 1);
			species->setSlidingDissipation(species->getDissipation());			
		} else if (study_num < 40) { //Case 39
			// set hertzian = true
		  std::cout << "Hertzian implementation has been changed" << std::endl;
			exit(-1);
			//set_Hertzian(true);
			///\todo Thomas: Hertzian does not  appear in the restart file
		} else if (study_num < 43) { //Case 40, 41, 42
			// set Foerster glass, Lorenz steel, Lorenz glass
			Mdouble eps[] = {0.97 , 0.95 , 0.972};
			Mdouble beta[]= {0.44 , 0.32 , 0.25 };
			Mdouble mu[]  = {0.092, 0.099, 0.177};
			species->setSlidingFrictionCoefficient(mu[study_num-40]);
			species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient
				(50.*getTimeStep(), eps[study_num-40], beta[study_num-40], 1);
		} else if (study_num < 47) { //Case 43, 44, 45, 46
			// set Silbert, Foerster glass, Lorenz steel, Lorenz glass with rolling friction
			Mdouble eps[] = {0.97 , 0.95 , 0.972};
			Mdouble beta[]= {0.44 , 0.32 , 0.25 };
			Mdouble mu[]  = {0.092, 0.099, 0.177};
			if (study_num!=43) {
			  species->setSlidingFrictionCoefficient(mu[study_num-44]);
			  species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient
					(50.*getTimeStep(), eps[study_num-44], beta[study_num-44], 1);
			}
            ///\todo TW: can we switch species from SlidingFriction to Friction only in the case rolling friction is needed?
			species->setRollingStiffness(0.4*species->getStiffness());
			species->setRollingFrictionCoefficient(0.05);
		} else if (study_num < 48) { //Case 47
			// set lambda = 1, mu_all = 0.5, vary mu_half
			Mdouble MuHalf[] = {0};
			species->setSlidingFrictionCoefficient(MuHalf[study_num-47]);
			setSlidingFrictionCoefficientBottom(0.5);
			setFixedParticleRadius(0.5);
			randomiseSpecies = true;
		} else if (study_num < 52) { //Case 48, 49, 50, 51
			// set vary eps
			Mdouble eps[] = {0.001, 0.01, 0.1, 1};
			species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient
					(50.*getTimeStep(), eps[study_num-48], eps[study_num-48], 1);
			species->setSlidingFrictionCoefficient(0.0);
		} else if (study_num < 53) { //Case 52
			///\todo turn on rolling friction only at the wall
		} else if (study_num < 54) { //Case 53
		    std::cout << "using mu=0.3, r=0.1" << std::endl;			
            species->setSlidingFrictionCoefficient(0.3);
            species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient
                     (50.*getTimeStep(),0.1,0.1,1);
 		} else if (study_num < 55) { //Case 54
		    std::cout << "using mu=0.3, r=0.88" << std::endl;			
            species->setSlidingFrictionCoefficient(0.3);
            species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient
                     (50.*getTimeStep(),0.88,0.88,1);
        } else {
			//If study_num is complete quit
		    std::cout << "Study is complete " << std::endl;
			exit(0);
		}
		//Note make sure h and a is defined
		if (study_num < 37 || (study_num>=53&&study_num<=55)) 
		{
			set_study(); 
		}
		else 
		{
		  std::stringstream name;
			name << "S" << study_num; 
			dataFile.setName(name.str().c_str());
			//set_data_filename();
		}

	}

	//Do not add or remove particles
	virtual void actionsBeforeTimeStep(){ };
		
	//Set up periodic walls, rough bottom, add flow particles
	void setupInitialConditions()
	{
        //fix_hgrid();
		particleHandler.setStorageCapacity(particleHandler.getNumberOfObjects()+getChuteLength()*getChuteWidth()*getZMax());//why is this line needed?
		
		createBottom();
        //~ write(std::cout,false);
		//cout << "correct fixed" << endl;
		if (speciesHandler.getNumberOfObjects()>1) {
			for (unsigned int i=0; i<particleHandler.getNumberOfObjects(); i++)
				if (particleHandler.getObject(i)->isFixed()) 
					particleHandler.getObject(i)->setSpecies(speciesHandler.getObject(1));
		}

        ///todo{I(Dinant) had to clear the WallHandler to prevent it from inserting the same wall twice, why?}
		wallHandler.clear();
		InfiniteWall w0;
		w0.setSpecies(speciesHandler.getObject(0));
		if (getFixedParticleRadius()) {
			w0.set(Vec3D(0,0,-1), Vec3D(0,0,-3.4* getMaxInflowParticleRadius()));
		} else {
			w0.set(Vec3D(0,0,-1), Vec3D(0,0,0));
		}
		wallHandler.copyAndAddObject(w0);

		PeriodicBoundary b0;
		b0.set(Vec3D(1.0,0.0,0.0), getXMin(), getXMax());
		boundaryHandler.copyAndAddObject(b0);
		b0.set(Vec3D(0.0,1.0,0.0), getYMin(), getYMax());
		boundaryHandler.copyAndAddObject(b0);
		
		add_flow_particles();

		std::cout << std::endl << "Status before solve:" << std::endl;
//		std::cout
//			<< "tc=" << getCollisionTime()
//			<< ", eps="	<< getRestitutionCoefficient()
//			<< ", vmax=" << getInflowParticle()->calculateMaximumVelocity(getSpecies())
//			<< ", inflowHeight/zMax=" << getInflowHeight()/getZMax()
//			<< std::endl << std::endl;
		//~ timer.set(t,tmax);

		//optimize number of buckets
		std::cout << "Nmax" << particleHandler.getStorageCapacity() << std::endl;
		//setHGridNumberOfBucketsToPower(particleHandler.getNumberOfObjects()*1.5);
	}

	//add flow particles
	void add_flow_particles() 
	{
	    //setHGridNumberOfBucketsToPower(particleHandler.getStorageCapacity());
		hGridActionsBeforeTimeLoop();
        checkAndDuplicatePeriodicParticles();
        hGridActionsBeforeTimeStep();
		unsigned int N=getChuteLength()*getChuteWidth()* getInflowHeight();
		particleHandler.setStorageCapacity(particleHandler.getNumberOfObjects()+1.5*N);
		Mdouble H = getInflowHeight();
		setZMax(1.0*getInflowHeight());
        //uncomment the following line to achieve a high packing fraction
        setInflowHeight(getZMin()+inflowParticle_.getRadius());

        writeRestartFile();
		//try to find new insertable particles
        while (getNCreated()<N){
			create_inflow_particle();
			if (checkParticleForInteraction(inflowParticle_)) {
                BaseParticle* p = particleHandler.copyAndAddObject(inflowParticle_);
                //duplicate particle
                for (BaseBoundary* it : boundaryHandler)
                    it->createPeriodicParticle(p, particleHandler);
                p = particleHandler.getLastObject();
                //duplicate duplicate particles (this is a hack which is needed as there are two boundaries, so the doubly periodic images are needed)
                for (BaseBoundary* it : boundaryHandler)
                    it->createPeriodicParticle(p, particleHandler);
                //hGridActionsBeforeTimeStep();
                increaseNCreated();
			} else setInflowHeight(getInflowHeight() + .0001* getMaxInflowParticleRadius());
		}
        for (unsigned int i = particleHandler.getNumberOfObjects(); i >= 1; i--)
        if (particleHandler.getObject(i - 1)->getPeriodicFromParticle() != nullptr)
        {
            while (particleHandler.getObject(i - 1)->getInteractions().size() > 0)
            {
                interactionHandler.removeObjectKeepingPeriodics(particleHandler.getObject(i - 1)->getInteractions().front()->getIndex());
            }
            particleHandler.removeObject(i - 1);
        }
		setInflowHeight(H);
		//setHGridNumberOfBucketsToPower();
		write(std::cout,false);
	}
	
	//defines type of flow particles	
	void create_inflow_particle()
	{
		inflowParticle_.setRadius(random.getRandomNumber(getMinInflowParticleRadius(),getMaxInflowParticleRadius()));
		//inflowParticle_.computeMass();
		
        //The position components are first stored in a Vec3D, because if you pass them directly into setPosition the compiler is allowed to change the order in which the numbers are generated
        Vec3D position;
		position.X = random.getRandomNumber(getXMin(),getXMax());
		position.Y = random.getRandomNumber(getYMin(),getYMax());
		position.Z = random.getRandomNumber(getZMin()+inflowParticle_.getRadius(),getInflowHeight());
        inflowParticle_.setPosition(position);
		inflowParticle_.setVelocity(Vec3D(getInflowVelocity(),0.0,0.0));
		if (randomiseSpecies) {
			int indSpecies = floor(random.getRandomNumber(0,speciesHandler.getNumberOfObjects()-1e-200));
			inflowParticle_.setSpecies(speciesHandler.getObject(indSpecies));
		}
	}

	//set approximate height of flow
	void set_H(Mdouble new_) {setInflowHeight(new_); setZMax(getInflowHeight());}
	Mdouble get_H() {return getInflowHeight();}

	void printTime() const {
	  std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime() 
		    << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
		    << ", N=" << std::setprecision(3) << std::left << std::setw(6) << particleHandler.getNumberOfObjects()
		    << ", theta=" << std::setprecision(3) << std::left << std::setw(6) << getChuteAngleDegrees()
			//<< ", time left=" << setprecision(3) << left << setw(6) << timer.getTime2Finish(t)
			//~ << ", finish by " << setprecision(3) << left << setw(6) << timer.getFinishTime(t)
		    << ". " << std::endl;
	}
	
	bool readNextArgument(int& i, int argc, char *argv[]) {
		if (!strcmp(argv[i],"-muBottom")) {
			setSlidingFrictionCoefficientBottom(atof(argv[i+1]));
			std::cout << "muB=" << getSlidingFrictionCoefficientBottom() << std::endl;
		} else if (!strcmp(argv[i],"-oldValues")) {
		  species->setSlidingDissipation(species->getDissipation());
		  std::cout << "getSlidingDissipation()=" << species->getSlidingDissipation() << std::endl;
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
	bool randomiseSpecies;
	SphericalParticle inflowParticle_;

    LinearViscoelasticFrictionSpecies* species;
    LinearViscoelasticFrictionMixedSpecies* baseSpecies;
};
