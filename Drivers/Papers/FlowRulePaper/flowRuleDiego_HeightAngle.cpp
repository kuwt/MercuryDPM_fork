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
#include <cstring>

#include "Chute.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticSpecies.h"

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
        species = speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());
        species->setDensity(6/constants::pi);
		species->setStiffness(2e5);
		species->setDissipation(25.0);
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
		dataFile.setFileType(FileType::ONE_FILE);
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

	void set_study(int study_num) {
	    std::cout << "using mu=0, r=0.5" << std::endl;			
	    species->setCollisionTimeAndRestitutionCoefficient
	             (50.*getTimeStep(),0.5,1);
	  	std::stringstream name;
		name << "H" << getInflowHeight() 
			 << "A" << getChuteAngleDegrees() 
			 << "L" << round(100.*getFixedParticleRadius()*2.)/100.;
		setName(name.str().c_str());
	}

	//Do not add or remove particles
	void actionsBeforeTimeStep() override { };
		
	//Set up periodic walls, rough bottom, add flow particles
	void setupInitialConditions() override {
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
            while (!particleHandler.getObject(i - 1)->getInteractions().empty())
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
			const unsigned int indSpecies = floor(random.getRandomNumber(0, speciesHandler.getNumberOfObjects() - 1e-200));
			inflowParticle_.setSpecies(speciesHandler.getObject(indSpecies));
		}
	}

	//set approximate height of flow
	void set_H(Mdouble new_) {setInflowHeight(new_); setZMax(getInflowHeight());}
	Mdouble get_H() {return getInflowHeight();}

	void printTime() const override {
	  std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime() 
		    << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
		    << ", N=" << std::setprecision(3) << std::left << std::setw(6) << particleHandler.getNumberOfObjects()
		    << ", theta=" << std::setprecision(3) << std::left << std::setw(6) << getChuteAngleDegrees()
			//<< ", time left=" << setprecision(3) << left << setw(6) << timer.getTime2Finish(t)
			//~ << ", finish by " << setprecision(3) << left << setw(6) << timer.getFinishTime(t)
		    << ". " << std::endl;
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

    LinearViscoelasticSpecies* species;
    LinearViscoelasticMixedSpecies* baseSpecies;
};

class FlowRule : public SilbertPeriodic {
public:   
    
	void launch(bool startfast=false) {
		launchNewRun("flowRule",startfast);
	}

    void run(std::vector<Mdouble> studyNumber, int argc, char *argv[])
	{
		//Set up a parameter study
		setInflowHeight(studyNumber[1]);
		setChuteAngle(studyNumber[2]);
		set_study(studyNumber[0]);          
		readArguments(argc, argv);
		//set_study();	
		
		//Save the info to disk
		writeRestartFile();

		//Check if the run has been done before
		//if (helpers::fileExists(dataFile.getName())) {
		//	//If it has move on to teh next run immedently
        //    std::cout << "Run " << getName() << " has already been done " << std::endl;
		//} else 
        {
			//launch this code
            std::stringstream com("");
			com << "echo started \tstudyNumber \t" 
				 << studyNumber[0] << " " 
				 << studyNumber[1] << " " 
				 << studyNumber[2] << " \tname \t" 
				 << getName() << " &>>ReportFlowRule";
			//std::cout << system(com.str().c_str()) << std::endl;
			std::cout << "started studyNumber " 
				 << studyNumber[0] << " " 
				 << studyNumber[1] << " " 
				 << studyNumber[2] << ", name " 
				  << getName() << std::endl;

			solve();

			com.str("");
			com << "echo finished \tstudyNumber \t" 
				 << studyNumber[0] << " " 
				 << studyNumber[1] << " " 
				 << studyNumber[2] << " \tname \t" 
				 << getName() << " &>>ReportFlowRule";
			//std::cout << system(com.str().c_str()) << std::endl;
		}
	}
};

int main(int argc, char *argv[])
{
	FlowRule problem;
	//set case, height, angle to given or default values
	std::vector<Mdouble> studyNumber;
	studyNumber.resize(3);

	//this line is needed for the code to work with demoparams
	//if (argc<4) exit(-1);
	
	if (argc>2) 
    {
        studyNumber[0]=0;
        studyNumber[1]=atof(argv[1]);
        studyNumber[2]=atof(argv[2]);
        problem.run(studyNumber,argc-2,argv+2);
    } 
    else
    {
        std::cout << "Not enough input arguments given (./flowRule_StudyHeightAngle $study $height $angle); " << std::endl
            << "using demo values (equivalent to ./flowRule_StudyHeightAngle 5 10 24 -tmax 0.01)" << std::endl;
        studyNumber[0]=5;
        studyNumber[1]=2;
        studyNumber[2]=24;
        problem.setTimeMax(0.01);
        problem.dataFile.setFileType(FileType::ONE_FILE);
        problem.setChuteLength(5);
        problem.setChuteWidth(5);
        //problem.setRoughBottomType(MULTILAYER);
        problem.run(studyNumber,1,argv);
        problem.setName("flowRuleSelfTest");
        problem.writeRestartFile();
    }
}
