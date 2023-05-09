#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Mercury3D.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"
#include <chrono>

class DrumRot : public Mercury3D
{
	public:

	void setupInitialConditions() override
	{
		radiusS1 = 0.0005; // 1mm diameter

		rhoS1 = 2500.0;

		massS1 = 4 / 3 * constants::pi * pow(radiusS1, 3.0) * rhoS1;

		double fillVolume = drumFillFraction*constants::pi*pow(drumRadius,2.0)*(std::abs(getYMax()-getYMin()));
		numS1 = (fillVolume )/(4./3. * constants::pi*pow(radiusS1,3.0)); // Drum volume - volume large particles
		//std::cout << "fillVolume" << fillVolume << "total particle volume" << numS1*4 / 3 * constants::pi * pow(radiusS1, 3.0) + numS2*4 / 3 * constants::pi * pow(radiusS2, 3.0) << std::endl;
        //std::cout<< "fillVolume = " << fillVolume << ", drumFillFraction = "<< drumFillFraction << ", S2Volume = " << S2Volume << std::endl;
        std::cout<< "Expected numS1 = "<<numS1<<std::endl;
		numS1ToBeInserted = numS1;

		//std::cout << " mass " << massS1 << " " << massS2 << std::endl;

		tc = 1 / 4000.0;
		//original value
		//tc = 0.005;

		speciesHandler.clear();

		auto speciesDrum = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
		auto speciesS1 =  speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());

		double RPSInitial = 0.0;

		speciesDrum->setDensity(rhoS1);
		speciesDrum->setCollisionTimeAndRestitutionCoefficient(tc, CORDrum, massS1);
		speciesDrum->setSlidingDissipation(speciesDrum->getDissipation()*2./7.);
		speciesDrum->setSlidingStiffness(speciesDrum->getStiffness()*2./7.);
		speciesDrum->setSlidingFrictionCoefficient(slidingFrictionDrum);
		speciesDrum->setRollingStiffness(speciesDrum->getStiffness()*2.0/7.0);
		speciesDrum->setRollingFrictionCoefficient(rollingFrictionDrum);
		speciesDrum->setRollingDissipation(speciesDrum->getDissipation()*2./7.);
		speciesDrum->setTorsionStiffness(speciesDrum->getStiffness()*2.0/7.0);
		speciesDrum->setTorsionFrictionCoefficient(torsionFrictionDrum);
		speciesDrum->setTorsionDissipation(speciesDrum->getDissipation()*2./7.);
		//

		//
		speciesS1->setDensity(rhoS1);
		speciesS1->setCollisionTimeAndRestitutionCoefficient(tc, CORS1, massS1);
		speciesS1->setSlidingDissipation(speciesS1->getDissipation()*2./7.);
		speciesS1->setSlidingStiffness(speciesS1->getStiffness()*2./7.);
		speciesS1->setSlidingFrictionCoefficient(slidingFriction1);
		speciesS1->setRollingStiffness(speciesS1->getStiffness()*2.0/7.0);
		speciesS1->setRollingFrictionCoefficient(rollingFriction1);
		speciesS1->setRollingDissipation(speciesS1->getDissipation()*2./7.);
		speciesS1->setTorsionStiffness(speciesS1->getStiffness()*2.0/7.0);
		speciesS1->setTorsionFrictionCoefficient(torsionFriction1);
		speciesS1->setTorsionDissipation(speciesS1->getDissipation()*2./7.);
		//

		auto speciesDrumAndS1 = speciesHandler.getMixedObject(speciesDrum,speciesS1);
		speciesDrumAndS1->setCollisionTimeAndRestitutionCoefficient(tc, ((CORS1 + CORDrum) / 2) , massS1, massS1);
		speciesDrumAndS1->setSlidingDissipation(speciesDrumAndS1->getDissipation()*2./7.);
		speciesDrumAndS1->setSlidingFrictionCoefficient( ((slidingFrictionDrum + slidingFriction1)/2));
		speciesDrumAndS1->setSlidingStiffness(speciesDrumAndS1->getStiffness()*2.0/7.0);
		speciesDrumAndS1->setRollingStiffness(speciesDrumAndS1->getStiffness()*2.0/7.0);
		speciesDrumAndS1->setRollingFrictionCoefficient(((rollingFrictionDrum + rollingFriction1)/2));
		speciesDrumAndS1->setRollingDissipation(speciesDrumAndS1->getDissipation()*2./7.);
		speciesDrumAndS1->setTorsionStiffness(speciesDrumAndS1->getStiffness()*2.0/7.0);
		speciesDrumAndS1->setTorsionFrictionCoefficient(((torsionFrictionDrum + torsionFriction1)/2));
		speciesDrumAndS1->setTorsionDissipation(speciesDrumAndS1->getDissipation()*2./7.);
		//


		//
		Vec3D drumCenter = {0.5*(getXMin() + getXMax()),
							0.5*(getYMin() + getYMax()),
							0.5*(getZMin() + getZMax())};

		wallHandler.clear();

		auto drumWall = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
		drumWall->setSpecies(speciesDrum);
		drumWall->setPosition(drumCenter);
		drumWall->setOrientation(Vec3D(0.0,1.0,0.0));
		drumWall->addObject(Vec3D(1,0,0), Vec3D(drumRadius,0.0,0.0));
		drumWall->setAngularVelocity(Vec3D(0.0,RPSInitial * 2.0 * constants::pi,0.0));

		InfiniteWall w0;
		w0.setSpecies(speciesDrum);
		w0.set(Vec3D(0.,-1.,0.),Vec3D(drumCenter.X,getYMin(),drumCenter.Z));
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D(0.,1.,0.),Vec3D(drumCenter.X,getYMax(),drumCenter.Z));
		wallHandler.copyAndAddObject(w0);

		SphericalParticle p0;
		double radius = 0.0;
		int numS1Inserted=0;
		Vec3D pos;
		double R, C;
		int failCounter = 0;

		while (numS1Inserted < numS1)
		{
            failCounter = 0;
            do
            {
                p0.setSpecies(speciesS1);
                p0.setRadius(radiusS1);
                R = random.getRandomNumber(2*radiusS1,drumRadius-2.0*p0.getRadius()); // Find a random radial position within the cylinder
                C = random.getRandomNumber(0.0*constants::pi,2.0*constants::pi); // Find a random angular position within the cylinder
                // Set the position in the cylinder with a random z-coordinate
                p0.setPosition(Vec3D(R*cos(C),random.getRandomNumber(getYMin()+p0.getRadius(),getYMax()-p0.getRadius()),R*sin(C)));
                p0.setVelocity(Vec3D(0.0,0.0,0.0)); // Set initial particle velocity
                failCounter++; // Adding 1 to failcounter
                    if (failCounter==1000)
                    {
                        break; // If we failed 1000 times to find a non-contact position we break the placement loop (makes sure that no infinite loops are created in combination with line 130)
                    }
            } while(!checkParticleForInteraction(p0));
            particleHandler.copyAndAddObject(p0);
            numS1Inserted++;
            hGridRebuild();
		}
		std::cout << "Number of S1 particles inserted" << numS1Inserted << std::endl;

		if ( (numS1Inserted == numS1) )
		{
			step = 2;
			checkTime = getTime() + 1.0;
			std::cout << "\n \n \n";
			std::cout << "Particles settling down, checkTime = " << checkTime << std::endl;
			std::cout << "--------------------------" << std::endl;
			std::cout << "\n\n\n";
		}
	}

	void actionsBeforeTimeStep() override
	{
		//wallHandler.getObject(0)->setOrientation(Vec3D(0.0,1.0,0.0));
		//wallHandler.getObject(1)->setOrientation(Vec3D(0.0,-1.0,0.0));
		//wallHandler.getObject(2)->setOrientation(Vec3D(0.0,1.0,0.0));

		if (step==2)
		{
			if (getTime() > checkTime)
			{
				std::cout << "Current KE" << getKineticEnergy() << std::endl;
				if (getKineticEnergy() < (10.0))
				{
					step = 3;
                    setTimeMax(100.0);
                    setSaveCount(2000);
					double drumStartTime = getTime();
					std::cout << "\n \n \n";
					std::cout << "Start the drum rotation, rpm = " << revolutionsPerSecond*60.0 << std::endl;
					std::cout << "--------------------------" << std::endl;
					std::cout << "\n\n\n";
					// rotate the drum
					wallHandler.getObject(0)->setAngularVelocity(Vec3D(0.0,revolutionsPerSecond*constants::pi*2.0,0.0));
					wallHandler.getObject(1)->setAngularVelocity(Vec3D(0.0,revolutionsPerSecond*constants::pi*2.0,0.0));
					wallHandler.getObject(2)->setAngularVelocity(Vec3D(0.0,revolutionsPerSecond*constants::pi*2.0,0.0));

					//wallHandler.getObject(0)->setOrientation(Vec3D(0.0,1.0,0.0));
					//wallHandler.getObject(1)->setOrientation(Vec3D(0.0,-1.0,0.0));
					//wallHandler.getObject(2)->setOrientation(Vec3D(0.0,1.0,0.0));
				}
				else
				{
					checkTime = getTime() + 0.1;
				}
			}
		}
	}

	void actionsOnRestart() override
	{

        //setTimeMax(200.0);
	    setSaveCount(2000);

	    std::cout << "\n \n \n";
		std::cout << "Set drum rotation speed (rpm) to " << revolutionsPerSecond*60.0 << std::endl;
		std::cout << "--------------------------" << std::endl;
		std::cout << "\n\n\n";

        wallHandler.getObject(0)->setAngularVelocity(Vec3D(0.0,revolutionsPerSecond*constants::pi*2.0,0.0)); //rad/s
		wallHandler.getObject(1)->setAngularVelocity(Vec3D(0.0,revolutionsPerSecond*constants::pi*2.0,0.0));
		wallHandler.getObject(2)->setAngularVelocity(Vec3D(0.0,revolutionsPerSecond*constants::pi*2.0,0.0));
	}

	void setDrumRadius (double radius)
	{
		drumRadius = radius;
	}

	void setRevolutionSpeed (double rpm)
	{
		revolutionsPerSecond = rpm/60.0;// non-dimensionalised based on 3 mm particles and g=9.81
	}

	void setFractionalPolydispersity(double fpd)
	{
		fractionalPolydispersity = fpd;
	}

	void setDrumFillFraction (double dff)
	{
		drumFillFraction = dff;
	}

	void setFrictionCoeff(double pwf, double ppf)
	{
		particleWallFriction = pwf;
		particleParticleFriction = ppf;
	}

	void setCOR (double drumCOR, double COR1)
	{
		CORDrum = drumCOR;
		CORS1 = COR1;
	}
	//a series of functions by which to easily set particles'
	//various frictional coefficients
	void setSlidingFriction (double drum, double f1)
	{
		slidingFrictionDrum = drum;
		slidingFriction1 = f1;
	}

	void setRollingFriction (double drum, double f1)
	{
		rollingFrictionDrum = drum;
		rollingFriction1 = f1;
	}

	void setTorsionFriction (double drum, double f1)
	{
		torsionFrictionDrum = drum;
		torsionFriction1 = f1;
	}

	double getDrumRadius()
	{
		return drumRadius;
	}

	//new functions to set the frequency and amplitude with which the drum is
	//vibrated

	private:

	double radiusS1;
	double rhoS1;
	double massS1;

	double CORDrum,CORS1,tc;

	double sizeRatio;
	double densityRatio;
	double drumFillFraction;
	double volumeFraction;

	double particleWallFriction,particleParticleFriction;

	int numS1;
	int numS1ToBeInserted;

	double drumRadius;
	double revolutionsPerSecond;
	double fractionalPolydispersity;

	double slidingFrictionDrum;
	double slidingFriction1;

	double rollingFrictionDrum;
	double rollingFriction1;

	double torsionFrictionDrum;
	double torsionFriction1;

	int step;
	double checkTime;

};

int main(int argc, char *argv[])
{

    // Start measuring elapsed time
    //std::chrono::time_point<std::chrono::system_clock> startClock, endClock;
    //startClock = std::chrono::system_clock::now();

	DrumRot problem;

	//setting locally-available variables to define key
	//parameters such that they can be automatically included in
	//the name of the file
	//*******************Excitation Properties*****************************************
	//Drum rotation rate (rpm)
	double rotRateBasal = 1.0;

	//*******************Frictional Properties*****************************************
	//sliding friction for species 1, 2 and wall
	double muSWall = 0.61;
	double muS1 = 0.19;

	//rolling friction for species 1, 2 and wall
	double muRWall = 0.06;
	double muR1 = 0.01;

	//torsion friction for species 1, 2 and wall
	double muTWall = 0.0;
	double muT1 = 0.000;

	double dimDrumRad = 45.0;
	double drumRad = dimDrumRad*0.0005; // 45 * smallest particle radius -> dimDrumRad particles fit in drum diameter
	//the dimensionless drum length (L/d_l)
	double dimDrumLength = 10.0;
	double drumLength = dimDrumLength*0.001; // dimDrumLength * smallest particle diameter

    // computes rpm based on froude number and drumRad
    problem.setTimeMax(150.0);
	problem.setTimeStep(1.0/(4000.0 * 50.0));
	double froudeNumber = 0.22;
	double rotRate = pow(froudeNumber*9.81/drumRad,0.5);//*(60.0/(2.0*3.1415926535));
    std::cout<< "rotRate = " << rotRate << std::endl;
	//Set the number of domains for parallel decomposition
	//problem.setNumberOfDomains({2,1,1});

	//*******************Setting Up Filename******************************************
	double dRatio;
	double rhoRatio;

	if (argc < 2)
    {
        std::stringstream nameStream;
        std::string nameBase = "SBead_";
        nameStream << nameBase << "Monodispersed" << "_R" << dimDrumRad << "_L" << dimDrumLength;
        std::string fullName = nameStream.str();
        problem.setName(fullName);
        problem.incrementRunNumberInFile();

    }

    if (argc>1)
    {
        //std::stringstream nameStream;
        //std::string fullName = nameStream.str();
        //fullName = argv[1];
        //std::cout << "restarting from " << fullName << std::endl;
        //problem.readRestartFile(fullName);
        //problem.setAppend(true);
    }
	//problem.random.randomise();

    //problem.autoNumber();

	problem.setDrumRadius(drumRad);// in meters

	problem.setXMin(-problem.getDrumRadius());
	problem.setYMin(-0.5*drumLength);
	problem.setZMin(-problem.getDrumRadius());

	problem.setXMax(problem.getDrumRadius());
	problem.setYMax(0.5*drumLength);// in meters
	problem.setZMax(problem.getDrumRadius());

	problem.setGravity(Vec3D(0.,0.,-9.81));
	problem.setCOR(0.831,0.831);//drumWall, species1, species2
	problem.setFractionalPolydispersity(0.0);//10% dispersity
	problem.setDrumFillFraction(0.23);// At 0.5 the drum is 3/4 filled.
	//redundant
	//problem.setFrictionCoeff(.6,.19);//(particle-wall-friction,part-part-friction)
	problem.setSlidingFriction(muSWall,muS1); //wall, species1, species2
	problem.setRollingFriction(muRWall,muR1); //wall, species1, species2
	problem.setTorsionFriction(muTWall,muT1); //wall, species1, species2

	problem.setRevolutionSpeed(rotRate);//rpm

	problem.setSaveCount(2000);
	problem.readArguments(argc,argv);

    problem.dataFile.setFileType(FileType::MULTIPLE_FILES);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.eneFile.setFileType(FileType::NO_FILE);

	problem.solve();

	return 0;
}
