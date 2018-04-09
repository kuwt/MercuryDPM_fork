#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Mercury3D.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"

class RotatingDrum : public Mercury3D
{
	public:
		
	RotatingDrum()
	{	
		radiusS1 = 0.003; // 3mm diameter
		fractionalPolydispersity = 0.0;
	}

	void setupInitialConditions()
	{

		radiusS2 = sizeRatio*radiusS1;
		
		rhoS1 = 2500.0;
		rhoS2 = densityRatio*rhoS1;

		massS1 = 4 / 3 * constants::pi * pow(radiusS1, 3.0) * rhoS1;
		massS2 = 4 / 3 * constants::pi * pow(radiusS2, 3.0) * rhoS2;

		double fillVolume = drumFillFraction*constants::pi*pow(drumRadius,2.0)*(std::abs(getYMax()-getYMin()));

		numS1 = volumeFraction*fillVolume/(4./3. * constants::pi*pow(radiusS1,3.0));
		numS2 = (1. - volumeFraction)*fillVolume/(4./3. * constants::pi*pow(radiusS2,3.0));
		//std::cout << "fillVolume" << fillVolume << "total particle volume" << numS1*4 / 3 * constants::pi * pow(radiusS1, 3.0) + numS2*4 / 3 * constants::pi * pow(radiusS2, 3.0) << std::endl;
		
		numS1ToBeInserted = numS1;
		numS2ToBeInserted = numS2;
	
		//std::cout << " mass " << massS1 << " " << massS2 << std::endl;
	
		tc = 0.005;
	
		speciesHandler.clear();
 
	    auto speciesDrum = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());	
		auto speciesS1 =  speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
		auto speciesS2 =  speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());

		//double CORInitial = CORDrum;
		double RPSInitial = 0.0;
		double particleWallFrictionInitial = particleWallFriction;
		double particleParticleFrictionInitial = particleParticleFriction;
		//		
		speciesDrum->setDensity(rhoS1);	
		speciesDrum->setCollisionTimeAndRestitutionCoefficient(tc, CORDrum, massS1);

		speciesDrum->setSlidingDissipation(speciesDrum->getDissipation()*2./7.);
		speciesDrum->setSlidingStiffness(speciesDrum->getStiffness()*2./7.);
		speciesDrum->setSlidingFrictionCoefficient(particleWallFrictionInitial);
		
		speciesDrum->setRollingStiffness(speciesDrum->getStiffness()*2.0/7.0);
		speciesDrum->setRollingFrictionCoefficient(0.5);
		speciesDrum->setRollingDissipation(speciesDrum->getDissipation()*2./7.);

		speciesDrum->setTorsionStiffness(speciesDrum->getStiffness()*2.0/7.0);
		speciesDrum->setTorsionFrictionCoefficient(0.01);
		speciesDrum->setTorsionDissipation(speciesDrum->getDissipation()*2./7.);
		//

		//
		speciesS1->setDensity(rhoS1);	
		speciesS1->setCollisionTimeAndRestitutionCoefficient(tc, CORS1, massS1);

		speciesS1->setSlidingDissipation(speciesS1->getDissipation()*2./7.);
		speciesS1->setSlidingStiffness(speciesS1->getStiffness()*2./7.);
		speciesS1->setSlidingFrictionCoefficient(particleParticleFrictionInitial);
		
		speciesS1->setRollingStiffness(speciesS1->getStiffness()*2.0/7.0);
		speciesS1->setRollingFrictionCoefficient(0.5);
		speciesS1->setRollingDissipation(speciesS1->getDissipation()*2./7.);

		speciesS1->setTorsionStiffness(speciesS1->getStiffness()*2.0/7.0);
		speciesS1->setTorsionFrictionCoefficient(0.01);
		speciesS1->setTorsionDissipation(speciesS1->getDissipation()*2./7.);
		//

		speciesS2->setDensity(rhoS2);	
		speciesS2->setCollisionTimeAndRestitutionCoefficient(tc, CORS2, massS2);

		speciesS2->setSlidingDissipation(speciesS2->getDissipation()*2./7.);
		speciesS2->setSlidingStiffness(speciesS2->getStiffness()*2./7.);
		speciesS2->setSlidingFrictionCoefficient(particleParticleFrictionInitial);
		
		speciesS2->setRollingStiffness(speciesS2->getStiffness()*2.0/7.0);
		speciesS2->setRollingFrictionCoefficient(0.5);
		speciesS2->setRollingDissipation(speciesS2->getDissipation()*2./7.);

		speciesS2->setTorsionStiffness(speciesS2->getStiffness()*2.0/7.0);
		speciesS2->setTorsionFrictionCoefficient(0.01);
		speciesS2->setTorsionDissipation(speciesS2->getDissipation()*2./7.);
	
		auto speciesDrumAndS1 = speciesHandler.getMixedObject(speciesDrum,speciesS1);

		speciesDrumAndS1->setCollisionTimeAndRestitutionCoefficient(tc, CORDrum, massS1, massS1);

		speciesDrumAndS1->setSlidingDissipation(speciesDrumAndS1->getDissipation()*2./7.);
		speciesDrumAndS1->setSlidingFrictionCoefficient(particleWallFrictionInitial);
		speciesDrumAndS1->setSlidingStiffness(speciesDrumAndS1->getStiffness()*2.0/7.0);

		speciesDrumAndS1->setRollingStiffness(speciesDrumAndS1->getStiffness()*2.0/7.0);
		speciesDrumAndS1->setRollingFrictionCoefficient(0.5);
		speciesDrumAndS1->setRollingDissipation(speciesDrumAndS1->getDissipation()*2./7.);

		speciesDrumAndS1->setTorsionStiffness(speciesDrumAndS1->getStiffness()*2.0/7.0);
		speciesDrumAndS1->setTorsionFrictionCoefficient(0.01);
		speciesDrumAndS1->setTorsionDissipation(speciesDrumAndS1->getDissipation()*2./7.);
		//
		auto speciesDrumAndS2 = speciesHandler.getMixedObject(speciesDrum,speciesS2);

		speciesDrumAndS2->setCollisionTimeAndRestitutionCoefficient(tc, CORDrum, massS1, massS2);

		speciesDrumAndS2->setSlidingDissipation(speciesDrumAndS2->getDissipation()*2./7.);
		speciesDrumAndS2->setSlidingFrictionCoefficient(particleWallFrictionInitial);
		speciesDrumAndS2->setSlidingStiffness(speciesDrumAndS2->getStiffness()*2.0/7.0);

		speciesDrumAndS2->setRollingStiffness(speciesDrumAndS2->getStiffness()*2.0/7.0);
		speciesDrumAndS2->setRollingFrictionCoefficient(0.5);
		speciesDrumAndS2->setRollingDissipation(speciesDrumAndS2->getDissipation()*2./7.);

		speciesDrumAndS2->setTorsionStiffness(speciesDrumAndS2->getStiffness()*2.0/7.0);
		speciesDrumAndS2->setTorsionFrictionCoefficient(0.01);
		speciesDrumAndS2->setTorsionDissipation(speciesDrumAndS2->getDissipation()*2./7.);
		//
		auto speciesS1AndS2 = speciesHandler.getMixedObject(speciesS1,speciesS2);

		speciesS1AndS2->setCollisionTimeAndRestitutionCoefficient(tc, CORS1, massS1, massS2);

		speciesS1AndS2->setSlidingDissipation(speciesS1AndS2->getDissipation()*2./7.);
		speciesS1AndS2->setSlidingFrictionCoefficient(particleParticleFrictionInitial);
		speciesS1AndS2->setSlidingStiffness(speciesS1AndS2->getStiffness()*2.0/7.0);

		speciesS1AndS2->setRollingStiffness(speciesS1AndS2->getStiffness()*2.0/7.0);
		speciesS1AndS2->setRollingFrictionCoefficient(.5);
		speciesS1AndS2->setRollingDissipation(speciesS1AndS2->getDissipation()*2./7.);

		speciesS1AndS2->setTorsionStiffness(speciesS1AndS2->getStiffness()*2.0/7.0);
		speciesS1AndS2->setTorsionFrictionCoefficient(0.01);
		speciesS1AndS2->setTorsionDissipation(speciesS1AndS2->getDissipation()*2./7.);

		Vec3D drumCenter = {0.5*(getXMin() + getXMax()),
							0.5*(getYMin() + getYMax()),
							0.5*(getZMin() + getZMax())};
		
		wallHandler.clear();
			
		auto drumWall = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
		drumWall -> setSpecies(speciesDrum);
		drumWall -> setPosition(drumCenter);
		drumWall -> setOrientation(Vec3D(0.0,1.0,0.0));
		drumWall -> addObject(Vec3D(1,0,0), Vec3D(drumRadius,0.0,0.0));
		drumWall -> setAngularVelocity(Vec3D(0.0,RPSInitial * 2.0 * constants::pi,0.0));
	
			
		InfiniteWall w0;
		w0.setSpecies(speciesDrum);

		w0.set(Vec3D(0.,-1.,0.),Vec3D(drumCenter.X,getYMin(),drumCenter.Z));
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D(0.,1.,0.),Vec3D(drumCenter.X,getYMax(),drumCenter.Z));
		wallHandler.copyAndAddObject(w0);


		BaseParticle P0;
		double radius = 0.0;
		int numS1Inserted=0;
		int numS2Inserted=0;
		Vec3D pos;
		double r, theta, y;
		int failCounter = 0;
		
		while( (numS1Inserted < numS1) || (numS2Inserted < numS2) )
        {
			int grn = random.getRandomNumber(1,numS1ToBeInserted+numS2ToBeInserted);

            if( grn > numS2ToBeInserted)
            {
            	radius = random.getRandomNumber((1. - fractionalPolydispersity)*radiusS1,(1. + fractionalPolydispersity)*radiusS1);//radius_1;
                P0.setSpecies(speciesS1);
                P0.setRadius(radius);
    
                failCounter = 0;
                do
                {
                    r = random.getRandomNumber(2.0*radius,drumRadius - 2.0*radius);
                    theta = random.getRandomNumber(0,constants::pi);
                    y = random.getRandomNumber(getYMin()+2.0*radius,getYMax()-2.0*radius);

                    pos.X = drumRadius + r*cos(theta);
                    pos.Y = y;
                    pos.Z = drumRadius + r*sin(theta);

                    P0.setPosition(pos);
                    P0.setVelocity(Vec3D(0.0,0.0,0.0));

                    failCounter++;
                    if (failCounter==1000) break;

                } while (checkParticleForInteraction(P0));
                
				numS1ToBeInserted--;
				numS1Inserted++;
            }
            else
            {
                radius = random.getRandomNumber((1. - fractionalPolydispersity)*radiusS2,(1. + fractionalPolydispersity)*radiusS2);//radius_2;
                P0.setSpecies(speciesS2);
                P0.setRadius(radius);

                failCounter = 0;
                do
                {
                    r = random.getRandomNumber(2.0*radius,drumRadius - 2.0*radius);
                    theta = random.getRandomNumber(constants::pi,constants::pi*2.);
                    y = random.getRandomNumber(getYMin()+2.0*radius,getYMax()-2.0*radius);

                    pos.X = drumRadius + r*cos(theta);
                    pos.Y = y;
                    pos.Z = drumRadius + r*sin(theta);

                    P0.setPosition(pos);
                    P0.setVelocity(Vec3D(0.0,0.0,0.0));

                    failCounter++;
                    if (failCounter==1000) break;

                } while (checkParticleForInteraction(P0));


				numS2ToBeInserted--;
				numS2Inserted++;
            }
            // For homogeneous mix
            /*
			failCounter = 0;			
          	do
          	{
				//r = random.getRandomNumber(drumCenter.X-10.0*radius,drumCenter.X+10.0*radius);
				//theta = random.getRandomNumber(0,constants::pi*2.);
				//y = random.getRandomNumber(getYMin()+2.0*radius,getYMax()-2.0*radius);

                r = random.getRandomNumber(2.0*radius,drumRadius-2.0*radius);
				theta = random.getRandomNumber(0,constants::pi*2.);
				y = random.getRandomNumber(getYMin()+2.0*radius,getYMax()-2.0*radius);

				pos.X = drumRadius + r*cos(theta);
				pos.Y = y;
				pos.Z = drumRadius + r*sin(theta);
          
				P0.setPosition(pos);
				P0.setVelocity(Vec3D(0.0,0.0,0.0));
	
				failCounter++;
				if (failCounter==1000) break;
	
			} while (checkParticleForInteraction(P0));
            */
			particleHandler.copyAndAddObject(P0);
            
            hGridRebuild();
        }

		std::cout << "Finished creating particles" << std::endl;
		std::cout << "Number of S1 particles inserted" << numS1Inserted << std::endl;
		std::cout << "Number of S2 particles inserted" << numS2Inserted << std::endl;

		//hGridRebuild();

		if ( (numS1ToBeInserted==0) && (numS2ToBeInserted==0) )
		{
			step = 2;
			std::cout << "\n \n \n";
			std::cout << "Particles settling down" << std::endl;
			std::cout << "--------------------------" << std::endl;
			std::cout << "\n\n\n";
			checkTime = getTime() + 5.0;
		}		
	}

	void actionsBeforeTimeStep()
	{
		wallHandler.getObject(0)->setOrientation(Vec3D(0.0,1.0,0.0));
		wallHandler.getObject(1)->setOrientation(Vec3D(0.0,1.0,0.0));
		wallHandler.getObject(2)->setOrientation(Vec3D(0.0,1.0,0.0));

		if (step==2)
		{
			if (getTime() > checkTime)
			{
				std::cout << "Current KE" << getKineticEnergy() << std::endl;	
				if (getKineticEnergy() < (10.0))
				{
					step = 3;
					double drumStartTime = getTime();
					std::cout << "\n \n \n";
					std::cout << "Start the drum rotation" << std::endl;
					std::cout << "--------------------------" << std::endl;
					std::cout << "\n\n\n";
					// rotate the drum
					wallHandler.getObject(0)->setAngularVelocity(Vec3D(0.0,revolutionsPerSecond*constants::pi*2.0,0.0));
					wallHandler.getObject(1)->setAngularVelocity(Vec3D(0.0,revolutionsPerSecond*constants::pi*2.0,0.0));
					wallHandler.getObject(2)->setAngularVelocity(Vec3D(0.0,revolutionsPerSecond*constants::pi*2.0,0.0));

					wallHandler.getObject(0)->setOrientation(Vec3D(0.0,1.0,0.0));
					wallHandler.getObject(1)->setOrientation(Vec3D(0.0,1.0,0.0));
					wallHandler.getObject(2)->setOrientation(Vec3D(0.0,1.0,0.0));

					/*
					for (int i = 0; i < particleHandler.getNumberOfObjects();i++)
					{
						BaseParticle* P0 = particleHandler.getObject(i);
						if (P0->getIndSpecies() == 1)
		                {
							P0->setSpecies(speciesS1);
							P0->setForce(Vec3D(0,0,0));
							P0->setTorque(Vec3D(0,0,0));
                 		}
						else
						{
							P0->setSpecies(speciesS2);
							P0->setForce(Vec3D(0,0,0));
							P0->setTorque(Vec3D(0,0,0));
						}
					}
					*/	
				}
				else
				{
					checkTime = getTime() + 1.0;
				}
			}
		}
	}

	void setDrumRadius (double radius)
	{
		drumRadius = radius;
	}
	
	void setRevolutionSpeed (double rpm)
	{
		revolutionsPerSecond = rpm/60.0;// non-dimensionalised based on 3 mm particles and g=9.81
	}

	void setSizeAndDensityRatio (double sr, double dr)
	{
		sizeRatio = sr;
		densityRatio = dr;
	}

	void setFractionalPolydispersity(double fpd)
	{
		fractionalPolydispersity = fpd;
	}

	void setDrumFillFraction (double dff)
	{
		drumFillFraction = dff;
	}

	void setSpeciesVolumeFraction(double vf)
	{
		volumeFraction = vf;
	}

	void setFrictionCoeff(double pwf, double ppf)
	{
		particleWallFriction = pwf;
		particleParticleFriction = ppf;
	}

	void setCOR (double drumCOR, double COR1, double COR2)
	{
		CORDrum = drumCOR;
		CORS1 = COR1;
		CORS2 = COR2;
	}

	double getDrumRadius()
	{
		return drumRadius;
	}

    double getSizeRatio ()
	{
		return sizeRatio;
	}
    
    double getLargestParticleRadius()
    {
        if (radiusS1 >= radiusS2)
        {
            return radiusS1;
        }
        else
        {
            return radiusS2;
        }
    }

	private:

	double radiusS1,radiusS2;
	double rhoS1, rhoS2;
	double massS1, massS2;

	double CORDrum,CORS1,CORS2,tc;

	double sizeRatio;	
	double densityRatio;
	double drumFillFraction;
	double volumeFraction;

	double particleWallFriction,particleParticleFriction;

	int numS1,numS2;
	int numS1ToBeInserted,numS2ToBeInserted;
	
	double drumRadius;
	double revolutionsPerSecond;
	double fractionalPolydispersity;

	int step;
	double checkTime;
};

int main(int argc, char *argv[])
{

	RotatingDrum problem;

	problem.setName("Bidisperse_2PolyD_05");
	//problem.autoNumber();
	//problem.setDrumRadius(0.125);// in meters
        problem.setDrumRadius(0.125);
        
	problem.setXMin(0.0);
	problem.setYMin(0.0);
	problem.setZMin(0.0);

	problem.setXMax(2. * problem.getDrumRadius());
	problem.setZMax(2. * problem.getDrumRadius());
	problem.setYMax(0.03);// in meters
	
	problem.setTimeMax(100.);
	problem.setTimeStep(0.005/50.);

	problem.setGravity(Vec3D(0.,0.,-9.81));
	problem.setCOR(0.97,0.831,0.831);//drumWall, species1, species2
	problem.setSizeAndDensityRatio(1.5,1.0);//size ratio and density ratio
	problem.setFractionalPolydispersity(0.02);//10% dispersity
	problem.setDrumFillFraction(0.25);// At 0.5 the drum is 3/4 filled. 
	problem.setSpeciesVolumeFraction(.5);//Species1 volume fraction
	problem.setFrictionCoeff(.6,.19);//(particle-wall-friction,part-part-friction)
	problem.setRevolutionSpeed(3.);//rpm

	problem.setSaveCount(10000);
	problem.setXBallsAdditionalArguments("-cmode 8 -solidf -v0");
	problem.readArguments(argc,argv);
        
        problem.dataFile.setFileType(FileType::MULTIPLE_FILES);
        problem.restartFile.setFileType(FileType::ONE_FILE);
        problem.fStatFile.setFileType(FileType::NO_FILE);
        problem.eneFile.setFileType(FileType::NO_FILE);

        problem.setWallsWriteVTK(FileType::MULTIPLE_FILES);
        problem.setParticlesWriteVTK(true);

	problem.solve();

	return 0;
}
