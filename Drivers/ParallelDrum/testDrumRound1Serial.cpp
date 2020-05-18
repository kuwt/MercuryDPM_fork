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
#include "Mercury3D.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"
#include <chrono>

class RotatingDrum : public Mercury3D
{
	public:
		
	RotatingDrum()
	{	
		radiusS1 = 0.0015; // 3mm diameter
		fractionalPolydispersity = 0.0;
	}

	void setupInitialConditions() override {

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
		
		tc = 1 / 800.0;
		//original value
		//tc = 0.005;
	
		speciesHandler.clear();
 
		auto speciesDrum = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());	
		auto speciesS1 =  speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
		auto speciesS2 =  speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());

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

		speciesS2->setDensity(rhoS2);	
		speciesS2->setCollisionTimeAndRestitutionCoefficient(tc, CORS2, massS2);

		speciesS2->setSlidingDissipation(speciesS2->getDissipation()*2./7.);
		speciesS2->setSlidingStiffness(speciesS2->getStiffness()*2./7.);
		speciesS2->setSlidingFrictionCoefficient(slidingFriction2);
		
		speciesS2->setRollingStiffness(speciesS2->getStiffness()*2.0/7.0);
		speciesS2->setRollingFrictionCoefficient(rollingFriction2);
		speciesS2->setRollingDissipation(speciesS2->getDissipation()*2./7.);

		speciesS2->setTorsionStiffness(speciesS2->getStiffness()*2.0/7.0);
		speciesS2->setTorsionFrictionCoefficient(torsionFriction2);
		speciesS2->setTorsionDissipation(speciesS2->getDissipation()*2./7.);
	
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
		auto speciesDrumAndS2 = speciesHandler.getMixedObject(speciesDrum,speciesS2);

		speciesDrumAndS2->setCollisionTimeAndRestitutionCoefficient(tc, ((CORDrum + CORS2) / 2), massS1, massS2);

		speciesDrumAndS2->setSlidingDissipation(speciesDrumAndS2->getDissipation()*2./7.);
		speciesDrumAndS2->setSlidingFrictionCoefficient(((slidingFrictionDrum + slidingFriction2)/2));
		speciesDrumAndS2->setSlidingStiffness(speciesDrumAndS2->getStiffness()*2.0/7.0);

		speciesDrumAndS2->setRollingStiffness(speciesDrumAndS2->getStiffness()*2.0/7.0);
		speciesDrumAndS2->setRollingFrictionCoefficient(((rollingFrictionDrum + rollingFriction2)/2));
		speciesDrumAndS2->setRollingDissipation(speciesDrumAndS2->getDissipation()*2./7.);

		speciesDrumAndS2->setTorsionStiffness(speciesDrumAndS2->getStiffness()*2.0/7.0);
		speciesDrumAndS2->setTorsionFrictionCoefficient(((torsionFrictionDrum + torsionFriction2)/2));
		speciesDrumAndS2->setTorsionDissipation(speciesDrumAndS2->getDissipation()*2./7.);
		//
		auto speciesS1AndS2 = speciesHandler.getMixedObject(speciesS1,speciesS2);

		speciesS1AndS2->setCollisionTimeAndRestitutionCoefficient(tc, ((CORS1 + CORS2) / 2), massS1, massS2);

		speciesS1AndS2->setSlidingDissipation(speciesS1AndS2->getDissipation()*2./7.);
		speciesS1AndS2->setSlidingFrictionCoefficient(((rollingFriction1 + rollingFriction2)/2));
		speciesS1AndS2->setSlidingStiffness(speciesS1AndS2->getStiffness()*2.0/7.0);

		speciesS1AndS2->setRollingStiffness(speciesS1AndS2->getStiffness()*2.0/7.0);
		speciesS1AndS2->setRollingFrictionCoefficient(((rollingFriction1 + rollingFriction2)/2));
		speciesS1AndS2->setRollingDissipation(speciesS1AndS2->getDissipation()*2./7.);

		speciesS1AndS2->setTorsionStiffness(speciesS1AndS2->getStiffness()*2.0/7.0);
		speciesS1AndS2->setTorsionFrictionCoefficient(((torsionFriction1 + torsionFriction2)/2));
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
	
		drumWall->setPrescribedPosition([this] (double time)
		{
			double t = time - 0;
			if (t > 0.0)
			{
				//return drumCenter + Vec3D(0.0,0.0,vibrationAmp * std::sin(t * 2.0 * vibrationFreq * constants::pi));
				return Vec3D (0.5*(getXMin() + getXMax()),
								0.5*(getYMin() + getYMax()),
								0.5*(getZMin() + getZMax()) + vibrationAmp * std::sin(t * 2.0 * vibrationFreq * constants::pi)
								);
			}
		else
			{
				//return drumCenter;
				return Vec3D (0.5*(getXMin() + getXMax()),
								0.5*(getYMin() + getYMax()),
								0.5*(getZMin() + getZMax()));
			}
		});
			
		InfiniteWall w0;
		w0.setSpecies(speciesDrum);

		w0.set(Vec3D(0.,-1.,0.),Vec3D(drumCenter.X,getYMin(),drumCenter.Z));
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D(0.,1.,0.),Vec3D(drumCenter.X,getYMax(),drumCenter.Z));
		wallHandler.copyAndAddObject(w0);


		SphericalParticle P0;
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

	void actionsBeforeTimeStep() override {

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
	//a series of functions by which to easily set particles'
	//various frictional coefficients	
	void setSlidingFriction (double drum, double f1, double f2)
	{
		slidingFrictionDrum = drum;		
		slidingFriction1 = f1;
		slidingFriction2 = f2;	
	}

	void setRollingFriction (double drum, double f1, double f2)
	{
		rollingFrictionDrum = drum;		
		rollingFriction1 = f1;
		rollingFriction2 = f2;	
	}

	void setTorsionFriction (double drum, double f1, double f2)
	{
		torsionFrictionDrum = drum;		
		torsionFriction1 = f1;
		torsionFriction2 = f2;	
	}

	double getDrumRadius()
	{
		return drumRadius;
	}
	
	//new functions to set the frequency and amplitude with which the drum is
	//vibrated
	void setVibrationFrequency (double f) {
		vibrationFreq = f;
	}
	void setVibrationAmplitude (double A) {
		vibrationAmp = A;
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
	
	double slidingFrictionDrum;		
	double slidingFriction1;
	double slidingFriction2;	
	
	double rollingFrictionDrum;		
	double rollingFriction1;
	double rollingFriction2;	

	double torsionFrictionDrum;		
	double torsionFriction1;
	double torsionFriction2;	

	int step;
	double checkTime;
	
	//New parameters to give the strength of applied vibrations
	double vibrationAmp;
	double vibrationFreq;
};

int main(int argc, char *argv[])
{

    // Start measuring elapsed time
    std::chrono::time_point<std::chrono::system_clock> startClock, endClock;
    startClock = std::chrono::system_clock::now();

	RotatingDrum problem;
	
	//setting locally-available variables to define key
	//parameters such that they can be automatically included in
	//the name of the file
	//*******************Excitation Properties*****************************************
	//Drum rotation rate (rpm)
	double rotRateBasal = 10.0;	
	
	//vibration amplitude and frequency
	double fVib = 0.0;
	double aVib = 0.0;
	
	//*******************Frictional Properties*****************************************
	//sliding friction for species 1, 2 and wall
	double muSWall = 0.6;
	double muS1 = 0.1;
	double muS2 = 0.1;
	
	//rolling friction for species 1, 2 and wall
	double muRWall = 0.06;
	double muR1 = 0.01;
	double muR2 = 0.01;
	
	//torsion friction for species 1, 2 and wall
	double muTWall = 0.0;
	double muT1 = 0.000;
	double muT2 = 0.000;
	
	//the density ratio of particles
	double rhoRatio = 2500.0/2500.0;
	//the size ratio of particles
	double dRatio = 1.5;
	
	//the fraction of species 1 particles
	double specFrac = 0.5;
	
	double drumRad = 20*1.5*0.0015;

	double rotRate = 10.0;
		
	double froudeNumber = (rotRate * 2.0 * 3.1415926535 / 60.0) * (rotRate * 2.0 * 3.1415926535 / 60.0) * drumRad / 9.81;

	//the dimensionless drum length (L/d_l)
	double dimDrumLength = 240.0;

	
	double drumLength = dimDrumLength*0.003*1.5;

	//Set the number of domains for parallel decomposition
	problem.setNumberOfDomains({2,1,1});

	  
	//*******************Setting Up Filename******************************************
	//setting up a stringstream to use in order to create an instructive filename
	std::stringstream nameStream;
	//the generic 'root' of the name for the file that will be applied to all files
	//irrespective of parameters
	//std::string nameBase = "binary";
	std::string nameBase = "binary";
	
	//Name stream for tests looking at size segregation and concentration
	nameStream << nameBase << "-length" << dimDrumLength;
       

	//Name stream for tests looking at effect of vibration
	/*nameStream << nameBase << "-rotationRate" << rotRate 
				<< "-frequency" << fVib << "-amplitude" << aVib;
	*/

	/*Name stream for tests looking for particle 'sinkage'
	 * nameStream << nameBase << "-mu_s" << muS1 << "," << muS2 << "," << muSWall << ","
							<< "-mu_r" << muR1 << "," << muR2 << "," << muRWall << ","
							<< "-mu_t" << muT1 << "," << muT2 << "," << muTWall
							<< "-velocity" << rotRate << "-densityRatio" << rhoRatio;
	*/
	std::string fullName = nameStream.str();
	problem.setName(fullName);
	//problem.autoNumber();
	problem.setDrumRadius(drumRad);// in meters

	problem.setXMin(0.0);
	problem.setYMin(0.0);
	problem.setZMin(0.0);

	problem.setXMax(2. * problem.getDrumRadius());
	problem.setZMax(2. * problem.getDrumRadius());
	problem.setYMax(drumLength);// in meters
	
	problem.setTimeMax(750.);
	problem.setTimeStep(1.0/(800.0 * 50.0));

	problem.setGravity(Vec3D(0.,0.,-9.81));
	problem.setCOR(0.97,0.97,0.97);//drumWall, species1, species2
	problem.setSizeAndDensityRatio(dRatio,rhoRatio);//size ratio and density ratio
	problem.setFractionalPolydispersity(0.05);//10% dispersity
	problem.setDrumFillFraction(0.3);// At 0.5 the drum is 3/4 filled. 
	problem.setSpeciesVolumeFraction(specFrac);//Species1 volume fraction
	//redundant
	//problem.setFrictionCoeff(.6,.19);//(particle-wall-friction,part-part-friction)
	problem.setSlidingFriction(muSWall,muS1,muS2); //wall, species1, species2
	problem.setRollingFriction(muRWall,muR1,muR2); //wall, species1, species2
	problem.setTorsionFriction(muTWall,muT1,muT2); //wall, species1, species2

	problem.setRevolutionSpeed(rotRate);//rpm
	
	//setting vibration parameters
	problem.setVibrationAmplitude(aVib);
	problem.setVibrationFrequency(fVib);

	problem.setSaveCount(20000);
	problem.setXBallsAdditionalArguments("-cmode 8 -solidf -v0");
	problem.readArguments(argc,argv);

	problem.solve();

    // Measure elapsed time
    endClock = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = endClock - startClock;
    logger(INFO, "Elapsed time for solving the PDE: % s", elapsed_seconds.count());

	return 0;
}
