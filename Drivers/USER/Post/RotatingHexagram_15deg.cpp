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

// Tutorial 7

/*
** This file is annotated with DoxyFile comments in order to show the code on
** the documentation - This is not needed for your real drivers.
** Please ignore these comments.
*/

//! [T7:headers]
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Mercury3D.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"
//! [T7:headers]

//! [T7:class]
class RotatingHexagon : public Mercury3D
{
  public:

      	RotatingHexagon()
	{	
		
		fractionalPolydispersity = 0.0;
	}

// SETUP INITIAL CONDITIONS        
  void setupInitialConditions()
  {
// INITIAL CALCULATIONS
                radS1 = 0.003; // 3mm radius
                radS2 = sizeRatio*radS1;
		
		rhoS1 = 2500.0;
		rhoS2 = densityRatio*rhoS1;

		massS1 = 4 / 3 * constants::pi * pow(radS1, 3.0) * rhoS1;
		massS2 = 4 / 3 * constants::pi * pow(radS2, 3.0) * rhoS2;

 //		double fillVolume = drumFillFraction*(6*pow(getZMax(),2.)*tan(constants::pi/6.)*(getYMax()-getYMin()));
                double fillVolume = drumFillFraction*(getYMax()-getYMin())*(6*pow(getZMax(),2.)*tan(constants::pi/6.)- 6*pow(tan(constants::pi/6)*getZMax(),2));

		numS1 = volumeFraction*fillVolume/(4./3. * constants::pi*pow(radS1,3.0));
		numS2 = (1. - volumeFraction)*fillVolume/(4./3. * constants::pi*pow(radS2,3.0));
		
		numS1ToBeInserted = numS1;
		numS2ToBeInserted = numS2;
	
		tc = 0.005;
                rad_excl = getZMax()*(1./cos(constants::pi/6));
// INITIAL CALCULATIONS *
                
// SETTING THE SPECIES
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
		speciesDrum->setRollingFrictionCoefficient(0.2);
		speciesDrum->setRollingDissipation(speciesDrum->getDissipation()*2./7.);

		speciesDrum->setTorsionStiffness(speciesDrum->getStiffness()*2.0/7.0);
		speciesDrum->setTorsionFrictionCoefficient(0.1);
		speciesDrum->setTorsionDissipation(speciesDrum->getDissipation()*2./7.);
		//

		//
		speciesS1->setDensity(rhoS1);	
		speciesS1->setCollisionTimeAndRestitutionCoefficient(tc, CORS1, massS1);

		speciesS1->setSlidingDissipation(speciesS1->getDissipation()*2./7.);
		speciesS1->setSlidingStiffness(speciesS1->getStiffness()*2./7.);
		speciesS1->setSlidingFrictionCoefficient(particleParticleFrictionInitial);
		
		speciesS1->setRollingStiffness(speciesS1->getStiffness()*2.0/7.0);
		speciesS1->setRollingFrictionCoefficient(0.2);
		speciesS1->setRollingDissipation(speciesS1->getDissipation()*2./7.);

		speciesS1->setTorsionStiffness(speciesS1->getStiffness()*2.0/7.0);
		speciesS1->setTorsionFrictionCoefficient(0.1);
		speciesS1->setTorsionDissipation(speciesS1->getDissipation()*2./7.);
		//

		speciesS2->setDensity(rhoS2);	
		speciesS2->setCollisionTimeAndRestitutionCoefficient(tc, CORS2, massS2);

		speciesS2->setSlidingDissipation(speciesS2->getDissipation()*2./7.);
		speciesS2->setSlidingStiffness(speciesS2->getStiffness()*2./7.);
		speciesS2->setSlidingFrictionCoefficient(particleParticleFrictionInitial);
		
		speciesS2->setRollingStiffness(speciesS2->getStiffness()*2.0/7.0);
		speciesS2->setRollingFrictionCoefficient(0.2);
		speciesS2->setRollingDissipation(speciesS2->getDissipation()*2./7.);

		speciesS2->setTorsionStiffness(speciesS2->getStiffness()*2.0/7.0);
		speciesS2->setTorsionFrictionCoefficient(0.1);
		speciesS2->setTorsionDissipation(speciesS2->getDissipation()*2./7.);
	
		auto speciesDrumAndS1 = speciesHandler.getMixedObject(speciesDrum,speciesS1);

		speciesDrumAndS1->setCollisionTimeAndRestitutionCoefficient(tc, CORDrum, massS1, massS1);

		speciesDrumAndS1->setSlidingDissipation(speciesDrumAndS1->getDissipation()*2./7.);
		speciesDrumAndS1->setSlidingFrictionCoefficient(particleWallFrictionInitial);
		speciesDrumAndS1->setSlidingStiffness(speciesDrumAndS1->getStiffness()*2.0/7.0);

		speciesDrumAndS1->setRollingStiffness(speciesDrumAndS1->getStiffness()*2.0/7.0);
		speciesDrumAndS1->setRollingFrictionCoefficient(0.2);
		speciesDrumAndS1->setRollingDissipation(speciesDrumAndS1->getDissipation()*2./7.);

		speciesDrumAndS1->setTorsionStiffness(speciesDrumAndS1->getStiffness()*2.0/7.0);
		speciesDrumAndS1->setTorsionFrictionCoefficient(0.1);
		speciesDrumAndS1->setTorsionDissipation(speciesDrumAndS1->getDissipation()*2./7.);
		//
		auto speciesDrumAndS2 = speciesHandler.getMixedObject(speciesDrum,speciesS2);

		speciesDrumAndS2->setCollisionTimeAndRestitutionCoefficient(tc, CORDrum, massS1, massS2);

		speciesDrumAndS2->setSlidingDissipation(speciesDrumAndS2->getDissipation()*2./7.);
		speciesDrumAndS2->setSlidingFrictionCoefficient(particleWallFrictionInitial);
		speciesDrumAndS2->setSlidingStiffness(speciesDrumAndS2->getStiffness()*2.0/7.0);

		speciesDrumAndS2->setRollingStiffness(speciesDrumAndS2->getStiffness()*2.0/7.0);
		speciesDrumAndS2->setRollingFrictionCoefficient(0.2);
		speciesDrumAndS2->setRollingDissipation(speciesDrumAndS2->getDissipation()*2./7.);

		speciesDrumAndS2->setTorsionStiffness(speciesDrumAndS2->getStiffness()*2.0/7.0);
		speciesDrumAndS2->setTorsionFrictionCoefficient(0.1);
		speciesDrumAndS2->setTorsionDissipation(speciesDrumAndS2->getDissipation()*2./7.);
		//
		auto speciesS1AndS2 = speciesHandler.getMixedObject(speciesS1,speciesS2);

		speciesS1AndS2->setCollisionTimeAndRestitutionCoefficient(tc, CORS1, massS1, massS2);

		speciesS1AndS2->setSlidingDissipation(speciesS1AndS2->getDissipation()*2./7.);
		speciesS1AndS2->setSlidingFrictionCoefficient(particleParticleFrictionInitial);
		speciesS1AndS2->setSlidingStiffness(speciesS1AndS2->getStiffness()*2.0/7.0);

		speciesS1AndS2->setRollingStiffness(speciesS1AndS2->getStiffness()*2.0/7.0);
		speciesS1AndS2->setRollingFrictionCoefficient(.2);
		speciesS1AndS2->setRollingDissipation(speciesS1AndS2->getDissipation()*2./7.);

		speciesS1AndS2->setTorsionStiffness(speciesS1AndS2->getStiffness()*2.0/7.0);
		speciesS1AndS2->setTorsionFrictionCoefficient(0.1);
		speciesS1AndS2->setTorsionDissipation(speciesS1AndS2->getDissipation()*2./7.);
// SETTING THE SPECIES *
                
// PLACING THE WALLS
//              wallHandler.clear();
			
//		auto drumWall = wallHandler.copyAndAddObject(InfiniteWalls());
//		drumWall -> setSpecies(speciesDrum);
//		drumWall -> setPosition(drumCenter);
//		drumWall -> setOrientation(Vec3D(0.0,1.0,0.0));
//		drumWall -> addObject(Vec3D(1,0,0), Vec3D(drumRadius,0.0,0.0));
//		drumWall -> setAngularVelocity(Vec3D(0.0,RPSInitial * 2.0 * constants::pi,0.0));

                wallHandler.clear();
                InfiniteWall w0;
                w0.setSpecies(speciesDrum);
                w0.setAngularVelocity(Vec3D(0.0,RPSInitial * 2.0 * constants::pi,0.0));
                
//                w0.set(Vec3D(1.0,0.0,-(rad_excl/getZMax() -tan(constants::pi/6))),Vec3D((rad_excl/getZMax())*getZMax(),0.0,0.0));
//                wallHandler.copyAndAddObject(w0);
//                w0.set(Vec3D(1.0,0.0,((rad_excl/getZMax())-tan(constants::pi/6))),Vec3D((rad_excl/getZMax())*getZMax(),0.0,0.0));
//                wallHandler.copyAndAddObject(w0);
//                w0.set(Vec3D(-1.0,0.0,((rad_excl/getZMax())-tan(constants::pi/6))),Vec3D((rad_excl/getZMax())*getZMin(),0.0,0.0));
//                wallHandler.copyAndAddObject(w0);                
//                w0.set(Vec3D(-1.0,0.0,-((rad_excl/getZMax())-tan(constants::pi/6))),Vec3D((rad_excl/getZMax())*getZMin(),0.0,0.0));
//                wallHandler.copyAndAddObject(w0);   
                
                //Upper and lower side of the hexagon, these are horizontal
                w0.set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,getZMax()));
                wallHandler.copyAndAddObject(w0);
                w0.set(Vec3D(0.0,0.0,-1.0),Vec3D(0.0,0.0,getZMin()));
                wallHandler.copyAndAddObject(w0);
                //Front and back wall
                w0.set(Vec3D(0.0,1.0,0.0),Vec3D(0.0,getYMax(),0.0));
                wallHandler.copyAndAddObject(w0);
                w0.set(Vec3D(0.0,-1.0,0.0),Vec3D(0.0,getYMin(),0.0));
                wallHandler.copyAndAddObject(w0);               
                
                //Points
                IntersectionOfWalls w1;
                w1.setSpecies(speciesDrum);

                std::vector<Vec3D> Points(3);
                //Top Triangle
                Points[0] = Vec3D(tan(constants::pi /6.)*getZMax(),0.0,getZMax());        
                Points[1] = Vec3D(0.0,0.0,(1-tan(constants::pi /12.))*getZMax());               
                Points[2] = Vec3D(-tan(constants::pi /6.)*getZMax(),0.0,getZMax());
                w1.createOpenPrism(Points);
                wallHandler.copyAndAddObject(w1);
                //Bottom Triangle
                Points[0] = Vec3D(tan(constants::pi /6.)*getZMin(),0.0,getZMin());        
                Points[1] = Vec3D(0.0,0.0,(1-tan(constants::pi /12.))*getZMin());               
                Points[2] = Vec3D(-tan(constants::pi /6.)*getZMin(),0.0,getZMin());
                w1.createOpenPrism(Points);
                wallHandler.copyAndAddObject(w1);
                //Right Top Triangle
                Points[0] = Vec3D(tan(constants::pi /6.)*getZMax(),0.0,getZMax());                  
                Points[1] = Vec3D(cos(constants::pi /6)*15./135.*(2.*getZMax()/tan(constants::pi /12)),0.0,sin(constants::pi /6)*15./135.*(2.*getZMax()/tan(constants::pi /12)));               
                Points[2] = Vec3D(getZMax()/cos(constants::pi /6),0.0,0.0);        
                w1.createOpenPrism(Points);
                wallHandler.copyAndAddObject(w1);                
                //Left Top Triangle
                Points[0] = Vec3D(tan(constants::pi /6.)*getZMin(),0.0,getZMax());                  
                Points[1] = Vec3D(-cos(constants::pi /6)*15./135.*(2.*getZMax()/tan(constants::pi /12)),0.0,sin(constants::pi /6)*15./135.*(2.*getZMax()/tan(constants::pi /12)));               
                Points[2] = Vec3D(getZMin()/cos(constants::pi /6),0.0,0.0);        
                w1.createOpenPrism(Points);
                wallHandler.copyAndAddObject(w1);    
                //Right Bottom Triangle
                Points[0] = Vec3D(tan(constants::pi /6.)*getZMax(),0.0,getZMin());                  
                Points[1] = Vec3D(cos(constants::pi /6)*15./135.*(2.*getZMax()/tan(constants::pi /12)),0.0,-sin(constants::pi /6)*15./135.*(2.*getZMax()/tan(constants::pi /12)));               
                Points[2] = Vec3D(getZMax()/cos(constants::pi /6),0.0,0.0);        
                w1.createOpenPrism(Points);
                wallHandler.copyAndAddObject(w1);                
                //Left Bottom Triangle
                Points[0] = Vec3D(tan(constants::pi /6.)*getZMin(),0.0,getZMin());                  
                Points[1] = Vec3D(-cos(constants::pi /6)*15./135.*(2.*getZMax()/tan(constants::pi /12)),0.0,-sin(constants::pi /6)*15./135.*(2.*getZMax()/tan(constants::pi /12)));               
                Points[2] = Vec3D(getZMin()/cos(constants::pi /6),0.0,0.0);        
                w1.createOpenPrism(Points);
                wallHandler.copyAndAddObject(w1);                  
                
                
// PLACING THE WALLS *

                
// PLACING THE PARTICLES
    // SETTING THE SPECIES FOR THE BASEPARTICLE
                SphericalParticle p0;

		int numS1Inserted=0;
		int numS2Inserted=0;
		Vec3D pos;
		double r, theta, y, x, z;
		int failCounter = 0;
    // SETTING THE SPECIES FOR THE BASEPARTICLE
                
    while( (numS1Inserted < numS1) || (numS2Inserted < numS2) )
        {
            int grn = random.getRandomNumber(1,numS1ToBeInserted+numS2ToBeInserted);

            if( grn > numS2ToBeInserted)
            {
            	//radius = random.getRandomNumber((1. - fractionalPolydispersity)*radS1,(1. + fractionalPolydispersity)*radS1);//radius_1;
                p0.setSpecies(speciesS1);
                p0.setRadius(radS1);
    
                failCounter = 0;
                do
                {
                    theta = random.getRandomNumber(0,constants::pi);
                    r = random.getRandomNumber(2.0*radS1,getZMax()-1.1*radS1);
                    y = random.getRandomNumber(getYMin()+2.0*radS1,getYMax()-2.0*radS1);

                    pos.X = r*cos(theta);
                    pos.Y = y;
                    pos.Z = r*sin(theta);

                    p0.setPosition(pos);
                    p0.setVelocity(Vec3D(0.0,0.0,0.0));

                    failCounter++;
                    if (failCounter==1000) break;
                } while (checkParticleForInteraction(p0));
				numS1ToBeInserted--;
				numS1Inserted++;
            }
            else
            {
                //radius = random.getRandomNumber((1. - fractionalPolydispersity)*radS2,(1. + fractionalPolydispersity)*radS2);//radius_2;
                p0.setSpecies(speciesS2);
                p0.setRadius(radS2);

                failCounter = 0;
                do
                {
//                    x = random.getRandomNumber(getZMin()+2.0*radS2,getZMax()-2.0*radS2);
//                    z = random.getRandomNumber(getZMin()+2.0*radS2,getZMax()-2.0*radS2);
                    
                    theta = random.getRandomNumber(constants::pi,constants::pi * 2.0);
                    r = random.getRandomNumber(2.0*radS2,getZMax()-1.1*radS2);
                    y = random.getRandomNumber(getYMin()+2.0*radS2,getYMax()-2.0*radS2);

                    pos.X = r*cos(theta);
                    pos.Y = y;
                    pos.Z = r*sin(theta);

                    p0.setPosition(pos);
                    p0.setVelocity(Vec3D(0.0,0.0,0.0));

                    failCounter++;
                    if (failCounter==1000) break;

                } while (checkParticleForInteraction(p0));
        			numS2ToBeInserted--;
				numS2Inserted++;
            }
 			if (failCounter==1000) break;           
			particleHandler.copyAndAddObject(p0);
            
            hGridRebuild();
        }

    std::cout << "Finished creating particles" << std::endl;
    std::cout << "Number of S1 particles inserted" << numS1Inserted << std::endl;
    std::cout << "Number of S2 particles inserted" << numS2Inserted << std::endl;

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
// PLACING THE PARTICLES *
           
// SETTING INITIAL CONDITIONS *
                


// ACTIONS BEFORE TIMESTEP
                
void actionsBeforeTimeStep()
	{
		wallHandler.getObject(0)->setOrientation(Vec3D(0.0,1.0,0.0));
		wallHandler.getObject(1)->setOrientation(Vec3D(0.0,1.0,0.0));
		wallHandler.getObject(2)->setOrientation(Vec3D(0.0,1.0,0.0));
                wallHandler.getObject(3)->setOrientation(Vec3D(0.0,1.0,0.0));
		wallHandler.getObject(4)->setOrientation(Vec3D(0.0,1.0,0.0));
		wallHandler.getObject(5)->setOrientation(Vec3D(0.0,1.0,0.0));
                wallHandler.getObject(6)->setOrientation(Vec3D(0.0,1.0,0.0));
		wallHandler.getObject(7)->setOrientation(Vec3D(0.0,1.0,0.0));


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
                                        wallHandler.getObject(3)->setAngularVelocity(Vec3D(0.0,revolutionsPerSecond*constants::pi*2.0,0.0));
					wallHandler.getObject(4)->setAngularVelocity(Vec3D(0.0,revolutionsPerSecond*constants::pi*2.0,0.0));
					wallHandler.getObject(5)->setAngularVelocity(Vec3D(0.0,revolutionsPerSecond*constants::pi*2.0,0.0));
                                        wallHandler.getObject(6)->setAngularVelocity(Vec3D(0.0,revolutionsPerSecond*constants::pi*2.0,0.0));
					wallHandler.getObject(7)->setAngularVelocity(Vec3D(0.0,revolutionsPerSecond*constants::pi*2.0,0.0));
                                        
					wallHandler.getObject(0)->setOrientation(Vec3D(0.0,1.0,0.0));
					wallHandler.getObject(1)->setOrientation(Vec3D(0.0,1.0,0.0));
					wallHandler.getObject(2)->setOrientation(Vec3D(0.0,1.0,0.0));
                                        wallHandler.getObject(3)->setOrientation(Vec3D(0.0,1.0,0.0));
					wallHandler.getObject(4)->setOrientation(Vec3D(0.0,1.0,0.0));
					wallHandler.getObject(5)->setOrientation(Vec3D(0.0,1.0,0.0));
					wallHandler.getObject(6)->setOrientation(Vec3D(0.0,1.0,0.0));
					wallHandler.getObject(7)->setOrientation(Vec3D(0.0,1.0,0.0));

				}
				else
				{
					checkTime = getTime() + 1.0;
				}
			}
		}
	}                
              
// ACTIONS BOFRE TIMESTEP *

// SEVERAL VOIDS
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
    

    
// SEVERAL VOIDS *
    
	private:

	double radS1,radS2;
	double rhoS1, rhoS2;
	double massS1, massS2;

	double CORDrum,CORS1,CORS2,tc;

	double sizeRatio;	
	double densityRatio;
	double drumFillFraction;
	double volumeFraction;
        double rad_encl, rad_excl;


	double particleWallFriction,particleParticleFriction;

	int numS1,numS2;
	int numS1ToBeInserted,numS2ToBeInserted;
	
	double drumRadius;
	double revolutionsPerSecond;
	double fractionalPolydispersity;

	int step;
	double checkTime;
        
  };
//! [T7:class]

int main(int argc, char *argv[])
{
  // Problem setup
  RotatingHexagon problem; // instantiate an object of class Tutorial 6

  problem.setName("RotatingHexagram_15deg");
//  problem.autoNumber();
  problem.setSystemDimensions(3);
  
  problem.setXMin(-0.0866025404*(5./3.));
  problem.setXMax(0.0866025404*(5./3.));
  problem.setYMin(-0.015);
  problem.setYMax(0.015);
  problem.setZMin(-0.075*(5./3.));  
  problem.setZMax(0.075*(5./3.));
  
    problem.setTimeMax(1000.);
    problem.setTimeStep(0.005/50.);
    problem.setDrumRadius(0.025);// in meters

    problem.setGravity(Vec3D(0.,0.,-9.81));
    problem.setCOR(0.97,0.831,0.831);//drumWall, species1, species2
    problem.setSizeAndDensityRatio(1.5,1.0);//size ratio and density ratio
    problem.setFractionalPolydispersity(0.0);//10% dispersity
    problem.setDrumFillFraction(0.05);// At 0.5 the drum is 3/4 filled. 
    problem.setSpeciesVolumeFraction(.5);//Species1 volume fraction
    problem.setFrictionCoeff(.6,.19);//(particle-wall-friction,part-part-friction)
    problem.setRevolutionSpeed(15);//rpm

    problem.setSaveCount(10000);
    problem.dataFile.setFileType(FileType::ONE_FILE);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.eneFile.setFileType(FileType::NO_FILE);

    problem.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    problem.setParticlesWriteVTK(true);

    problem.setXBallsAdditionalArguments("-cmode 8 -solidf -v0");
    problem.readArguments(argc,argv);

    problem.solve();

  return 0;
}
