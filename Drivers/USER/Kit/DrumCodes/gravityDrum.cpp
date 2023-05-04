//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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


/// \todo Make the scale factor of moving the surface of the drum particles closer together a setable paramter
/// \todo Make it four stages with a counter
/// Stage 1 : Create the drum
/// Stage 2 : Get the particles in : At the moment due to how to this checked only low particles numbers can be placed in.
/// Stage 3 : Settle the particles
/// Repeat stage 2 until all required particles are in
/// Stage 4 : Settle to a very low KE
/// Stage 5 : Rotate the drum
/// \todo write a restarter for this code
/// \todo insert upto rmin not just the middle 30%
/// \todo Vary the shape down the length
/// \todo place more particle in high r and then can reduce the overlap a bit more
/// \bug This is now using precribed particles; but, if the cor=0.1 they still stick. These needs to be checked to see if there is a bug in the precribed particles.
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Chute.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"


/// \brief Class for simulation rotating drums whose outer walls are made of particles
class RotatingDrum : public Chute
{
public:
    /// \brief Default constructor call. Sets the diameter of the small particles to be one.
    RotatingDrum() : radius_s(0.5)
    {
        densityRatio=1.0;
	fractionalPolydispersity = 0.0;
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////
    /// /brief setupInitial conditions. Basics does step 1 only; creating the walls of the drum
    ////////////////////////////////////////////////////////////////////////////////////////
    void setupInitialConditions() override
    {      
        // Set the step counter to 1
        step=1;
        
        
        //Do default no overlap.
        //scaleFactor=2.0;
        
        // Set the total number of partiles to be inserted to the number of particles which should be in the drum (the drum is currently empty)
        numSmallToBeInserted=numSmall;
        numLargeToBeInserted=numLarge;
        
        // Set the particles sizes for the chute code, is this still needed?
        setInflowParticleRadius(radius_s, radius_l);
        
        // Set the non-dim density
        double rho1 = 1.0 * 6.0 / constants::pi;
        double rho2 = densityRatio*rho1;
        //double rho2 = 1.6 * 6.0 / constants::pi;

        // Set the collision time and coefficient of resitution.
        //double tc = 1.0 / 200.0;
        
        
        //Computer the mass of the partices
        double mass_small = 4 / 3 * constants::pi * pow(radius_s, 3.0) * rho1;
        double mass_large = 4 / 3 * constants::pi * pow(radius_l, 3.0) * rho2;
      
        speciesDrumWall = new LinearViscoelasticFrictionSpecies;
        speciesParticles1 = new LinearViscoelasticFrictionSpecies;
        speciesParticles2 = new LinearViscoelasticFrictionSpecies;
        
        speciesDrumWall->setDensity(rho1);
        speciesDrumWall->setCollisionTimeAndRestitutionCoefficient(tc, wallCOR, mass_small); //MD::setCollisionTimeAndRestitutionCoefficient(tc, r,mass_small);
        speciesDrumWall->setSlidingDissipation(speciesDrumWall->getDissipation()); //  Set the tangential dissipation equal to the normal disipation for small-small collsions
        speciesDrumWall->setSlidingStiffness(speciesDrumWall->getStiffness()*2.0/7.0);
        speciesDrumWall->setSlidingFrictionCoefficient(4.0);
        //////
        //
        speciesParticles1->setDensity(rho1);
        speciesParticles1->setCollisionTimeAndRestitutionCoefficient(tc, smallCOR, mass_small);
        speciesParticles1->setSlidingDissipation(speciesParticles1->getDissipation()); //  Set the tangential dissipationequal to the normal disipation for large-large collision
        speciesParticles1->setSlidingStiffness(speciesParticles1->getStiffness()*2.0/7.0);
        speciesParticles1->setSlidingFrictionCoefficient(1.);

	speciesParticles1->setRollingStiffness(speciesParticles1->getStiffness()*2.0/5.0);
	speciesParticles1->setRollingFrictionCoefficient(0.002);
	speciesParticles1->setRollingDissipation(speciesParticles1->getDissipation());

	speciesParticles1->setTorsionStiffness(speciesParticles1->getStiffness()*2.0/5.0);
	speciesParticles1->setTorsionFrictionCoefficient(0.00001);
	speciesParticles1->setTorsionDissipation(speciesParticles1->getDissipation());
        //
        //////////
        speciesParticles2->setDensity(rho2);
        speciesParticles2->setCollisionTimeAndRestitutionCoefficient(tc, largeCOR, mass_large);
        speciesParticles2->setSlidingDissipation(speciesParticles2->getDissipation());
        speciesParticles2->setSlidingStiffness(speciesParticles2->getStiffness()*2.0/7.0);
        speciesParticles2->setSlidingFrictionCoefficient(1.);
        
	speciesParticles2->setRollingStiffness(speciesParticles1->getStiffness()*2.0/5.0);
	speciesParticles2->setRollingFrictionCoefficient(0.002);
	speciesParticles2->setRollingDissipation(speciesParticles2->getDissipation());

	speciesParticles2->setTorsionStiffness(speciesParticles1->getStiffness()*2.0/5.0);
	speciesParticles2->setTorsionFrictionCoefficient(0.00001);
	speciesParticles2->setTorsionDissipation(speciesParticles2->getDissipation());
        
        speciesHandler.addObject(speciesDrumWall);
        speciesHandler.addObject(speciesParticles1);
        speciesHandler.addObject(speciesParticles2);
        speciesMixedDrumAnd1 = speciesHandler.getMixedObject(speciesDrumWall, speciesParticles1);
        speciesMixedDrumAnd2 = speciesHandler.getMixedObject(speciesDrumWall, speciesParticles2);
        speciesMixed1And2 = speciesHandler.getMixedObject(speciesParticles1, speciesParticles2);
        
        
        //  Set the contact time (tc), resitution coefficeient (r) and density (rho) for small for all particles
        
   
        //speciesMixedDrumAnd1->setCollisionTimeAndRestitutionCoefficient(tc, (0.5*(wallCOR+smallCOR)), mass_small, mass_small);
	speciesMixedDrumAnd1->setCollisionTimeAndRestitutionCoefficient(tc, wallCOR, mass_small, mass_small);
        speciesMixedDrumAnd1->setSlidingDissipation(speciesMixedDrumAnd1->getDissipation()); //  Set the tangential dissipation equal to the normal disipation for mixed collision
	//Kit Hack
	//speciesMixedDrumAnd1->setSlidingDissipation(speciesDrumWall->getDissipation());
        speciesMixedDrumAnd1->setSlidingFrictionCoefficient(4.0);
        speciesMixedDrumAnd1->setSlidingStiffness(speciesMixedDrumAnd1->getStiffness()*2.0/7.0);
	//Kit Hack
	//speciesMixedDrumAnd1->setSlidingStiffness(speciesDrumWall->getStiffness()*2.0/7.0);     
	speciesMixedDrumAnd1->setRollingStiffness(speciesMixedDrumAnd1->getStiffness()*2.0/5.0);
	//Kit Hack
	//speciesMixedDrumAnd1->setRollingStiffness(speciesDrumWall->getStiffness()*2.0/5.0);
	speciesMixedDrumAnd1->setRollingFrictionCoefficient(0.0020);
	speciesMixedDrumAnd1->setRollingDissipation(speciesDrumWall->getDissipation());

	speciesMixedDrumAnd1->setTorsionStiffness(speciesParticles1->getStiffness()*2.0/5.0);
	//Kit Hackl
	//speciesMixedDrumAnd1->setTorsionStiffness(speciesDrumWall->getStiffness()*2.0/5.0);
	speciesMixedDrumAnd1->setTorsionFrictionCoefficient(0.00001);
	speciesMixedDrumAnd1->setTorsionDissipation(speciesDrumWall->getDissipation());

        //speciesMixedDrumAnd2->setCollisionTimeAndRestitutionCoefficient(tc, (0.5*(wallCOR+largeCOR)), mass_small, mass_large);
	speciesMixedDrumAnd2->setCollisionTimeAndRestitutionCoefficient(tc, wallCOR, mass_small, mass_large);
        speciesMixedDrumAnd2->setSlidingDissipation(speciesMixedDrumAnd2->getDissipation()); //  Set the tangential dissipation equal to the normal disipation
	//Kit Hack
	//speciesMixedDrumAnd2->setSlidingDissipation(speciesDrumWall->getDissipation());
        speciesMixedDrumAnd2->setSlidingFrictionCoefficient(4.0);
        speciesMixedDrumAnd2->setSlidingStiffness(speciesMixedDrumAnd2->getStiffness()*2.0/7.0);
        
	speciesMixedDrumAnd2->setRollingStiffness(speciesMixedDrumAnd2->getStiffness()*2.0/5.0);
	//Kit Hack
	//speciesMixedDrumAnd2->setRollingStiffness(speciesDrumWall->getStiffness()*2.0/5.0);
	speciesMixedDrumAnd2->setRollingFrictionCoefficient(0.0020);
	speciesMixedDrumAnd2->setRollingDissipation(speciesDrumWall->getDissipation());

	speciesMixedDrumAnd2->setTorsionStiffness(speciesMixedDrumAnd2->getStiffness()*2.0/5.0);
	//Kit Hack
	//speciesMixedDrumAnd2->setTorsionStiffness(speciesDrumWall->getStiffness()*2.0/5.0);
	speciesMixedDrumAnd2->setTorsionFrictionCoefficient(0.00001);
	speciesMixedDrumAnd2->setTorsionDissipation(speciesDrumWall->getDissipation());

        speciesMixed1And2->setCollisionTimeAndRestitutionCoefficient(tc, (0.5*(largeCOR+smallCOR)), mass_small, mass_large);
        speciesMixed1And2->setSlidingDissipation(speciesMixed1And2->getDissipation()); //  Set the tangential dissipation equal to the normal disipation
        speciesMixed1And2->setSlidingFrictionCoefficient(1.);
        speciesMixed1And2->setSlidingStiffness(speciesMixed1And2->getStiffness()*2.0/7.0);
                
	speciesMixed1And2->setRollingStiffness(speciesMixed1And2->getStiffness()*2.0/5.0);
	speciesMixed1And2->setRollingFrictionCoefficient(0.002);
	speciesMixed1And2->setRollingDissipation(speciesMixed1And2->getDissipation());

	speciesMixed1And2->setTorsionStiffness(speciesParticles1->getStiffness()*2.0/5.0);
	speciesMixed1And2->setTorsionFrictionCoefficient(0.00001);
	speciesMixed1And2->setTorsionDissipation(speciesParticles1->getDissipation());

        
        
        setXMin(0.0);
        setZMin(0.0);
        
        
        
        scaleFactor=8.0;
        // scalefactor 2 is fine for n-polygon but for stars it most be higher. 3.0 was not higher enought so when to 4;
        
        setXMax(drumRadius*constants::pi*2.0*scaleFactor);
        //Kit Hack -- originally no `scaleFactor' term in the below... Not sure if should/shouldn't be...
	setZMax(drumRadius*constants::pi*2.0*scaleFactor);
        Chute::setupInitialConditions();
        
        // Remove the two existing boundaries insertion and periodic and put the peroidic back if perodic else put
        boundaryHandler.clear();
        wallHandler.clear();
        
        if (getIsPeriodic())
        {
            PeriodicBoundary b0;
            b0.set(Vec3D(0.0, 1.0, 0.0), getYMin(), getYMax());
            boundaryHandler.copyAndAddObject(b0);
        }
        
      
        
        //Now add two extra solid walls at the end
        InfiniteWall w0;
        w0.set(Vec3D(0.0, -1.0, 0.0), getMin() + Vec3D(0,2 * getFixedParticleRadius(),0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 1.0, 0.0), getMax() - Vec3D(0,2 * getFixedParticleRadius(),0));
        wallHandler.copyAndAddObject(w0);
        
        
      
        
        setXMax(drumRadius*2.0);
        setZMax(drumRadius*2.0);
        
       
        //Now map the bed particles to the drum
        std::cout << "\n \n \n";
        std::cout << "STEP 1 : Creating the Drum " << std::endl;
        std::cout << "---------------------------" << std::endl;
        std::cout << "\n \n \n";
        
        Vec3D position;
        double r,r0,theta,y;
        
        BedParticlesInitialLocation.resize(particleHandler.getNumberOfObjects());
        for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
        {
            BaseParticle* P0 = particleHandler.getObject(i);
            position=P0->getPosition();
            theta=(position.X)/(drumRadius*scaleFactor);
            y=position.Y;
            r0=-(drumRadius*constants::pi*2.0-position.Z)/(drumRadius*constants::pi*2.0);
            
            //If simple drum
            if (nPolygon.size()==1)
            {
                r=getStarRCoord(nPolygon[0],mPolygon[0],r0,theta);
            }
            //else split drum
            else
            {
            
                if (y<0.5*getYMax())
                {
                    r=getStarRCoord(nPolygon[0],mPolygon[0],r0,theta);
                
                }
                else
                {
                    r=getStarRCoord(nPolygon[1],mPolygon[1],r0,theta);
                
                }
            }
            
            BedParticlesInitialLocation[i].X=theta;
            BedParticlesInitialLocation[i].Y=y;
            BedParticlesInitialLocation[i].Z=r;

            
            
            position.X=(r*cos(theta)+1+getXMin())*(getXMax()-getXMin())/2.0;
            position.Y=y;
            position.Z=(r*sin(theta)+1+getZMin())*(getZMax()-getZMin())/2.0;
            
            P0->setPosition(position);
            
        
            
        }
        
        //Now fill in the gap : If it is a 2 split drum
        if (nPolygon.size()==2)
        {
            BaseParticle* P0 = particleHandler.getObject(0);
            BaseParticle* P1 = P0->copy();
            int numWallParticles=particleHandler.getNumberOfObjects();
	        int numPartitionWallParticles = (particleHandler.getNumberOfObjects() / 12 ) * (particleHandler.getNumberOfObjects() / 12 );
            BedParticlesInitialLocation.resize(numWallParticles+numPartitionWallParticles);
            for (int i=0;i<numPartitionWallParticles;i++)
            {
            
                position.X=random.getRandomNumber(0.0,2.0*constants::pi);
                position.Y=getYMax()/2.0;
                position.Z=getZMin();
            
                theta=(position.X);
                y=position.Y;
                r0=-(drumRadius*constants::pi*2.0-position.Z)/(drumRadius*constants::pi*2.0);
                double rmin=getStarRCoord(5,2,r0,theta);
                double rmax=getStarRCoord(5,1,r0,theta);
            
                double r=random.getRandomNumber(rmin,rmax);
                
                BedParticlesInitialLocation[i+numWallParticles].X=theta;
                BedParticlesInitialLocation[i+numWallParticles].Y=y;
                BedParticlesInitialLocation[i+numWallParticles].Z=r;
            
                position.X=(r*cos(theta)+1+getXMin())*(getXMax()-getXMin())/2.0;
                position.Y=y;
                position.Z=(r*sin(theta)+1+getZMin())*(getZMax()-getZMin())/2.0;
            
            
                P1->setPosition(position);
            
                particleHandler.copyAndAddObject(P1);
            
                }
            delete P1;
        }
        
        // Set the step to step 2 and set the time to check if relaxed to 0.1 seconds from now.
        step=2;
        
        
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Create and inserts in the drum new partciles random at non-overlapping locations in the drum. If also check the locations is not already filled by a particle
    ////////////////////////////////////////////////////////////////////////////////////////////
    void createParticles()
    {
        

        
        std::cout << "\n \n \n";
        std::cout << "STEP 2: Inserting particles " << std::endl;
        std::cout << "---------------------------" << std::endl;
        std::cout << "\n \n \n";
        
        std::cout << "Number of large particles to be inserted:" << numLargeToBeInserted << std::endl;
        std::cout << "Number of small particles to be inserted:" << numSmallToBeInserted << std::endl;
        
        // CREATE THE PARTICLES
        SphericalParticle P0;
        
        Vec3D position;
        double r, theta, y;
        
        int failCounter=0;
       
        
        while ((numSmallToBeInserted > 0) || (numLargeToBeInserted > 0))
        {
            
            //random to see if want to generate a large or small particles, helps makes the initial conditions homogenious
            if (random.getRandomNumber(1.0, numLargeToBeInserted + numSmallToBeInserted) > (numLargeToBeInserted))
            {
	      //Kit hack - adding a random element to the radius of each particle in order to ensure
	      //a degree of polydispersity and hence eliminated unrealistic crystalline packing
	      P0.setRadius( radius_s + random.getRandomNumber( (-1 * fractionalPolydispersity * radius_s),(fractionalPolydispersity * radius_s) ) );
                P0.setSpecies(speciesParticles1);
                numSmallToBeInserted--;
            }
            else
            {
                P0.setRadius( radius_l + random.getRandomNumber( (-1 * fractionalPolydispersity * radius_l),( fractionalPolydispersity * radius_l ) ) );
                P0.setSpecies(speciesParticles2);
                numLargeToBeInserted--;
            }
            //randomise particle position, zero intial velocity
            failCounter=0;
            do
            {
                r=random.getRandomNumber(-0.3,0.3);
                theta=random.getRandomNumber(0,constants::pi*2.0);
                y=random.getRandomNumber(getYMin(),getYMax());
                
                position.X=(r*cos(theta)+1+getXMin())*(getXMax()-getXMin())/2.0;
                position.Y=y;
                position.Z=(r*sin(theta)+1+getZMin())*(getZMax()-getZMin())/2.0;
                
                P0.setPosition(position);
                P0.setVelocity(Vec3D(0.0, 0.0, 0.0));
                
                failCounter++;
                
                if (failCounter==1000) break;
                
              
            }
            while (!checkParticleForInteraction(P0));
            
            if (failCounter==1000) break;
                
            particleHandler.copyAndAddObject(P0);
            
        }
        
        std::cout << "Finished creating particles" << std::endl;
        std::cout << "Number of large particles still to be inserted:" << numLargeToBeInserted << std::endl;
        std::cout << "Number of small particles still to be inserted:" << numSmallToBeInserted << std::endl;
        
    
        
        if ((numSmallToBeInserted==0) && (numLargeToBeInserted==0))
            {
                // If you are here are partices are inserted and you are moving to step 4; relax to very low KE
                step=4;
                
                std::cout << "\n \n \n";
                std::cout << "STEP 4: Relaxing particles " << std::endl;
                std::cout << "---------------------------" << std::endl;
                std::cout << "\n \n \n";
                checkTime=getTime()+1.0;
                
            }
        else
            {
                // If you are here you are the dum is not full and you are settling the inserted particles
                
                step=3;
                
                std::cout << "\n \n \n";
                std::cout << "STEP 3: Settling the inserted particles " << std::endl;
                std::cout << "---------------------------" << std::endl;
                std::cout << "\n \n \n";
                checkTime=getTime()+1.0;
    
            }
        
    }
    
    
   
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// /brief actionsBeforeTimeStep. This does stage 2: insert particles, stage 3: setlle particles, stage 4: relax particles and stage 5: start the drum rotating.
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    void actionsBeforeTimeStep() override
    {
        if (step==2)
            createParticles();
        
        
        if (step==4 || step==3)
            {
            if (getTime() > checkTime)
                {
                    std::cout << "Current KE " << getKineticEnergy() << std::endl;
                    if (getKineticEnergy() < (particleHandler.getNumberOfObjects()/100.0))
                    {
                        if (step==4)
                        {
                            checkTime=getTime()+1.0;
                            if (getKineticEnergy() < (particleHandler.getNumberOfObjects()/1000.0))
                            {
                                step=5;
                                drumStartTime=getTime();
                                
                                std::cout << "\n \n \n";
                                std::cout << "STEP 5: Starting the drum rotation " << std::endl;
                                std::cout << "---------------------------" << std::endl;
                                std::cout << "\n \n \n";
                                checkTime=getTime()+1.0;
                                
                            }
                        }
                        else
                        {
                            step=2;
                        }
                    }
                    else
                    {
                        checkTime=getTime()+1.0;
                    }
                }
            }
        
            
        
        // Make a list of the partickes in the bed
        if (BedParticles.empty())
        {
            for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
            {
                BaseParticle* P0 = particleHandler.getObject(i);
                if (P0->getIndSpecies() == 0)
                {
                    P0->setVelocity(Vec3D(0.0,0.0,0.0));
                    BedParticles.push_back(P0);
                }
            }
            std::cout << "Finished storing bed particles" << std::endl;
        }
        
        
        
        // Move the particles which make up the drum
        if (step==5)
        {
            for (int i = 0; i < BedParticles.size(); i++)
            {
                
                BaseParticle* P0 = BedParticles[i];
                if (P0->getIndSpecies() == 0)
                {
                    /*
                    Vec3D position;
                    position = P0->getPosition();
                    
                    double theta, y, r;
                    
                    double Time=getTime();
                    
                    theta=BedParticlesInitialLocation[i].X;
                    y=BedParticlesInitialLocation[i].Y;
                    r=BedParticlesInitialLocation[i].Z;
                    
                    theta=theta+rpm*30.0*(Time-drumStartTime)/constants::pi; //rpm = linearSpeed/2*pi * 60.0

                    position.X=(r*cos(theta)+1+getXMin())*(getXMax()-getXMin())/2.0;
                    position.Y=y;
                    position.Z=(r*sin(theta)+1+getZMin())*(getZMax()-getZMin())/2.0;
                  
    //                 P0->setPosition(position);
    //                  */
                    
                     P0->fixParticle();
                     P0->setPrescribedPosition([this,i](double time)
                                             {
                                               
                                                  double Time=getTime();
                                                 
                                                  double theta, y, r;
                                                 
                                                  theta=BedParticlesInitialLocation[i].X;
                                                  y=BedParticlesInitialLocation[i].Y;
                                                  r=BedParticlesInitialLocation[i].Z;
                                                 
                                                  theta=theta+rpm*30.0*(Time-drumStartTime)/constants::pi; //rpm = linearSpeed/2*pi * 60.0
                                                 
                                            
                                                  return(Vec3D((r*cos(theta)+1+getXMin())*(getXMax()-getXMin())/2.0,
                                                               y,
                                                               (r*sin(theta)+1+getZMin())*(getZMax()-getZMin())/2.0));
                                              });

                 }
             }
         }
     }
    
    
    // /////////////
    // /// Function to map a coordinate on to a star polygon
    // ////////////
     double getStarRCoord(int n,int m, double r0, double theta0)
     {
        
         double theta,r;
         theta=theta0;
         r=r0;
        
         if (n == 0)
        {
            
             r=r0;
         }
        
         else
            
         {
            
             if (m==1)
             {
                 r=r0 * cos(constants::pi/n)/cos( ( fmod(theta, (2*constants::pi/n)) ) - constants::pi/n);
                
             }
             else
             {
                 double tmp=0.0;
                 double tmp2=1.0;
                 for (int i = 0 ; i<n; i++)
                 {
                     tmp=cos(constants::pi*m/n)/cos( ( fmod(theta, (2*constants::pi*m/n)) ) - constants::pi*m/n*i);
                    
                     if (tmp<=1.0) r=fmax(tmp,r);
                    
                 }
                 r=r0*r;
                
             }
            
         }
     return r;
     }
    
    // ////////////////
    
		void setRPM(double new_speed)
		{
		  rpm = new_speed;
		}
    void set_particle_numbers(int new_num_small, int new_num_large)
    {
        
        if (new_num_small > 0)
        {
            numSmall = new_num_small;
        }
        else
        {           
            std::cerr << "Please give a positive numnber if small particles" << std::endl;
        }
        
        if (new_num_large > 0)
        {
            numLarge = new_num_large;
        }
        else
        {           
            std::cerr << "Please give a positive numnber if small particles" << std::endl;
        }
        
    }
    void set_particle_numbers(int new_num_small)
    {
        if (new_num_small > 0)
        {          
            numSmall = new_num_small;
            numLarge = std::pow(radius_s / radius_l, 3) * numSmall * particleVolRatio;
        }
        else
        {
            std::cerr << "Please give a positve number of particels" << std::endl;
        }
        
    }
    
    void set_radiusLarge(double new_large_radius)
    {
        
        if (new_large_radius > 0)
        {           
            radius_l = new_large_radius;    
        }
        else
        {           
            std::cerr << "Radius must be greater than zero" << std::endl;
        }
    }
    
    void set_particle_number_volRatio(double new_volume_ratio)
    {
        particleVolRatio = new_volume_ratio;
    }
//Auto set now so shoud not be set by user
//    void setDrumStartTime(double time)
//    {
//        drumStartTime=time;
//    }
    
    void setDrumRadius(double radius)
    {
        drumRadius=radius;
    }
    
    /// \todo Should check n and m are coprime or both 0.
    void setStarShape(int new_n, int new_m)
    {
        nPolygon.resize(1);
        mPolygon.resize(1);
        nPolygon[0]=new_n;
        mPolygon[0]=new_m;
    }
    
    ///Set teh star shape for split drums
    void setStarShape(int n1,int m1,int n2,int m2)
    {
        nPolygon.resize(2);
        mPolygon.resize(2);
        nPolygon[0]=n1;
        mPolygon[0]=m1;
        nPolygon[1]=n2;
        mPolygon[1]=m2;
    }

    /// Set the particle size of the wall particles in dimeter and overlap between
    void setWallParameters(double particlesSize, double new_overlap)
    {
       // if (new_overlap*particlesSize<1.0)
        {
            setFixedParticleRadius(particlesSize*0.5);
            scaleFactor=new_overlap;
        }
        //else
        {
            
          //  std::cerr << "Requested paramters not possible" <<std::endl;
          //  exit(-1);
        }
    }
    
    void setDensityRatio(double newDensity)
    {
        
        densityRatio = newDensity;
        
    }
    
    void setCoefficientOfRestitutionLarge(double newCOR)
    {
        
      largeCOR=newCOR;
        
    }

  void setCoefficientOfRestitutionSmall(double newCOR)
    {
        
      smallCOR=newCOR;
        
    }

  void setCoefficientOfRestitutionWall(double newCOR)
  {
    
    wallCOR=newCOR;
    
  }

  void setFractionalPolydispersity(double poly) {
    //setting fractionalPolydispersity to the **MAGNITUDE** of the value input
    //such that a negative value will not break the random function!
    fractionalPolydispersity = sqrt(poly*poly); 
  }

unsigned int num_restart_small;
unsigned int num_restart_large;
double tc;
    
private:
  std::vector<BaseParticle*> BedParticles;
  std::vector<Vec3D> BedParticlesInitialLocation;
  double rpm;
  double radius_l;
  double densityRatio;
  double smallCOR;
  double largeCOR;
  double wallCOR;
  const double radius_s;
  int numSmall;
  int numLarge;
  int numSmallToBeInserted;
  int numLargeToBeInserted;
  int step;
    
    //The numbers which store the shapes
    std::vector<int> nPolygon;
    std::vector<int> mPolygon;
    
    
    // This is the factor the particles in the drum are overlapped by 1.0 means no overlap; 2.0 means 50% overlap etc.
    double scaleFactor;
    
    double particleVolRatio;
    double drumStartTime;
    double drumRadius;
    double checkTime;
    //Kit hack: defining a parameter to vary the degree of randomised polydispersity possessed by particles such that unrealistic 
    //crystallisation may be avoided.
    //The number provided gives the maximal (+/-) size variation of particles as a fraction of the (mean) size declared by the user. 
  double fractionalPolydispersity;
  LinearViscoelasticFrictionSpecies* speciesDrumWall;
  LinearViscoelasticFrictionSpecies* speciesParticles1;
  LinearViscoelasticFrictionSpecies* speciesParticles2;
  LinearViscoelasticFrictionMixedSpecies* speciesMixedDrumAnd1;
  LinearViscoelasticFrictionMixedSpecies* speciesMixedDrumAnd2;
  LinearViscoelasticFrictionMixedSpecies* speciesMixed1And2;
  
    
    
  
    
};

int main(int argc, char *argv[])
{
    //Print description
    std::cout << std::endl << "Description: A quasi-2D moving-bed channel with walls on the left and right boundary." << std::endl;
    
    // Problem parameters
    RotatingDrum problem;
    problem.setName("t38-SF8-tc1-over-2000");
    problem.autoNumber();
    problem.setTimeMax(1000.0);
    problem.tc = 1. / 2000.0;
    //problem.tc = 1. / 200.0;
    problem.setTimeStep(problem.tc / (50.0));

    problem.set_radiusLarge(0.75);
    problem.setDensityRatio(1.0);
    //problem.setCoefficientOfRestitution(0.9);
    problem.setCoefficientOfRestitutionLarge(0.8);
    problem.setCoefficientOfRestitutionSmall(0.8);
    problem.setCoefficientOfRestitutionWall(0.8);
    problem.set_particle_number_volRatio(1.0); //volume ratio of large to small
    //problem.set_particle_numbers(16700, 16700);
    problem.set_particle_numbers(150,150);

    problem.setChuteAngleAndMagnitudeOfGravity(0.0, 1.0);
    problem.setRPM(0.1);
    
    
    // Chute properties : Simply remove the first line to add side walls.
   //problem.makeChutePeriodic();
    //problem.setYMax(40.5);
    problem.setYMax(8.0);
    //problem.setDrumRadius(35.5);
    problem.setDrumRadius(15.0);
    //defining the polydispersity of the system. The value given is the fraction of the 'mean', user-defined
    //particle radius by which the size of each particle species may vary in the maximal case
    //i.e. for each individual particle, r = r_mean +/- (fractionalPolydispersity * r_mean);
    problem.setFractionalPolydispersity(0.1);
    
    
    // scalefactor 2 is fine for n-polygon but for stars it most be higher. 3.0 was not higher enought so when to 4;
    /// \bug The scale factor set here gets ignore and the only used one is hard wired in the code. Just track this bug down later.
    problem.setWallParameters(1.75,2.9);
    
    //Swap the next two lines to swap between the different type of rought bottoms. Note MONOLAYER_ORDERED does not give enough friction for the circular case.
    //problem.setRoughBottomType(MULTILAYER);
    problem.setRoughBottomType(MONOLAYER_DISORDERED);
    //problem.setRoughBottomType(MONOLAYER_ORDERED);

    
    problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(2000, problem.getTimeMax(), problem.getTimeStep()));
    problem.setXBallsAdditionalArguments("-cmode 8");
    problem.readArguments(argc, argv);
    
    problem.setStarShape(0,0);


    
    if (argc > 4)
    {
        problem.num_restart_large=atoi(argv[3]);
        problem.num_restart_small=atoi(argv[4]);
    }
    else
    {
        problem.num_restart_large=0;
        problem.num_restart_small=0;
    }
    
    std::cout << problem.num_restart_small <<" " <<problem.num_restart_large << std::endl;

    problem.solve();
    
}
      
