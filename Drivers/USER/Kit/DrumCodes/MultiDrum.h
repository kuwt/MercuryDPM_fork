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

/**
* \todo Fix the other wall creation algorithm (no chute)
* \todo Make the initial step a settable command line parameter for easy restarts
* \todo Add the option to work both with an fixed angular density or fixed linear
*       on the walls.
*/
/// \todo Make the scale factor of moving the surface of the drum particles closer together a setable paramter
/// \todo insert upto rmin not just the middle 30%
/// \todo Vary the shape down the length
/// \todo place more particle in high r and then can reduce the overlap a bit more
/// \bug This is now using precribed particles; but, if the cor=0.1 they still stick. These needs to be checked to see if there is a bug in the precribed particles.
#include <Logger.h>
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Chute.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

constexpr static int SEGMENT_DIVIDER = (75 * 75);

/// \brief Class for simulation rotating drums whose outer walls are made of particles
class RotatingDrum : public Chute
{
public:
    /// \brief Default constructor call. Sets the diameter of the small particles to be one.
    RotatingDrum() : radius_s(0.5)
    {
        densityRatio=1.0;
        isPoly=false;
    }

    /**
    * \brief Serialize the drum to the restartfile
    * Serializes the drum to the restartfile. This is kind of tricky because the restartfile isn't
    * as sturdy as we would like.
    *
    * This is also step dependant to make it easier to hack together restartfiles in weird cases :)
    *
    * Because Chute::write / Chute::read is broken, we skip these and chain directly to the DPMBase::
    * calls.
    *
    * \args os An std::ostream to write to
    * \args writeAllParticles Whether or not to write a full restart (the answer is probably yes.)
    */
    void write(std::ostream& os, bool writeAllParticles = true) const override
    {
      os << "RPM " << rpm;
      os << " step " << step << '\n';
      switch (step) {
        case 1:
          logger(WARN, "You are restarting from step 1, this is not supported. Just restart your sim.");
        case 2:
        case 3:
        case 4:
          os << " numSmall " << numSmall;
          os << " numLarge " << numLarge;
          os << " numSmallToBeInserted " << numSmallToBeInserted;
          os << " numLargeToBeInserted " << numLargeToBeInserted;
          os << " isPoly " << isPoly;
          os << " radius_l " << radius_l;
          os << " COR " << COR;
          os << " densityRatio " << densityRatio << '\n';
        case 5:
        case 6:
          break;
        default:
          logger(FATAL, "Bug in Mercury! write() called before setupInitialConditions, or step is wrong (%)", step);
          break;
      }
      
      os << "ShapeCount " << nPolygon.size() << '\n';
      for (std::size_t i = 0 ;i < nPolygon.size(); i++) {
          os << " base " << nPolygon[i] << " skip " << mPolygon[i] << " len " << lPolygon[i] << '\n';
      }
      os << "VolRatio " << particleVolRatio;
      os << " drumStartTime " << drumStartTime;
      os << " drumRadius " << drumRadius;
      os << " checkTime " << checkTime << std::endl;
      DPMBase::write( os, writeAllParticles );
    }
    
    /**
    * \brief Deserialize the drum from the restartfile
    * Deserializes the drum from the restartfile.
    * \sa RotatingDrum::write(std::ostream&,bool)
    *
    * \args in An std::istream where to read from
    */
    void read(std::istream& in) override
    {
       particleHandler.clear();
       wallHandler.clear();
       std::string _dummy;
       in >> _dummy >> rpm;    //Read the RPM

       in >> _dummy >> step;   //And the current phase...
       logger( INFO, "Step = %", step );

       switch (step) {
         case 1:
         case 2:
         case 3:
         case 4:
           in >> _dummy >> numSmall;
           in >> _dummy >> numLarge;
           in >> _dummy >> numSmallToBeInserted;
           in >> _dummy >> numLargeToBeInserted;
           in >> _dummy >> isPoly;
           in >> _dummy >> radius_l;
           in >> _dummy >> COR;
           in >> _dummy >> densityRatio;
           break;
         case 5:
         case 6:
           break;
         default:
           logger(FATAL, "Step is not okay! (%)", step);
       }
       
       std::size_t count;
       in >> _dummy >> count;
       nPolygon.resize(count);
       mPolygon.resize(count);
       lPolygon.resize(count);
       for (std::size_t i = 0; i < count; i++) {
         in >> _dummy >> nPolygon[i];
         in >> _dummy >> mPolygon[i];
         in >> _dummy >> lPolygon[i];
       }
       
       in >> _dummy >> particleVolRatio;
       in >> _dummy >> drumStartTime;
       in >> _dummy >> drumRadius;
       in >> _dummy >> checkTime;
       
       DPMBase::read( in );
    }
    
    /**
    * \brief restores the data after a restart
    * Restores species for all the particles, as well as prescribed positions and such.
    * Also restores the step and continues exactly where we were just at.
    */
    void actionsOnRestart() override
    {
      speciesDrumWall      = dynamic_cast<LinearViscoelasticFrictionSpecies*>(speciesHandler.getObject(0));
      speciesParticles1    = dynamic_cast<LinearViscoelasticFrictionSpecies*>(speciesHandler.getObject(1));
      speciesParticles2    = dynamic_cast<LinearViscoelasticFrictionSpecies*>(speciesHandler.getObject(2));
      speciesMixedDrumAnd1 = speciesHandler.getMixedObject(speciesDrumWall, speciesParticles1);
      speciesMixedDrumAnd2 = speciesHandler.getMixedObject(speciesDrumWall, speciesParticles2);
      speciesMixed1And2    = speciesHandler.getMixedObject(speciesParticles1, speciesParticles2);
      
      if (step == 6)
      {
        step = 5; //Prescribe the positions...
      }

      double timeMoved;
      if (drumStartTime == 0) {
        timeMoved = 0;
      } else {
        timeMoved = getTime() - drumStartTime;
      }
      
      for (BaseParticle* p : particleHandler)
      {
        if (p->getIndSpecies() == 0) {
          bedParticles.push_back(p);
          
          double theta, r, y;
          
          Vec3D partPos = p->getPosition();
          
          //Transform back, take into account any movement of the drum
          //to keep the particles where they should be.
          
          partPos.X = partPos.X * 2.0 / ( getXMax() - getXMin() ) - getXMin() - 1;
          partPos.Z = partPos.Z * 2.0 / ( getZMax() - getZMin() ) - getZMin() - 1;
          
          y = partPos.Y;
          r = std::sqrt(partPos.X * partPos.X + partPos.Z * partPos.Z);
          theta = std::atan2( partPos.Z, partPos.X );
          
          theta=theta-rpm*30.0*(timeMoved)/constants::pi; //rpm = linearSpeed/2*pi * 60.0
          
          partPos.X = theta;
          partPos.Y = y;
          partPos.Z = r;
          
          bedParticlesInitialLocation.push_back(partPos);
        }
      }

    }
    
    /**************************************************************************
    * \brief setupInitial conditions. Basics does step 1 only; creating the walls of the drum
    *
    * Creates a chute, rescales it to cover the walls of the drum, chains setupInitial
    * and then rescales our system to the wanted dimensions. This is only step 1.
    *
    * \sa Chute::setupInitialConditions()
    **************************************************************************/
    void setupInitialConditions() override
    {     
        bool exit = false; 
        if ( nPolygon.size() != mPolygon.size() ||
             nPolygon.size() != lPolygon.size() ) {
          std::cerr << "nPolygon.size() != mPolygon.size() != lPolygon.size()!"
                  << "\nThis should never, ever, happen."
                  << "\ndid you modify this code?" << std::endl;
          exit = true;
        }
        if ( nPolygon.size() == 0 ) {
          std::cerr << "Oops! You forgot to set the drumtype." << std::endl;
          exit = true;
        }
        
        bool allPositive = true;
        bool allNegative = true;
        for ( double d : lPolygon ) {
          if (d < 0)
            allPositive = false;
          if (d > 0)
            allNegative = false;
        }
        if (!allPositive && !allNegative) {
          std::cerr << "You have both segments with negative and positive length."
                  << "\nNegative means figure it out yourself,"
                  << "\nWhile positive is a fixed length."
                  << "\nThese are mutaly exclusive."
                  << std::endl;
          exit = true;
        }
        
        if (exit) {
          std::exit(5);
        }
        // Set the step counter to 1
        step=1;
        
        drumStartTime = 0;
        //Do default no overlap.
        //scaleFactor=1.0;
        
        // Set the total number of partiles to be inserted to the number of particles which should be in the drum (the drum is currently empty)
        numSmallToBeInserted=numSmall;
        numLargeToBeInserted=numLarge;
        
        // Set the particles sizes for the chute code, is this still needed?
        setInflowParticleRadius(radius_s, radius_l);
        
        // Set the non-dim density
        double rho1 = 6.0 / constants::pi;
	double rho2 = 6.0 / constants::pi;
        //double rho2 = densityRatio*rho1;
        
        // Set the collision time and coefficient of resitution.
        double tc = 1.0 / 200.0;
        
        
        //Computer the mass of the partices
        double mass_small = 4 / 3 * constants::pi * pow(radius_s, 3.0) * rho1;
        double mass_large = 4 / 3 * constants::pi * pow(radius_l, 3.0) * rho2;
        
        speciesDrumWall = new LinearViscoelasticFrictionSpecies;
        speciesParticles1 = new LinearViscoelasticFrictionSpecies;
        speciesParticles2 = new LinearViscoelasticFrictionSpecies;
        
        speciesDrumWall->setDensity(rho1);
        speciesDrumWall->setCollisionTimeAndRestitutionCoefficient(tc, COR, mass_small); //MD::setCollisionTimeAndRestitutionCoefficient(tc, r,mass_small);
        speciesDrumWall->setSlidingDissipation(speciesDrumWall->getDissipation()); //  Set the tangential dissipation equal to the normal disipation for small-small collsions
        speciesDrumWall->setSlidingStiffness(speciesDrumWall->getStiffness()*2.0/7.0);
        speciesDrumWall->setSlidingFrictionCoefficient(0.5);
        //////
        //
        speciesParticles1->setDensity(rho1);
        speciesParticles1->setCollisionTimeAndRestitutionCoefficient(tc, COR, mass_small);
        speciesParticles1->setSlidingDissipation(speciesParticles1->getDissipation()); //  Set the tangential dissipationequal to the normal disipation for large-large collision
        speciesParticles1->setSlidingStiffness(speciesParticles1->getStiffness()*2.0/7.0);
        speciesParticles1->setSlidingFrictionCoefficient(0.5);
        //
        //////////
        speciesParticles2->setDensity(rho2);
        speciesParticles2->setCollisionTimeAndRestitutionCoefficient(tc, COR, mass_large);
        speciesParticles2->setSlidingDissipation(speciesParticles2->getDissipation());
        speciesParticles2->setSlidingStiffness(speciesParticles2->getStiffness()*2.0/7.0);
        speciesParticles2->setSlidingFrictionCoefficient(0.5);
        
        
        speciesHandler.addObject(speciesDrumWall);
        speciesHandler.addObject(speciesParticles1);
        speciesHandler.addObject(speciesParticles2);
        speciesMixedDrumAnd1 = speciesHandler.getMixedObject(speciesDrumWall, speciesParticles1);
        speciesMixedDrumAnd2 = speciesHandler.getMixedObject(speciesDrumWall, speciesParticles2);
        speciesMixed1And2 = speciesHandler.getMixedObject(speciesParticles1, speciesParticles2);
        
        
        //  Set the contact time (tc), resitution coefficeient (r) and density (rho) for small for all particles
        
   
        speciesMixedDrumAnd1->setCollisionTimeAndRestitutionCoefficient(tc, COR, mass_small, mass_small);
        speciesMixedDrumAnd1->setSlidingDissipation(speciesMixedDrumAnd1->getDissipation()); //  Set the tangential dissipation equal to the normal disipation for mixed collision
        speciesMixedDrumAnd1->setSlidingFrictionCoefficient(0.5);
        speciesMixedDrumAnd1->setSlidingStiffness(speciesMixedDrumAnd1->getStiffness()*2.0/7.0);
                
        speciesMixedDrumAnd2->setCollisionTimeAndRestitutionCoefficient(tc, COR, mass_small, mass_large);
        speciesMixedDrumAnd2->setSlidingDissipation(speciesMixedDrumAnd2->getDissipation()); //  Set the tangential dissipation equal to the normal disipation
        speciesMixedDrumAnd2->setSlidingFrictionCoefficient(0.5);
        speciesMixedDrumAnd2->setSlidingStiffness(speciesMixedDrumAnd2->getStiffness()*2.0/7.0);
                
        speciesMixed1And2->setCollisionTimeAndRestitutionCoefficient(tc, COR, mass_small, mass_large);
        speciesMixed1And2->setSlidingDissipation(speciesMixed1And2->getDissipation()); //  Set the tangential dissipation equal to the normal disipation
        speciesMixed1And2->setSlidingFrictionCoefficient(0.5);
        speciesMixed1And2->setSlidingStiffness(speciesMixed1And2->getStiffness()*2.0/7.0);
                

        
        
        setXMin(0.0);
        setZMin(0.0);
        
        
        //Put this back to 2.9
        scaleFactor=14;
        // scalefactor 2 is fine for n-polygon but for stars it most be higher. 3.0 was not higher enought so when to 4;
        
        
        //---------------------[ Magic goes here ]-------------------
        // So, in order to create the walls, we take a chute with a bed of particles.
        // Then, we just remap the coordinates of the fixed particles to form a nice drum!
        // First, we need to beat the chute with a stick until it spawns the particles in a
        // convenient configuration. We do this by setting the X & Z maxima to convenient
        // constants, like something with 2pi!
        //
        // After that, we have to figure out the correct scalings so we have a constant
        // particle density! We throw in a couple of walls or periodic boundaries based on
        // settings.
        
        //So, we temporary set the XMax and ZMax for the chute setupInitialConditions!
        setXMax(drumRadius*constants::pi*2.0*scaleFactor);
        setZMax(drumRadius*constants::pi*2.0);
        Chute::setupInitialConditions();
        
        //But we don't like the walls of the Chute! So, we remove them again!
        boundaryHandler.clear();
        wallHandler.clear();
        
        //Well, are we periodic?
        if (getIsPeriodic())
        {
            //Aw yiss! Mothafuckin PeriodicBoundaries!
            PeriodicBoundary b0;
            b0.set(Vec3D(0.0, 1.0, 0.0), getYMin(), getYMax());
            boundaryHandler.copyAndAddObject(b0);
        } else {
            //Nope.gif -> We set up two walls at the ends!
            InfiniteWall w0;
            w0.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0, getYMin(), 0));
            wallHandler.copyAndAddObject(w0);
            w0.set(Vec3D(0.0,  1.0, 0.0), Vec3D(0, getYMax(), 0));
            wallHandler.copyAndAddObject(w0);
        }
        
        //Now, we've fooled Chute into creating the bed. Now we move the
        //X & Z maxima to sane numbers        
        setXMax(drumRadius*2.0);
        setZMax(drumRadius*2.0);
        
        //If we're all positive, and our Y doesn't match, rescale it.
        if (allPositive)
        {
          double totalLen = 0;
          for ( double l : lPolygon )
          {
            totalLen += l;
          } 
          
          double diff = (getYMax() - getYMin() - totalLen) / totalLen;
          if (std::abs(diff) > 1e-6) {
            logger(WARN, "You system doesn't add up! I'm rescaling yMax to make your system fit!\n"
                            "(segments added up: %, while ymax - ymin: %!)",
                            totalLen, getYMax() - getYMin() );            
            setYMax(getYMin() + totalLen);
          }
        } else if (allNegative) {
          double totalLen = 0;
          for ( double l : lPolygon )
          {
            totalLen -= l;
          }
          
          //Now, figure out unit length.
          double unitLen = -(getYMax() - getYMin()) / totalLen;
          for ( double& l : lPolygon ) //and rescale them!
          {
            l *= unitLen;
          }
        }
        
        //Move the particles into our actual drum shape!
        std::cout << "\n \n \n";
        std::cout << "STEP 1 : Creating the Drum " << std::endl;
        std::cout << "---------------------------" << std::endl;
        std::cout << "\n \n \n";

        // So, we just use some temporary variables to make our live slightly easier        
        Vec3D position;
        double r,r0,theta,y;
        
        //The part of the drum which this particle falls in
        int drumPart;
        double lengthSum;
        
        /*
        ** So, what we do here is essentially the following mapping
        **
        **   ############            ___  
        **  /##########/#    --->   /   \   
        ** +----------+## y        |\___/| 
        ** | #########|##          ||___|| (r,y,theta)
        ** | #########|##          ||   ||       
        ** |/         |/ z          \___/
        ** +----------+
        **     x    (x,y,z)
        **  With the particles who were stuck, now on the wall!
        */
        
        
        //We want to store the initial locations to once again, make our life
        // a lot easier.
        bedParticlesInitialLocation.resize(particleHandler.getNumberOfObjects());
        //Iterate over all the particles so far, as they're all part of the wall!
        for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
        {
            //Retrieve the particle...
            BaseParticle* part = particleHandler.getObject(i);
            position= part->getPosition();
            //X gets transformed into the angle, radius follows from math
            //and Y gets into Y!
            theta = (position.X)/(drumRadius*scaleFactor);
            y     =  position.Y;
            r0    = -1 + position.Z/(drumRadius*constants::pi*2.0);
            
            //Divide it into even parts!
            
//            drumPart = (int)((y / getYMax()) * nPolygon.size());
            lengthSum = 0;
            drumPart = -1; //because we start at 0...
            do {
              drumPart++;
              lengthSum += lPolygon[drumPart];
            } while ( lengthSum < y );
            
            
            if (drumPart >= nPolygon.size()) {
              std::cerr << "Numerical rounding is getting the best of you!" << std::endl;
              drumPart = nPolygon.size() - 1;
            }
            
            r = getStarRCoord(nPolygon[drumPart],mPolygon[drumPart], r0, theta);
            
            bedParticlesInitialLocation[i].X=theta;
            bedParticlesInitialLocation[i].Y=y;
            bedParticlesInitialLocation[i].Z=r;

            position = mapToCartesian(theta, y, r);
            
            part->setPosition(position);
        }
        
        //Now, we need to fill in the holes between the seperate parts!
        int numWallParticles = particleHandler.getNumberOfObjects();
        
        //This seems to be a good constant. Why? Nobody knows!
        int numSegmentWallParticles = (numWallParticles * numWallParticles) / ( SEGMENT_DIVIDER );
        
        //We have N-1 seperator walls, where N = nPolygon.size()..
        bedParticlesInitialLocation.resize(numWallParticles + numSegmentWallParticles * (nPolygon.size() - 1));
        
        BaseParticle* original = particleHandler.getObject(0);
 
        for (int i = 1; i < nPolygon.size(); i++) {
            BaseParticle* wallPart = original->copy(); //We copy!
            //and now we tune!...
            
            lengthSum = 0;
            for (int j = 0; j < i; j++) {
              lengthSum += lPolygon[j];
            }
            
            for (int j = 0; j < numSegmentWallParticles; j++) {
                theta = random.getRandomNumber(0.0, 2.0 * constants::pi);
                y = lengthSum; //getYMax() / nPolygon.size() * i;
                r0 = -1 + getZMin() / ( drumRadius * constants::pi * 2.0);
             
                //Just fill the space between the two segments! (hence rmin / rmax) 
                double rfirst  = getStarRCoord(nPolygon[i-1], mPolygon[i-1], r0, theta);
                double rsecond = getStarRCoord(nPolygon[i]  , mPolygon[i]  , r0, theta);
                double r = random.getRandomNumber(rfirst, rsecond);
              
                Vec3D & initialPos = bedParticlesInitialLocation[(i-1)*numSegmentWallParticles + numWallParticles + j];
                initialPos.X = theta;
                initialPos.Y = y;
                initialPos.Z = r;
              
                position = mapToCartesian(theta, y, r);
                wallPart->setPosition(position);
                particleHandler.copyAndAddObject(wallPart);
            }
            delete wallPart;
            
        }
        
        // Set the step to step 2 and set the time to check if relaxed to 0.1 seconds from now.
        step=2;
        
        
    }
    
    
    /**************************************************************************
    * \brief Convert cylindrical into carthesian
    * Utility function to map from cylindrical coordinates
    * to our nicely defined cartesian coordinates, taking our
    * system dimensions into account
    *
    * \param r radius
    * \param y distance on Y-axis
    * \param theta angular component
    * \returns Vec3D with the now carthesian coordinates!
    *
    **************************************************************************/
    Vec3D mapToCartesian(double theta, double y, double r) const
    {
        Vec3D retVal;
        retVal.X=(r*cos(theta)+1+getXMin())*(getXMax()-getXMin())/2.0;
        retVal.Y=y;
        retVal.Z=(r*sin(theta)+1+getZMin())*(getZMax()-getZMin())/2.0;
      
        return retVal;

    }
    
    /**************************************************************************
    * \brief Create and inserts particles.
    *
    * Creates and inserts particles in the drum new partciles random at
    * non-overlapping locations in the drum. If also check the locations is 
    * not already filled by a particle
    *
    * This function will only be called from step 2.
    **************************************************************************/
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
        double r, r0, theta, y;
        
        //Initial radius for scaling purposes
        r0 = -1 + getZMin() / ( drumRadius * constants::pi * 2.0);

        //Count the failures in case we do something stupid, so we won't be looping infinite.
        int failCounter=0;
        
        while ((numSmallToBeInserted > 0) || (numLargeToBeInserted > 0))
        {
            
            //random to see if want to generate a large or small particles, helps makes the initial conditions homogenious
            if (random.getRandomNumber(1.0, numLargeToBeInserted + numSmallToBeInserted) > (numLargeToBeInserted))
            {
                P0.setRadius(radius_s);
                P0.setSpecies(speciesParticles1);
                numSmallToBeInserted--;
            }
            else
            {
                P0.setRadius(radius_l);
                P0.setSpecies(speciesParticles2);
                numLargeToBeInserted--;
            }
            /// \todo horrible code I know but it works will clean up later
            if (isPoly)
            {
                P0.setRadius(random.getRandomNumber(radius_s,radius_l));
            }
            //randomise particle position, zero intial velocity
            failCounter=0;
            do
            {
                theta = random.getRandomNumber(0,constants::pi*2.0);
                y     = random.getRandomNumber(getYMin(),getYMax());
                
                int drumPart = (int)((y / getYMax()) * nPolygon.size());
                double maxR = getStarRCoord(nPolygon[drumPart], mPolygon[drumPart], r0, theta);
                r     = random.getRandomNumber(0, maxR);
        
                                
                position = mapToCartesian(theta, y, r);
                                
                P0.setPosition(position);
                P0.setVelocity(Vec3D(0.0, 0.0, 0.0));
                
                failCounter++;
                
                if (failCounter==1000) break;
                //Don't insert behind walls :)
                if (maxR < r)
                  continue;
              
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
   
    /**************************************************************************
    * /brief Depending on the stage, check and move stages.
    *
    * Depending on stage do the following:
    * 1. Can't happen (after SetupInitialConditions we are stage 2)
    * 2. Insert more particles. (moves to 3 or 4)
    * 3. Settle particles       (moves to 4)
    * 4. Relax particles        (moves to 2 or 5)
    * 5. Start rotation         (moves to 6)
    * 6. Sit back and relax.    (doesn't move)
    **************************************************************************/
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
        if (bedParticles.empty())
        {
            for (int i = 0; i < particleHandler.getNumberOfObjects(); i++)
            {
                BaseParticle* P0 = particleHandler.getObject(i);
                if (P0->getIndSpecies() == 0)
                {
                    P0->setVelocity(Vec3D(0.0,0.0,0.0));
                    bedParticles.push_back(P0);
                }
            }
            std::cout << "Finished storing bed particles" << std::endl;
        }
        
        
        
        // Move the particles which make up the drum
        if (step==5)
        {
            std::cout << "Moving drum..." << std::endl;
            std::cout << drumStartTime << std::endl;
            for (int i = 0; i < bedParticles.size(); i++)
            {
                
                BaseParticle* P0 = bedParticles[i];
                if (P0->getIndSpecies() == 0)
                {
                    P0->fixParticle();
                    P0->setPrescribedPosition([this,i](double time)
                                             {
                                               
                                                 double Time=getTime();
                                                 
                                                 double theta, y, r;
                                                 
                                                 theta = bedParticlesInitialLocation[i].X;
                                                 y     = bedParticlesInitialLocation[i].Y;
                                                 r     = bedParticlesInitialLocation[i].Z;
                                                 
                                                 theta=theta+rpm*30.0*(Time-drumStartTime)/constants::pi; //rpm = linearSpeed/2*pi * 60.0
                                                 
                                            
                                                 return mapToCartesian(theta, y, r);
                                             });

                }
            }
            step = 6;
        }
    }
    
    
    /////////////
    /// Function to map a coordinate on to a star polygon
    ////////////
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
                    tmp=cos(constants::pi*m/n)/cos( ( fmod(theta, (2*constants::pi*m/n)) ) - i * constants::pi*2*(m-1)/n);
                    
                    if (tmp<=1.0) r=fmax(tmp,r);
                    
                }
                r=r0*r;
                
            }
            
        }
    return r;
    }
    
    ////////////////
    
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
    
    void setDrumRadius(double radius)
    {
        drumRadius=radius;
    }
    
    /// \todo Should check n and m are coprime or both 0.
    MERCURY_DEPRECATED void setStarShape(int new_n, int new_m)
    {
        nPolygon.resize(1);
        mPolygon.resize(1);
        lPolygon.resize(1);
        nPolygon[0]=new_n;
        mPolygon[0]=new_m;
        lPolygon[0]= -1;
    }
    
    /**
    * \brief clears the current drum state
    * In case you did something wrong, clears the drum state.
    * please beware that this will not influence the simulation
    * after step 1.
    */
    void clearSegments()
    {
        nPolygon.clear();
        mPolygon.clear();
        lPolygon.clear();
    }
    
    /**
    * \brief add a segment
    *
    * Adds a segment with n outer vertices, m skips per connection
    * and length l. If all the lengths are -1, it will be scaled
    * such that the whole drum is the size of zmax.
    * Please beware that this does not influence the simulation
    * after setupInitialConditions
    *
    * \args n outer vertices
    * \args m skips for connecting
    * \args l length
    *
    * \bug Should check for coprimeness of n and m
    * \bug Some higher order cases fail for a weird reason
    * \sa setupInitialConditions()
    * \sa clearSegments()
    */
    void addSegment(int n, int m, double l = -1)
    {
        nPolygon.push_back(n);
        mPolygon.push_back(m);
        lPolygon.push_back(l);
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
    
    void setCoefficientOfRestitution(double newCOR)
    {
        
        COR=newCOR;
        
    }
    
    void makePoly()
    {
        isPoly=true;
    }
    
    unsigned int num_restart_small;
    unsigned int num_restart_large;
    
private:
    std::vector<BaseParticle*> bedParticles;
    std::vector<Vec3D> bedParticlesInitialLocation;
    double rpm;
    double radius_l;
    double densityRatio;
    double COR;
    const double radius_s;
    int numSmall;
    int numLarge;
    int numSmallToBeInserted;
    int numLargeToBeInserted;
    int step;
    
    bool isPoly;
    
    //The numbers which store the shapes
    std::vector<int>    nPolygon;
    std::vector<int>    mPolygon;
    std::vector<double> lPolygon;
    
    
    // This is the factor the particles in the drum are overlapped by 1.0 means no overlap; 2.0 means 50% overlap etc.
    double scaleFactor;
    
    double particleVolRatio;
    double drumStartTime;
    double drumRadius;
    double checkTime;
    
    LinearViscoelasticFrictionSpecies* speciesDrumWall;
    LinearViscoelasticFrictionSpecies* speciesParticles1;
    LinearViscoelasticFrictionSpecies* speciesParticles2;
    LinearViscoelasticFrictionMixedSpecies* speciesMixedDrumAnd1;
    LinearViscoelasticFrictionMixedSpecies* speciesMixedDrumAnd2;
    LinearViscoelasticFrictionMixedSpecies* speciesMixed1And2;
};

