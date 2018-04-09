//Copyright (c) 2013-2017, The MercuryDPM Developers Team. All rights reserved.
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

#include "Mercury3D.h"
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include "Coil.h"

/** This code creates a cylindrical container, inserts particles and lets them settle.
*/
class Booster : public Mercury3D{
public:

    //set default values
    Booster()
    {
        setName("Booster");
        setFileType(FileType::ONE_FILE);
        setSaveCount(200);
        setTimeMax(1e20); //run forever
        revolutionsPerSecond = 0.0;
        drumInclination = 0.0;
        polyDisp_=false;
        fillingFraction_=0.5;

        //creating species for particle-particle interactions and particle-wall interactions
        particleSpecies = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        wallSpecies = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        wallParticleSpecies = speciesHandler.getMixedObject(wallSpecies, particleSpecies);

        drumInclination = 5.0*constants::pi/180.0;
        revolutionsPerSecond = 8.0/60.0;
        particleDiameter = 7.2e-2*2.0;
        particleSpecies->setDensity(2000.0);
        collisionTime = 0.01; //soft particles, but relative overlap due to gravity, o/d ~ (tc/tg)^2 = 0.0136, is still small
        restitutionCoefficient = 0.001; // very dissipative particles
        setXBallsAdditionalArguments(" -v0 -solidf  -w 1100 -s 0.76 -noborder 3");
        setTimeMax(5.0*60.0/8.0); //run 5 revolutions
    }

    //set slave variables
    void setupInitialConditions()
	{
        double volumeInserted=0;
        setGravity(9.8*Vec3D(std::sin(drumInclination),0.0,-std::cos(drumInclination)));

        Mdouble effectiveMass;
        Mdouble scaledFillHeight;
        //set contact properties:
        //- calculate stiffness and dissipation
        if (polyDisp_)
        {
            effectiveMass = 0.5*particleSpecies->getMassFromRadius(0.5*particleDiameter/255.0*255.0);
            scaledFillHeight=fillingFraction_/2.3338;
            
        }
        else
        {
            effectiveMass = 0.5*particleSpecies->getMassFromRadius(0.5*particleDiameter);
            scaledFillHeight=fillingFraction_;
        }
        
        //std::cout << "Mass " << 2.0*effectiveMass << std::endl;
        helpers::KAndDisp kAndDisp = helpers::computeKAndDispFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(collisionTime, restitutionCoefficient, effectiveMass);

        //- set stiffness, dissipation and timestep
        particleSpecies->setStiffness(kAndDisp.k);
        particleSpecies->setDissipation(kAndDisp.disp);
        particleSpecies->setSlidingStiffness(2.0/7.0*kAndDisp.k);
        particleSpecies->setSlidingDissipation(2.0/7.0*kAndDisp.disp);
        particleSpecies->setSlidingFrictionCoefficient(0.5);
        if (polyDisp_)
        {
            // Note this is 30^1.5 or more precisionly ((max size)/(min size))^1.5
	    // Actuall not (600/255)^1.5, think about it
            collisionTime=collisionTime/164.3168;
            setSaveCount(dataFile.getSaveCount()*164.3168);
        }
        setTimeStep(0.02*collisionTime);
      

        //set wall properties
        wallParticleSpecies->setStiffness(particleSpecies->getStiffness());
        wallParticleSpecies->setDissipation(particleSpecies->getDissipation());
        wallParticleSpecies->setSlidingStiffness(particleSpecies->getSlidingStiffness());
        wallParticleSpecies->setSlidingDissipation(particleSpecies->getSlidingDissipation());
        wallParticleSpecies->setSlidingFrictionCoefficient(particleSpecies->getSlidingFrictionCoefficient());

        Mdouble lidWidth = 216.73e-3;
        Mdouble upperWidth = 976.57e-3;
        Mdouble centreWidth = 2000.99e-3;
        Mdouble bottomWidth = 322.68e-3;
        Mdouble tripodWidth = bottomWidth; //check!
        Mdouble lidRadius = 300e-3;
        Mdouble centreRadius = 642.5e-3;
        Mdouble bottomRadius = 120e-3;

        setXMax(lidWidth+upperWidth+centreWidth+bottomWidth);
        setXMin(0.0);
        setYMax(centreRadius);
        setYMin(-getYMax());
        setZMax(getYMax());
        setZMin(getYMin());

        std::cout << "creating centre wall" << std::endl;
        //centre wall: cylindrical wall in the centre, r=centreRadius, lidWidth<x<lidWidth+centreWidth
        sideWall = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
        sideWall->setPosition(Vec3D(0.0,0.0,0.0) );
        sideWall->setOrientation(Vec3D(1.0,0.0,0.0) );
        sideWall->addObject(Vec3D(1,0,0), Vec3D(centreRadius,0.0,0.0) );
        sideWall->setSpecies(wallSpecies);
        sideWall->setAngularVelocity(Vec3D(revolutionsPerSecond * 2.0 * constants::pi,0.0,0.0));

        std::cout << "creating top wall" << std::endl;
        //top wall: conical wall near left end, lidRadius<r<centreRadius, 0<x<lidWidth (actually extends also to negative x)
        upperWall = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
        upperWall->setPosition(Vec3D(0.0,0.0,0.0) );
        upperWall->setOrientation(Vec3D(1.0,0.0,0.0) );
        upperWall->addObject(Vec3D(upperWidth-lidWidth,0,lidRadius-centreRadius), Vec3D(lidRadius,0.0,lidWidth) );
        upperWall->setSpecies(wallSpecies);
        upperWall->setAngularVelocity(Vec3D(revolutionsPerSecond * 2.0 * constants::pi,0.0,0.0));

        std::cout << "creating lower wall" << std::endl;
        //lower wall: conical wall near right end, bottomRadius<r<centreRadius, xMax-bottomWidth<x<xMax
        lowerWall = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
        lowerWall->setPosition(Vec3D(0.0,0.0,0.0) );
        lowerWall->setOrientation(Vec3D(1.0,0.0,0.0) );
        //the 30 is because the bending of the lower wall only starts after 30 mm (see sketch)
        ///\todo tw Walls still require that the normal is of unit length!!
        Vec3D normal = Vec3D(bottomWidth-30e-3,0,centreRadius-bottomRadius); normal.normalize();
        lowerWall->addObject(normal, Vec3D(bottomRadius,0.0,getXMax()) );
        lowerWall->setSpecies(wallSpecies);
        lowerWall->setAngularVelocity(Vec3D(revolutionsPerSecond * 2.0 * constants::pi,0.0,0.0));

        std::cout << "creating tripod base" << std::endl;
        //tripod
        tripodBase = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
        tripodBase->setPosition(Vec3D(0.0,0.0,0.0) );
        tripodBase->setOrientation(Vec3D(1.0,0.0,0.0) );
        //the 30 is because the bending of the lower wall only starts after 30 mm (see sketch)
        ///\todo tw Walls still require that the normal is of unit length!!
        tripodBase->addObject(Vec3D(-1.0,0.0,0.0)/*(r,0,x)*/, Vec3D(bottomRadius,0.0,0.0) );
        tripodBase->addObject(Vec3D( 0.0,0.0,1.0)/*(r,0,x)*/, Vec3D(0.0,0.0,getXMax()-tripodWidth) );
        tripodBase->setSpecies(wallSpecies);
        tripodBase->setAngularVelocity(Vec3D(revolutionsPerSecond * 2.0 * constants::pi,0.0,0.0));

        addBlades();

        // check if all particle diameter is set correctly
        if (particleDiameter==0.0)
        {
            std::cerr << "Error in booster::setupbooster: The particle diameter is not set." << std::endl;
            exit(-1);
        }

        // particles are set based on particle diameter
        BaseParticle P;
        P.setSpecies(particleSpecies);

        Vec3D pos;
        unsigned int numberOfParticlesInserted = 0;
        ///\todo TW 0.35 was 0.4 before
        unsigned int numberOfParticles = std::round(scaledFillHeight*(getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin())/mathsFunc::cubic(particleDiameter)/(4.0/3.14));
        hGridRebuild();
        //std::cout << std::endl << "Status before inserting particles:" << std::endl;
        //write(std::cout,false);
        
        // Create and insert the particles.
        std::cout << std::endl << "Inserting " << std::floor(numberOfParticles) << " particles" << std::endl;
        if (polyDisp_)
        {
            //This is crap code (as it is full of magic numbers) for FC and could be generised and hence used for a lot of problems
            
            //Step 1 : setup
            std::vector<int> particlesPerBin(6);
            std::vector<double> sizeRanges(7);
            /* Orginal data from compant
             * 10 percent less than 50 micron
             * 50 percent less than 225 micron
             * 90 percent less than 390 micron
             * 99 percent less than 570 micron
             *
             * Will add the extra fake data to complete teh information
             * 0 percent less than 20 micro
             * 100 percent less than 600 micro
             */
            particlesPerBin[0]=0.1*numberOfParticles;
            particlesPerBin[1]=0.4*numberOfParticles;
            particlesPerBin[2]=0.5*numberOfParticles;
            particlesPerBin[3]=0.4*numberOfParticles;
            particlesPerBin[4]=0.09*numberOfParticles;
            particlesPerBin[5]=0.01*numberOfParticles;
            
            /* Now set up the size ranges
             * Bin 0 has range 0-1
             * Bin 1 has range 1-2
             * Bin 2 has range 2-3
             * Bin 3 has range 3-4
             * Bin 4 has rnage 4-5
             * Bin 6 has range 5-6
             * Hence range 4 is the medium particles size and the one the users sets.
             */
            
        
            
            // Now rest are from the experimental data choosen to get the correct scaled size
            sizeRanges[0]=particleDiameter/225.0*20.0;
            sizeRanges[1]=particleDiameter/225.0*40.0;
            sizeRanges[2]=particleDiameter/225.0*50.0;
            sizeRanges[3]=particleDiameter/225.0*225.0;
            sizeRanges[4]=particleDiameter/225.0*390.0;
            sizeRanges[5]=particleDiameter/225.0*570.0;
            sizeRanges[6]=particleDiameter/225.0*600.0;
            
            //sizeRanges[0]=particleDiameter/225.0*225.0;
            //sizeRanges[1]=particleDiameter/225.0*225.0;
            //sizeRanges[2]=particleDiameter/225.0*225.0;
            //sizeRanges[3]=particleDiameter/225.0*225.0;
            //sizeRanges[4]=particleDiameter/225.0*390.0;
            //sizeRanges[5]=particleDiameter/225.0*570.0;
            //sizeRanges[6]=particleDiameter/225.0*600.0;
            
            std::cout << "About to create particles" << std::endl;
            int numberOfParticlesToInsert;
            do
            {
            
                //Step 2 Pick which type of particles to insert
               
                numberOfParticlesToInsert=0;
                for (int i=0; i<particlesPerBin.size();i++)
                {
                   // std::cout << "Need to insert " << particlesPerBin[i] << " particles between size " <<sizeRanges[i] <<" and " << sizeRanges[i+1] << std::endl;
                    numberOfParticlesToInsert += particlesPerBin[i];
                }
            
                int rand=random.getRandomNumber(1.0,numberOfParticlesToInsert);
            
                int pickedBin=-1;
            
           
                while (rand>0)
                {
                    pickedBin++;
                    //std::cout << pickedBin << ":" << rand<<std::endl;
                    rand=rand-particlesPerBin[pickedBin];
                
                }
                //std::cout << "Inserting particle in bin " <<pickedBin << std::endl;
        
            
                P.setRadius(0.5*random.getRandomNumber(sizeRanges[pickedBin],sizeRanges[pickedBin+1]));
            
                pos.X = random.getRandomNumber(getXMin()+P.getRadius(), getXMax()-P.getRadius());
                pos.Y = random.getRandomNumber(getYMin()+P.getRadius(), getYMax()-P.getRadius());
                pos.Z = random.getRandomNumber(getZMin()+P.getRadius(), getZMax()-P.getRadius());
                P.setPosition(pos);
                
              
                //std::cout << checkParticleForInteraction(P);
               // if (checkParticleForInteraction(P))
                {
                  //  std::cout << P.getRadius() << std::endl;
                    particleHandler.copyAndAddObject(P);
                    volumeInserted+=P.getRadius()*P.getRadius()*P.getRadius()*4.0/3.0*constants::pi;

                
                    ++numberOfParticlesInserted;
                    particlesPerBin[pickedBin]=particlesPerBin[pickedBin]-1;
                    //std::cout << "." << std::flush;
                }
            } while (numberOfParticlesToInsert>1);
        
        }
        else
        {
            P.setRadius(0.5*particleDiameter*random.getRandomNumber(0.95,1.05));
            
            while (numberOfParticlesInserted < numberOfParticles)
            {
                pos.X = random.getRandomNumber(getXMin()+P.getRadius(), getXMax()-P.getRadius());
                pos.Y = random.getRandomNumber(getYMin()+P.getRadius(), getYMax()-P.getRadius());
                pos.Z = random.getRandomNumber(getZMin()+P.getRadius(), getZMax()-P.getRadius());
                P.setPosition(pos);
                //std::cout << checkParticleForInteraction(P);
                //if (checkParticleForInteraction(P))
                //{
                    particleHandler.copyAndAddObject(P);
                    volumeInserted+=P.getRadius()*P.getRadius()*P.getRadius()*4.0/3.0*constants::pi;
                    P.setRadius(0.5*particleDiameter*random.getRandomNumber(0.95,1.05));
                
                    ++numberOfParticlesInserted;
                    //std::cout << "." << std::flush;
                //}
            }
        }
        std::cout << std::endl;
        std::cout << "Total volume insterted = " << volumeInserted << std::endl;
        std::cout << std::endl << "Starting simulation, timeMax= " << getTimeMax() << std::endl;
    }

    virtual void addBlades()
    {   }

    void printTime() const
    {
        std::cout
            << "t " << getTime()
            << " angle " <<  revolutionsPerSecond * getTime() * 2.0 * constants::pi
            << " ene " << getKineticEnergy()/getElasticEnergy() << std::endl;
    }

    LinearViscoelasticSlidingFrictionSpecies* particleSpecies;
    LinearViscoelasticSlidingFrictionMixedSpecies* wallParticleSpecies;
    LinearViscoelasticSlidingFrictionSpecies* wallSpecies;
    AxisymmetricIntersectionOfWalls* upperWall, * sideWall, *lowerWall, *tripodBase;
    InfiniteWall* baseWall;

    void actionsBeforeTimeStep()
    {
        sideWall->setOrientation(Vec3D(1.0,0.0,0.0) );
        upperWall->setOrientation(Vec3D(1.0,0.0,0.0) );
        lowerWall->setOrientation(Vec3D(1.0,0.0,0.0) );
        tripodBase->setOrientation(Vec3D(1.0,0.0,0.0) );
    }
    
    
    void makePolydispersed()
    {
        polyDisp_=true;
    }
    
    void setFillingFraction(double fillFraction)
    {
        fillingFraction_=fillFraction;
    }
    
    void actionsOnRestart()
    {
       
        
        //- set stiffness, dissipation and timestep
        particleSpecies = dynamic_cast<LinearViscoelasticSlidingFrictionSpecies*>(speciesHandler.getObject(0));
        wallSpecies =dynamic_cast<LinearViscoelasticSlidingFrictionSpecies*>(speciesHandler.getObject(1));
        wallParticleSpecies = dynamic_cast<LinearViscoelasticSlidingFrictionMixedSpecies*>(speciesHandler.getMixedObject(0,1));
        sideWall = dynamic_cast<AxisymmetricIntersectionOfWalls*>(wallHandler.getObject(0));
        upperWall = dynamic_cast<AxisymmetricIntersectionOfWalls*>(wallHandler.getObject(1));
        lowerWall = dynamic_cast<AxisymmetricIntersectionOfWalls*>(wallHandler.getObject(2));
        tripodBase = dynamic_cast<AxisymmetricIntersectionOfWalls*>(wallHandler.getObject(3));
        
       // addBlades();
    }


    double particleDiameter;
    double drumInclination;
    double revolutionsPerSecond;
    Mdouble collisionTime; //softness of particles
    Mdouble restitutionCoefficient; // dissipativeness of particles
    
private:
    bool polyDisp_; //flag to turn on true size distribution.
    double fillingFraction_;
};

class BoosterWithAngledBlade : public Booster {
public:

    BoosterWithAngledBlade(Mdouble bladeThickness, Mdouble bladeHeight, Mdouble bladeInclination, unsigned int bladeNumber)
        : Booster(), bladeThickness_(bladeThickness), bladeHeight_(bladeHeight), bladeInclination_(bladeInclination), bladeNumber_(bladeNumber)
    {
        setName("BoosterWithAngledBladeBladeNumber3");
    }

    void addBlades()
    {
        //blade
        bladesDown.reserve(bladeNumber_);
        bladesRight.reserve(bladeNumber_);
        bladesUp.reserve(bladeNumber_);
        bladesLeft.reserve(bladeNumber_);

        IntersectionOfWalls defaultBlade;
        defaultBlade.setSpecies(wallSpecies);

        Mdouble lidWidth = 216.73e-3;
        Mdouble upperWidth = 976.57e-3;
        Mdouble centreWidth = 2000.99e-3;

        double x = lidWidth+upperWidth;
        double dx = 0.5*centreWidth/static_cast<double>(bladeNumber_-0.5);
        for (unsigned int i=0; i<bladeNumber_; ++i)
        {
            //double x = (static_cast<double>(i)+0.5)*getXMax()/ static_cast<double>(bladeNumber_);
            bladesDown.push_back(wallHandler.copyAndAddObject(defaultBlade));
            bladesDown.back()->setPosition(Vec3D(x,0.0,0.0));
            bladesUp.push_back(wallHandler.copyAndAddObject(defaultBlade));
            bladesUp.back()->setPosition(Vec3D(x,0.0,0.0));
            x += dx;
            bladesRight.push_back(wallHandler.copyAndAddObject(defaultBlade));
            bladesRight.back()->setPosition(Vec3D(x,0.0,0.0));
            bladesLeft.push_back(wallHandler.copyAndAddObject(defaultBlade));
            bladesLeft.back()->setPosition(Vec3D(x,0.0,0.0));
            x += dx;
        }

        actionsAfterTimeStep();
    }

    virtual void actionsAfterTimeStep()
    {
        Mdouble angle = revolutionsPerSecond * getTime() * 2.0 * constants::pi;
        for (IntersectionOfWalls* blade : bladesDown)
            setBladeAngle(blade,angle);
        for (IntersectionOfWalls* blade : bladesRight)
            setBladeAngle(blade,angle+constants::pi/2.0);
        for (IntersectionOfWalls* blade : bladesUp)
            setBladeAngle(blade,angle+2.0*constants::pi/2.0);
        for (IntersectionOfWalls* blade : bladesLeft)
            setBladeAngle(blade,angle+3.0*constants::pi/2.0);
    }

    virtual void setBladeAngle(IntersectionOfWalls* blade, Mdouble angle)
    {
        blade->clear();
        Mdouble s = std::sin(angle), c = std::cos(angle);
        Vec3D anchorPosition = Vec3D(blade->getPosition().X, getZMax()*s, getZMax()*c); //clockwise as seen from lid
        Vec3D longAxis = Vec3D(0.0, s, c);
        Vec3D shortAxis = Vec3D(1./std::tan(bladeInclination_), -c, s); shortAxis.normalize();
        Vec3D mediumAxis = Vec3D::cross(shortAxis,longAxis);
        blade->addObject(longAxis, anchorPosition -bladeHeight_* longAxis);
        blade->addObject(-longAxis, anchorPosition +bladeHeight_* longAxis);
        blade->addObject( shortAxis, anchorPosition -0.5*bladeThickness_* shortAxis);
        blade->addObject(-shortAxis, anchorPosition +0.5*bladeThickness_* shortAxis);
    }

    void write(std::ostream& os, bool x) const
    {
        Mercury3D::write(os, x);
        os  << " bladeNumber " << bladeNumber_ << std::endl
        << " bladeInclination " << bladeInclination_ << std::endl
        << " bladeThickness " << bladeThickness_ << std::endl
        << " bladeHeight " << bladeHeight_ << std::endl;
    }
protected:
    // pointers to the blades, split into groups denoting the four initial positions of the anchor points
    std::vector<IntersectionOfWalls*> bladesDown, bladesRight, bladesUp, bladesLeft;
public:
    // inclination of the blade in the xy plane when the blade is at it's lowest point (when setBladeAngle::angle=0)
    Mdouble bladeInclination_;
    Mdouble bladeThickness_;
    Mdouble bladeHeight_;
    unsigned int bladeNumber_;
};

class BoosterWithAngledBladeUnstacked : public BoosterWithAngledBlade {
public:

    BoosterWithAngledBladeUnstacked(Mdouble bladeThickness, Mdouble bladeHeight, Mdouble bladeInclination, unsigned int bladeNumber)
        : BoosterWithAngledBlade(bladeThickness,bladeHeight,bladeInclination,bladeNumber)
    {
        setName("BoosterWithAngledBladeUnstacked");
    }

    void addBlades()
    {
        //blade
        bladesDown.reserve(bladeNumber_);
        bladesRight.reserve(bladeNumber_);
        bladesUp.reserve(bladeNumber_);
        bladesLeft.reserve(bladeNumber_);

        IntersectionOfWalls defaultBlade;
        defaultBlade.setSpecies(wallSpecies);

        Mdouble lidWidth = 216.73e-3;
        Mdouble upperWidth = 976.57e-3;
        Mdouble centreWidth = 2000.99e-3;

        double x = lidWidth+upperWidth;
        double dx = centreWidth/static_cast<double>(bladeNumber_+1);
        for (unsigned int i=0; i<bladeNumber_; ++i)
        {
            x += dx;
            //double x = (static_cast<double>(i)+0.5)*getXMax()/ static_cast<double>(bladeNumber_);
            bladesDown.push_back(wallHandler.copyAndAddObject(defaultBlade));
            bladesDown.back()->setPosition(Vec3D(x,0.0,0.0));
            bladesUp.push_back(wallHandler.copyAndAddObject(defaultBlade));
            bladesUp.back()->setPosition(Vec3D(x,0.0,0.0));
            bladesRight.push_back(wallHandler.copyAndAddObject(defaultBlade));
            bladesRight.back()->setPosition(Vec3D(x,0.0,0.0));
            bladesLeft.push_back(wallHandler.copyAndAddObject(defaultBlade));
            bladesLeft.back()->setPosition(Vec3D(x,0.0,0.0));
        }

        actionsAfterTimeStep();
    }
};

class BoosterWithSmallAngledBlade : public BoosterWithAngledBlade {
public:

    BoosterWithSmallAngledBlade(Mdouble bladeThickness, Mdouble bladeHeight, Mdouble bladeWidth, Mdouble bladeInclination, unsigned int bladeNumber)
        : BoosterWithAngledBlade(bladeThickness, bladeHeight, bladeInclination, bladeNumber), bladeWidth_(bladeWidth)
    {
        setName("BoosterWithSmallAngledBladeFixedBladeWidthAndInclination");
    }

    void setBladeAngle(IntersectionOfWalls* blade, Mdouble angle)
    {
        blade->clear();
        Mdouble s = std::sin(angle), c = std::cos(angle);
        Vec3D anchorPosition = Vec3D(blade->getPosition().X, getZMax()*s, getZMax()*c); //clockwise as seen from lid
        Vec3D longAxis = Vec3D(0.0, s, c);
        Vec3D shortAxis = Vec3D(1./std::tan(bladeInclination_), -c, s); shortAxis.normalize();
        Vec3D mediumAxis = Vec3D::cross(shortAxis,longAxis);
        blade->addObject(longAxis, anchorPosition -bladeHeight_* longAxis);
        blade->addObject(-longAxis, anchorPosition +bladeHeight_* longAxis);
        blade->addObject( shortAxis, anchorPosition -0.5*bladeThickness_* shortAxis);
        blade->addObject(-shortAxis, anchorPosition +0.5*bladeThickness_* shortAxis);
        blade->addObject( mediumAxis, anchorPosition -0.5*bladeWidth_* mediumAxis);
        blade->addObject(-mediumAxis, anchorPosition +0.5*bladeWidth_* mediumAxis);
    }

    void write(std::ostream& os, bool x) const
    {
        Mercury3D::write(os, x);
        os  << " bladeInclination " << bladeInclination_ << std::endl
        << " bladeNumber " << bladeNumber_ << std::endl
        << " bladeHeight " << bladeHeight_ << std::endl
        << " bladeWidth " << bladeWidth_ << std::endl
        << " bladeThickness " << bladeThickness_ << std::endl;
    }

public:
    Mdouble bladeWidth_;
};

class BoosterWithStackedBlade : public BoosterWithSmallAngledBlade {
public:

    BoosterWithStackedBlade(Mdouble bladeThickness, Mdouble bladeHeight, Mdouble bladeDistance, Mdouble bladeWidth, Mdouble bladeInclination, unsigned int bladeNumber)
        : BoosterWithSmallAngledBlade(bladeThickness, bladeHeight, bladeWidth, bladeInclination, bladeNumber), bladeDistance_(bladeDistance)
    {
        setName("BoosterWithStackedBladeFixedBladeWidthAndInclination");
    }

    void addBlades()
    {
        //blade
        bladesDown.reserve(bladeNumber_);
        bladesRight.reserve(bladeNumber_);
        bladesUp.reserve(bladeNumber_);
        bladesLeft.reserve(bladeNumber_);
        stackBladesDown.reserve(bladeNumber_);
        stackBladesRight.reserve(bladeNumber_);
        stackBladesUp.reserve(bladeNumber_);
        stackBladesLeft.reserve(bladeNumber_);

        IntersectionOfWalls defaultBlade;
        defaultBlade.setSpecies(wallSpecies);

        Mdouble lidWidth = 216.73e-3;
        Mdouble upperWidth = 976.57e-3;
        Mdouble centreWidth = 2000.99e-3;

        double x = lidWidth+upperWidth;
        double dx = 0.5*centreWidth/static_cast<double>(bladeNumber_-0.5);
        for (unsigned int i=0; i<bladeNumber_; ++i)
        {
            //double x = (static_cast<double>(i)+0.5)*getXMax()/ static_cast<double>(bladeNumber_);
            bladesDown.push_back(wallHandler.copyAndAddObject(defaultBlade));
            bladesDown.back()->setPosition(Vec3D(x,0.0,0.0));
            bladesUp.push_back(wallHandler.copyAndAddObject(defaultBlade));
            bladesUp.back()->setPosition(Vec3D(x,0.0,0.0));
            stackBladesDown.push_back(wallHandler.copyAndAddObject(defaultBlade));
            stackBladesDown.back()->setPosition(Vec3D(x,0.0,0.0));
            stackBladesUp.push_back(wallHandler.copyAndAddObject(defaultBlade));
            stackBladesUp.back()->setPosition(Vec3D(x,0.0,0.0));
            x += dx;
            bladesRight.push_back(wallHandler.copyAndAddObject(defaultBlade));
            bladesRight.back()->setPosition(Vec3D(x,0.0,0.0));
            bladesLeft.push_back(wallHandler.copyAndAddObject(defaultBlade));
            bladesLeft.back()->setPosition(Vec3D(x,0.0,0.0));
            stackBladesRight.push_back(wallHandler.copyAndAddObject(defaultBlade));
            stackBladesRight.back()->setPosition(Vec3D(x,0.0,0.0));
            stackBladesLeft.push_back(wallHandler.copyAndAddObject(defaultBlade));
            stackBladesLeft.back()->setPosition(Vec3D(x,0.0,0.0));
            x += dx;
        }

        actionsAfterTimeStep();
    }
    
    void setBladeAngle(IntersectionOfWalls* blade, Mdouble angle)
    {
        blade->clear();
        Mdouble s = std::sin(angle), c = std::cos(angle);
        Vec3D anchorPosition = Vec3D(blade->getPosition().X, getZMax()*s, getZMax()*c); //clockwise as seen from lid
        Vec3D longAxis = Vec3D(0.0, s, c);
        Vec3D shortAxis = Vec3D(1./std::tan(bladeInclination_), -c, s); shortAxis.normalize();
        Vec3D mediumAxis = Vec3D::cross(shortAxis,longAxis);
        blade->addObject( longAxis, anchorPosition -bladeHeight_* longAxis);
        //blade->addObject(-longAxis, anchorPosition +bladeHeight_* longAxis);
        blade->addObject( shortAxis, anchorPosition -0.5*bladeThickness_* shortAxis);
        blade->addObject(-shortAxis, anchorPosition +0.5*bladeThickness_* shortAxis);
        blade->addObject( mediumAxis, anchorPosition -0.5*bladeWidth_* mediumAxis);
        blade->addObject(-mediumAxis, anchorPosition +0.5*bladeWidth_* mediumAxis);
    }
 
    void setStackBladeAngle(IntersectionOfWalls* blade, Mdouble angle)
    {
        blade->clear();
        Mdouble s = std::sin(angle), c = std::cos(angle);
        Vec3D anchorPosition = Vec3D(blade->getPosition().X, getZMax()*s, getZMax()*c); //clockwise as seen from lid
        Vec3D longAxis = Vec3D(0.0, s, c);
        Vec3D shortAxis = Vec3D(-1./std::tan(bladeInclination_), -c, s); shortAxis.normalize();
        Vec3D mediumAxis = Vec3D::cross(shortAxis,longAxis);
        blade->addObject( longAxis, anchorPosition -(2.0*bladeHeight_+bladeDistance_)* longAxis);
        blade->addObject(-longAxis, anchorPosition -(    bladeHeight_+bladeDistance_)* longAxis);
        //blade->addObject(longAxis, anchorPosition  + (     bladeHeight_)* longAxis);
        //blade->addObject(-longAxis, anchorPosition  + (2.0*bladeHeight_)* longAxis);
        blade->addObject( shortAxis, anchorPosition -0.5*bladeThickness_* shortAxis);
        blade->addObject(-shortAxis, anchorPosition +0.5*bladeThickness_* shortAxis);
        blade->addObject( mediumAxis, anchorPosition -0.5*bladeWidth_* mediumAxis);
        blade->addObject(-mediumAxis, anchorPosition +0.5*bladeWidth_* mediumAxis);
    }
    
    void actionsAfterTimeStep()
    {
        Mdouble angle = revolutionsPerSecond * getTime() * 2.0 * constants::pi;
        for (IntersectionOfWalls* blade : bladesDown)
            setBladeAngle(blade,angle);
        for (IntersectionOfWalls* blade : bladesRight)
            setBladeAngle(blade,angle+constants::pi/2.0);
        for (IntersectionOfWalls* blade : bladesUp)
            setBladeAngle(blade,angle+2.0*constants::pi/2.0);
        for (IntersectionOfWalls* blade : bladesLeft)
            setBladeAngle(blade,angle+3.0*constants::pi/2.0);
        for (IntersectionOfWalls* blade : stackBladesDown)
            setStackBladeAngle(blade,angle);
        for (IntersectionOfWalls* blade : stackBladesRight)
            setStackBladeAngle(blade,angle+constants::pi/2.0);
        for (IntersectionOfWalls* blade : stackBladesUp)
            setStackBladeAngle(blade,angle+2.0*constants::pi/2.0);
        for (IntersectionOfWalls* blade : stackBladesLeft)
            setStackBladeAngle(blade,angle+3.0*constants::pi/2.0);
    }
    
    

private:
    std::vector<IntersectionOfWalls*> stackBladesDown, stackBladesRight, stackBladesUp, stackBladesLeft;
    Mdouble bladeDistance_;
};

class BoosterWithOuterCoil : public Booster {
public:

    BoosterWithOuterCoil(Mdouble coilLength, Mdouble coilRadius, Mdouble coilAngularVelocity, Mdouble coilWindings)
        : Booster(), coilLength_(coilLength), coilRadius_(coilRadius), coilAngularVelocity_(coilAngularVelocity), coilWindings_(coilWindings)
    {
        setName("BoosterWithCoil");
        coilThickness_ = 7.2e-2*2.0;
        coilOrientation_ = Vec3D(1, 0, 0);
    }

    void addBlades()
    {
        std::cout << "creating outer coil" << std::endl;
        outerCoil = wallHandler.copyAndAddObject(Coil());
        setBladeAngle(0.0);
    }

    virtual void setBladeAngle(Mdouble angle)
    {
        Vec3D coilTangent = Vec3D(0, std::cos(angle), std::sin(angle));
        outerCoil->setPosition(Vec3D(getXMax() - coilLength_, 0, 0));
        outerCoil->set(coilOrientation_, coilTangent, coilLength_, coilRadius_, coilWindings_, coilThickness_);
    }

    void actionsAfterTimeStep()
    {
        Mdouble angle = revolutionsPerSecond * getTime() * 2.0 * constants::pi;
        setBladeAngle(angle);
    }
   

public:
    Coil* outerCoil;
    Mdouble coilLength_;
    Mdouble coilRadius_;
    Mdouble coilAngularVelocity_;
    Mdouble coilWindings_;
    Mdouble coilThickness_;
    Vec3D coilOrientation_;
};

class BoosterWithDoubleCoil : public BoosterWithOuterCoil {
public:

    BoosterWithDoubleCoil(Mdouble coilLength, Mdouble coilRadius, Mdouble coilAngularVelocity, Mdouble coilWindings, Mdouble innerCoilRadius)
        : BoosterWithOuterCoil(coilLength, coilRadius, coilAngularVelocity, coilWindings), innerCoilRadius_(innerCoilRadius)
    {
        setName("BoosterWithDoubleCoil");
    }

    void addBlades()
    {
        BoosterWithOuterCoil::addBlades();

        std::cout << "creating inner coil" << std::endl;
        innerCoil = wallHandler.copyAndAddObject(Coil());
        Vec3D coilTangent = Vec3D(0, 1, 0);
        innerCoil->setPosition(Vec3D(getXMax(), 0, 0));
        innerCoil->set(-coilOrientation_, coilTangent, coilLength_, innerCoilRadius_, coilWindings_, coilThickness_);
    }

public:
    Coil* innerCoil;
    Mdouble innerCoilRadius_;
};
