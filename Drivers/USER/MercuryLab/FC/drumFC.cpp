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

/** This code creates a cylindrical container, inserts particles and lets them settle.
*/
class Drum : public Mercury3D{
public:

    //set default values
    Drum()
    {
        setName("drum");
        setFileType(FileType::ONE_FILE);
        setSaveCount(200);
        setTimeMax(1e20); //run forever
        rotationsPerSecond = 0.0;
        inclination = 0.0;

        //creating species for particle-particle interactions and particle-wall interactions
        particleSpecies = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        wallSpecies = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        wallParticleSpecies = speciesHandler.getMixedObject(wallSpecies, particleSpecies);

        collisionTime = std::numeric_limits<double>::quiet_NaN();
        restitutionCoefficient = std::numeric_limits<double>::quiet_NaN();
    }

    //set slave variables
    void setupInitialConditions()
	{
        setGravity(9.8*Vec3D(std::sin(inclination),0.0,-std::cos(inclination)));

        //set contact properties:
        //- calculate stiffness and dissipation
        Mdouble effectiveMass = 0.5*particleSpecies->getMassFromRadius(0.5*particleDiameter);
        std::cout << "Mass" << 2.0*effectiveMass << std::endl;
        //helpers::KAndDisp kAndDisp = helpers::computeKAndDispFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(collisionTime, restitutionCoefficient, effectiveMass);
        particleSpecies->setCollisionTimeAndRestitutionCoefficient(collisionTime, restitutionCoefficient, 2.0*effectiveMass);
        //- set stiffness, dissipation and timestep
        //particleSpecies->setStiffness(kAndDisp.k);
        //particleSpecies->setDissipation(kAndDisp.disp);
        particleSpecies->setSlidingStiffness(2.0/7.0*particleSpecies->getStiffness());
        particleSpecies->setSlidingDissipation(2.0/7.0*particleSpecies->getDissipation());
        particleSpecies->setSlidingFrictionCoefficient(0.5);
        setTimeStep(0.02*collisionTime);

        //set wall properties
        wallParticleSpecies->setStiffness(particleSpecies->getStiffness());
        wallParticleSpecies->setDissipation(particleSpecies->getDissipation());
        wallParticleSpecies->setSlidingStiffness(particleSpecies->getSlidingStiffness());
        wallParticleSpecies->setSlidingDissipation(particleSpecies->getSlidingDissipation());
        wallParticleSpecies->setSlidingFrictionCoefficient(particleSpecies->getSlidingFrictionCoefficient());

        Mdouble lidWidth = 216.73e-3;
        Mdouble upperWidth = 976.57e-3;
        Mdouble centerWidth = 2000.99e-3;
        Mdouble bottomWidth = 322.68e-3;
        Mdouble lidRadius = 300e-3;
        Mdouble centerRadius = 642.5e-3;
        Mdouble bottomRadius = 120e-3;

        setXMax(lidWidth+upperWidth+centerWidth+bottomWidth);
        setXMin(0.0);
        setYMax(centerRadius);
        setYMin(-getYMax());
        setZMax(getYMax());
        setZMin(getYMin());

        // base wall
        baseWall = wallHandler.copyAndAddObject(InfiniteWall());
        baseWall->set(Vec3D(1.0,0.0,0.0), Vec3D(getXMax(),0.0,0.0) );
        baseWall->setSpecies(wallSpecies);
        baseWall->setAngularVelocity(Vec3D(rotationsPerSecond * 2.0 * constants::pi,0.0,0.0));

        //center wall
        sideWall = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
        sideWall->setPosition(Vec3D(0.0,0.0,0.0) );
        sideWall->setOrientation(Vec3D(1.0,0.0,0.0) );
        sideWall->addObject(Vec3D(1,0,0), Vec3D(centerRadius,0.0,0.0) );
        sideWall->setSpecies(wallSpecies);
        sideWall->setAngularVelocity(Vec3D(rotationsPerSecond * 2.0 * constants::pi,0.0,0.0));

        //top wall
        upperWall = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
        upperWall->setPosition(Vec3D(0.0,0.0,0.0) );
        upperWall->setOrientation(Vec3D(1.0,0.0,0.0) );
        upperWall->addObject(Vec3D(upperWidth-lidWidth,0,lidRadius-centerRadius), Vec3D(lidRadius,0.0,lidWidth) );
        upperWall->setSpecies(wallSpecies);
        upperWall->setAngularVelocity(Vec3D(rotationsPerSecond * 2.0 * constants::pi,0.0,0.0));

        //lower wall
        lowerWall = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
//        lowerWall = new AxisymmetricIntersectionOfWalls();
        lowerWall->setPosition(Vec3D(0.0,0.0,0.0) );
        lowerWall->setOrientation(Vec3D(1.0,0.0,0.0) );
        //the 30 is because the bending of the lower wall only starts after 30 mm (see sketch)
        ///\todo tw Walls still require that the normal is of unit length!!
        Vec3D normal = Vec3D(bottomWidth-30e-3,0,centerRadius-bottomRadius); normal.normalize();
        lowerWall->addObject(normal, Vec3D(bottomRadius,0.0,getXMax()) );
        lowerWall->setSpecies(wallSpecies);
        lowerWall->setAngularVelocity(Vec3D(rotationsPerSecond * 2.0 * constants::pi,0.0,0.0));

        addBlades();

        // check if all particle diameter is set correctly
        if (particleDiameter==0.0)
        {
            std::cerr << "Error in Drum::setupDrum: The particle diameter is not set." << std::endl;
            exit(-1);
        }

        // particles are set based on particle diameter
        BaseParticle P;
        P.setSpecies(particleSpecies);

        Vec3D pos;
        Mdouble numberOfParticlesInserted = 0;
        Mdouble numberOfParticles = 0.4*(getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin())/mathsFunc::cubic(particleDiameter)/(4.0/3.14);
        hGridRebuild();
        std::cout << std::endl << "Status before inserting particles:" << std::endl;
        write(std::cout,false);
        std::cout << std::endl << "Inserting " << std::floor(numberOfParticles) << " particles" << std::endl;
        P.setRadius(0.5*particleDiameter*random.getRandomNumber(0.95,1.05));
        while (numberOfParticlesInserted < numberOfParticles)
        {
            pos.X = random.getRandomNumber(getXMin()+P.getRadius(), getXMax()-P.getRadius());
            pos.Y = random.getRandomNumber(getYMin()+P.getRadius(), getYMax()-P.getRadius());
            pos.Z = random.getRandomNumber(getZMin()+P.getRadius(), getZMax()-P.getRadius());
            P.setPosition(pos);

            //std::cout << checkParticleForInteraction(P);
            if (checkParticleForInteraction(P))
            {
                particleHandler.copyAndAddObject(P);
                P.setRadius(0.5*particleDiameter*random.getRandomNumber(0.95,1.05));
                ++numberOfParticlesInserted;
                std::cout << "." << std::flush;
            }
        }
        std::cout << std::endl;
        std::cout << std::endl << "Starting simulation, timeMax= " << getTimeMax() << std::endl;
    }

    virtual void addBlades()
    {   }

    void printTime() const
    {
        std::cout << "t=" << getTime() << " Ene " << getKineticEnergy()/getElasticEnergy() << std::endl;
    }

    LinearViscoelasticSlidingFrictionSpecies* particleSpecies;
    LinearViscoelasticSlidingFrictionMixedSpecies* wallParticleSpecies;
    LinearViscoelasticSlidingFrictionSpecies* wallSpecies;
    AxisymmetricIntersectionOfWalls* upperWall, * sideWall, *lowerWall;
    InfiniteWall* baseWall;

    double particleDiameter;
    double inclination;
    double rotationsPerSecond;
    Mdouble collisionTime; //softness of particles
    Mdouble restitutionCoefficient; // dissipativeness of particles
};

class DrumWithStraightBlade : public Drum {
public:

    DrumWithStraightBlade(Mdouble bladeThickness, Mdouble bladeHeight)
    : Drum(), bladeThickness_(bladeThickness), bladeHeight_(bladeHeight)
    {   }

    void addBlades()
    {
        //blade
        blade0 = wallHandler.copyAndAddObject(IntersectionOfWalls());
        blade1 = wallHandler.copyAndAddObject(IntersectionOfWalls());
        blade2 = wallHandler.copyAndAddObject(IntersectionOfWalls());
        blade3 = wallHandler.copyAndAddObject(IntersectionOfWalls());
        blade0->setSpecies(wallSpecies);
        blade1->setSpecies(wallSpecies);
        blade2->setSpecies(wallSpecies);
        blade3->setSpecies(wallSpecies);
        setBladeAngle(blade0,0.0);
        setBladeAngle(blade1,constants::pi/2.0);
        setBladeAngle(blade2,2.0*constants::pi/2.0);
        setBladeAngle(blade3,3.0*constants::pi/2.0);
    }

    virtual void setBladeAngle(IntersectionOfWalls* blade, Mdouble angle)
    {
        blade->clear();
        Vec3D longAxis = Vec3D(0.0, std::sin(angle), -std::cos(angle));
        Vec3D shortAxis = Vec3D(0.0, -longAxis.Z, longAxis.Y);
        blade->addObject(  longAxis, (getZMax()-bladeHeight_)*longAxis);
        blade->addObject( shortAxis, -0.5*bladeThickness_*shortAxis);
        blade->addObject(-shortAxis,  0.5*bladeThickness_*shortAxis);
        sideWall->setOrientation(Vec3D(1.0,0.0,0.0) );
        upperWall->setOrientation(Vec3D(1.0,0.0,0.0) );
        lowerWall->setOrientation(Vec3D(1.0,0.0,0.0) );
    }

    void actionsAfterTimeStep()
    {
        Mdouble angle = rotationsPerSecond * getTime() * 2.0 * constants::pi;
        setBladeAngle(blade0,angle);
        setBladeAngle(blade1,angle+constants::pi/2.0);
        setBladeAngle(blade2,angle+2.0*constants::pi/2.0);
        setBladeAngle(blade3,angle+3.0*constants::pi/2.0);
    }

protected:
    IntersectionOfWalls* blade0, * blade1, * blade2, * blade3;
    Mdouble bladeThickness_;
    Mdouble bladeHeight_;
};

class DrumWithAngledBlade : public DrumWithStraightBlade {
public:

    DrumWithAngledBlade(Mdouble bladeThickness, Mdouble bladeHeight, Mdouble bladeInclination)
    : DrumWithStraightBlade(bladeThickness,bladeHeight), bladeInclination_(bladeInclination)
    {   }

    void setBladeAngle(IntersectionOfWalls* blade, Mdouble angle)
    {
        blade->clear();
        Vec3D longAxis = Vec3D(0.0, std::sin(angle), -std::cos(angle));
        Vec3D shortAxis = Vec3D(std::tan(bladeInclination_), -longAxis.Z, longAxis.Y);
        shortAxis.normalize();
        blade->addObject(  longAxis, (getZMax()-bladeHeight_)*longAxis);                                  //  _0_  YZ-plane
        blade->addObject( shortAxis, -0.5*bladeThickness_*shortAxis+Vec3D(0.5*getXMax(),0,0)); // |   |
        blade->addObject(-shortAxis,  0.5*bladeThickness_*shortAxis+Vec3D(0.5*getXMax(),0,0)); // |2  |1
        sideWall->setOrientation(Vec3D(1.0,0.0,0.0) );
        upperWall->setOrientation(Vec3D(1.0,0.0,0.0) );
        lowerWall->setOrientation(Vec3D(1.0,0.0,0.0) );
    }

private:
    // inclination of the blade in the xy plane when the blade is at it's lowest point (when setBladeAngle::angle=0)
    Mdouble bladeInclination_;
};

class DrumWithDoubleAngledBlade : public DrumWithStraightBlade {
public:

    DrumWithDoubleAngledBlade(Mdouble bladeThickness, Mdouble bladeHeight, Mdouble verticalBladeInclination, Mdouble horizontalBladeInclination)
        : DrumWithStraightBlade(bladeThickness,bladeHeight), verticalBladeInclination_(verticalBladeInclination), horizontalBladeInclination_(horizontalBladeInclination)
    {   }

    void setBladeAngle(IntersectionOfWalls* blade, Mdouble angle)
    {
        blade->clear();
        Vec3D anchorPosition = Vec3D(0.5*getXMax(), getZMax()*std::sin(angle), -getZMax()*std::cos(angle));
        Vec3D longAxis = Vec3D(0.0, std::sin(angle-verticalBladeInclination_), -std::cos(angle-verticalBladeInclination_));
        Vec3D shortAxis = Vec3D(std::tan(horizontalBladeInclination_), -longAxis.Z, longAxis.Y);
        shortAxis.normalize();
        blade->addObject(longAxis, anchorPosition -bladeHeight_* longAxis);
        blade->addObject( shortAxis, anchorPosition -0.5*bladeThickness_* shortAxis);
        blade->addObject(-shortAxis, anchorPosition +0.5*bladeThickness_* shortAxis);
        sideWall->setOrientation(Vec3D(1.0,0.0,0.0) );
        upperWall->setOrientation(Vec3D(1.0,0.0,0.0) );
        lowerWall->setOrientation(Vec3D(1.0,0.0,0.0) );
    }


private:
    // inclination of the blade in the xz plane when the blade is at it's lowest point (when setBladeAngle::angle=0)
    Mdouble verticalBladeInclination_;
    // inclination of the blade in the xy plane when the blade is at it's lowest point (when setBladeAngle::angle=0)
    Mdouble horizontalBladeInclination_;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{

    //Uncomment this (and comment the other constructors) if you want a straight blade
    Mdouble bladeThickness = 0.04*642.5e-3; // default: 4% of the center radius of 642.5e-3 m
    Mdouble bladeHeight = 0.1*642.5e-3; // default: 10% of the center radius of 642.5e-3 m
    DrumWithStraightBlade drum(bladeThickness,bladeHeight);

//    //Uncomment this (and comment the other constructors) if you want an angled blade
//    Mdouble bladeInclination = 5.71059314*constants::pi/180.0; // default angle of 5.71059314 deg such that tan(angle)=0.1 (it's easy to check the normals then)
//    Mdouble bladeThickness = 0.2*642.5e-3;
//    Mdouble bladeHeight = 0.5*642.5e-3;
//	DrumWithAngledBlade drum(bladeThickness, bladeHeight, bladeInclination);

//    //Uncomment this (and comment the other constructors) if you want a double angled blade
//    Mdouble verticalBladeInclination = 5.71059314*constants::pi/180.0;
//    Mdouble horizontalBladeInclination = 0.0*constants::pi/180.0; // default angle of 5.71059314 deg such that tan(angle)=0.1 (it's easy to check the normals then)
//    Mdouble bladeThickness = 0.1*642.5e-3;
//    Mdouble bladeHeight = 0.2*642.5e-3;
//	DrumWithDoubleAngledBlade drum(bladeThickness, bladeHeight, verticalBladeInclination, horizontalBladeInclination);

    drum.inclination = 5.0*constants::pi/180.0;
    drum.rotationsPerSecond = 0.1;
    drum.particleDiameter = 7.2e-2;
    drum.particleSpecies->setDensity(2000.0);
    drum.collisionTime = 0.01; //soft particles
    drum.restitutionCoefficient = 0.001; // dissipative particles
    drum.setXBallsAdditionalArguments(" -v0 -solidf  -w 1250 -s 0.76 -noborder 3");
    drum.solve();
    return true;
}
