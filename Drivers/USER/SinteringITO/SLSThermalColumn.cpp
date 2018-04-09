// This file is part of MercuryDPM.
// 
// MercuryDPM is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// MercuryDPM is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with MercuryDPM.  If not, see <http://www.gnu.org/licenses/>.
// 
// Copyright 2016 The Mercury Developers Team
// For the list of developers, see <http://www.MercuryDPM.org/Team>


//This code does a column with a given number of particles. After 1 second as the KE is low enough it sets the top partilce to a constant tempture and sees how to propagates.

#include "Mercury3D.h"
#include <Particles/ThermalParticle.h>
#include <Species/NormalForceSpecies/ThermalSpecies.h>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include <Species/HertzianViscoelasticFrictionSpecies.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Walls/InfiniteWall.h>
#include <map>
#include <random>
using constants::pi;
using mathsFunc::cubic;

typedef Species<ThermalSpecies<LinearViscoelasticNormalSpecies>,SlidingFrictionSpecies> ThermalLinearViscoelasticSlidingFrictionSpecies;
typedef MixedSpecies<ThermalSpecies<LinearViscoelasticNormalSpecies>,SlidingFrictionSpecies> ThermalLinearViscoelasticSlidingFrictionMixedSpecies;


/**
 * Simulates the spread of heat in an SLS powder bed.
 */
#include "Mercury3D.h"
#include "Species/HertzianViscoelasticFrictionSpecies.h"

class PowderBed : public Mercury3D
{
public:

    PowderBed (Mdouble domainLength, Mdouble domainWidth, Mdouble domainDepth, ParticleSpecies& species, Mdouble initialTemperature, Mdouble finalTemperature)
    : initialTemperature_(initialTemperature), finalTemperature_(finalTemperature)
    {
        
        TemperatureFile_.open("Temperature.data");
        //set name, gravity, output options
        setName("SLSThermalColumn");
        setGravity({0,0,-9.8});
        setXBallsAdditionalArguments("-solidf -v0 -cmode 8 -cmaxset 100 ");
        setTimeMax(10.0);
        logger(INFO,"Name of output file: %",getName());

        logger(INFO,"Defining domain as [-l/2,l/2]x[-w/2,w/2]x[0,d] with l=%, w=%, d=%",domainLength, domainWidth,domainDepth);
        setMin({-0.5*domainLength,-0.5*domainWidth,0});
        setMax({0.5*domainLength,0.5*domainWidth,domainDepth});

        logger(INFO,"Adding contact law");
        speciesHandler.copyAndAddObject(species); //particle-paricle interactions

        logger(INFO,"Setting up periodic boundary conditions in x and y");
        PeriodicBoundary p;
        p.set({1,0,0},getXMin(),getXMax());
        boundaryHandler.copyAndAddObject(p);
        p.set({0,1,0},getYMin(),getYMax());
        boundaryHandler.copyAndAddObject(p);

        logger(INFO,"Creating a base wall at z=0");
        InfiniteWall w;
        w.setSpecies(speciesHandler.getObject(0));
        w.set({0,0,-1},{0,0,getZMin()});
        wallHandler.copyAndAddObject(w);
    }

    void setGaussianDistribution (Mdouble meanRadius, Mdouble stdRadius) {
        //set up random number generator
        std::random_device rd;
        std::mt19937 gen(rd());

        //set up normal distribution
        std::normal_distribution<> d(meanRadius, stdRadius);

        //the volume to be added (assumes a volume fraction of 0.65)
        Mdouble addVolume = 0.65*(getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin());

        //add particles until the volume to be added is zero
        logger(INFO,"Adding particles ...");
        ThermalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(meanRadius);
        p.setTemperature(initialTemperature_);
        Mdouble fillHeight = getZMin();
        while (addVolume>0) {
            Mdouble x = random.getRandomNumber(getXMin(), getXMax());
            Mdouble y = random.getRandomNumber(getYMin(), getYMax());
            Mdouble z = random.getRandomNumber(getZMin(), fillHeight);
            p.setPosition({x, y, z});
            // check if particle can be inserted
            if (checkParticleForInteraction(p)) {
                particleHandler.copyAndAddObject(p);
                addVolume -= p.getVolume();
                do {
                    p.setRadius(d(gen));
                } while (p.getRadius()<0.1*meanRadius || p.getRadius()>0.9*meanRadius); //reject too small or large radii
                if (particleHandler.getNumberOfObjects()%100==0) std::cout << '.' << std::flush;
                if (particleHandler.getNumberOfObjects()%1000==0) std::cout << ' ';
                if (particleHandler.getNumberOfObjects()%10000==0) std::cout << addVolume << '\n';
            } else {
                fillHeight += 0.01*meanRadius; //increase fill height (slowly to insert particles as low as possible)
            }
        }
        logger(INFO," Inserted % particles",particleHandler.getNumberOfObjects());
    }

    void createThreeParticleStack (Mdouble radius1, Mdouble radius2, Mdouble radius3) {
        ThermalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(radius1);
        p.setPosition({0,0,radius1});
        p.setTemperature(0);
        particleHandler.copyAndAddObject(p);
        p.setRadius(radius2);
        p.setPosition({0,0,radius1+radius1+radius2});
        p.setTemperature(1);
        particleHandler.copyAndAddObject(p);
        p.setRadius(radius2);
        p.setPosition({0,0,radius1+radius1+radius2+radius2+radius3});
        p.setTemperature(2);
        particleHandler.copyAndAddObject(p);
        logger(INFO," Inserted % particles",particleHandler.getNumberOfObjects());

    }
    
    void createColumn (unsigned int numParticles, Mdouble radius, Mdouble temperatureBottom, Mdouble temperatureRest)
    {
        ThermalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(radius);
        p.setPosition({0,0,radius});
        p.setTemperature(temperatureBottom);
        particleHandler.copyAndAddObject(p);
        
        for (int i=1; i<numParticles;i++)
        {
            p.setPosition({0,0,radius+2*i*radius});
            p.setTemperature(temperatureRest);
            particleHandler.copyAndAddObject(p);
            
        }
        
        logger(INFO," Inserted % particles",particleHandler.getNumberOfObjects());
        
    }

    void setUniformDistribution (Mdouble minRadius, Mdouble maxRadius) {
        //set up random number generator
        std::random_device rd;
        std::mt19937 gen(rd());

        //set up normal distribution
        std::uniform_real_distribution<> d(minRadius, maxRadius);

        //the volume to be added
        Mdouble addVolume = 0.65*(getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin());

        //add particles until the volume to be added is zero
        logger(INFO,"Adding particles ...");
        ThermalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        Mdouble meanRadius = 0.5*(minRadius + maxRadius);
        p.setRadius(meanRadius);
        p.setTemperature(initialTemperature_);
        Mdouble fillHeight = getZMin();
        while (addVolume>0) {
            Mdouble x = random.getRandomNumber(getXMin(), getXMax());
            Mdouble y = random.getRandomNumber(getYMin(), getYMax());
            Mdouble z = random.getRandomNumber(getZMin(), fillHeight);
            p.setPosition({x, y, z});
            // check if particle can be inserted
            if (checkParticleForInteraction(p)) {
                particleHandler.copyAndAddObject(p);
                addVolume -= p.getVolume();
                p.setRadius(d(gen));
                if (particleHandler.getNumberOfObjects()%100==0) std::cout << '.' << std::flush;
                if (particleHandler.getNumberOfObjects()%1000==0) std::cout << ' ';
                if (particleHandler.getNumberOfObjects()%10000==0) std::cout << addVolume << '\n';
            } else {
                fillHeight += 0.01*meanRadius; //increase fill height (slowly to insert particles as low as possible)
            }
        }
        logger(INFO," Inserted % particles",particleHandler.getNumberOfObjects());
    }

    //remove particles according to (max.) outflow rate
    void actionsAfterTimeStep() override
    {
        //make decisions based on settled state
        static ThermalParticle* heatedParticle = nullptr;
        
        if (getTime()<1.0) return;

        
        TemperatureFile_ << (getTime()-1.0) <<"\t";
        
        for (int i=0;i<particleHandler.getNumberOfObjects();i++)
        {
            double particleTemperture =dynamic_cast<const ThermalParticle*>(particleHandler.getObject(i))->getTemperature();
            TemperatureFile_ << particleTemperture << "\t";
            
        }
        TemperatureFile_ << "\n";
        
        if (heatedParticle == nullptr)
        {
            // Functions after this statement only get executed every 100th timestep (if counter==100)
            static unsigned counter = 0;
            if (++counter != 100) return;
            else counter = 0;
            
            
            if (getKineticEnergy() < 5e-10 )
            {
                //set heated particle
                Mdouble minDistance = std::numeric_limits<Mdouble>::infinity();
                for (auto p : particleHandler) {
                    Mdouble distance = Vec3D::getDistance(p->getPosition(),Vec3D(0,0,getZMax()));
                    if (distance < minDistance) {
                        minDistance = distance;
                        heatedParticle = dynamic_cast<ThermalParticle*>(p);
                    }
                }
                logger(INFO, "Starting to heat the material");
            }
        } else {
            heatedParticle->setTemperature(finalTemperature_);
        }
    }

    // display time and eneRatio every time the files are printed
    void printTime() const override
    {
        int p1 = 0;
        int p2 = 1;
        int p3 = 2;
        
        std::cout
         << "t " << std::setprecision(3) << std::left << std::setw(8)
         << getTime()
         << " EneRatio " << std::setprecision(3) << std::left << std::setw(8)
        << getKineticEnergy()/getElasticEnergy() <<std::endl
         << " T0 " << std::setprecision(3) << std::left << std::setw(6)
        << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(p1))->getTemperature()
         << "At Position:" <<particleHandler.getObject(p1)->getPosition()
        <<std::endl
         << " T1 " << std::setprecision(3) << std::left << std::setw(6)
         << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(p2))->getTemperature()
           << "At Position:" <<particleHandler.getObject(p2)->getPosition()
        << std::endl
         << " T2 " << std::setprecision(3) << std::left << std::setw(6)
         << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(p3))->getTemperature()
           << "At Position:" <<particleHandler.getObject(p3)->getPosition()
         << std::endl;
    }

    double getInfo(const BaseParticle& p) const
    {
        return (dynamic_cast<const ThermalParticle &>(p).getTemperature()-initialTemperature_)/(finalTemperature_-initialTemperature_);
        //return dynamic_cast<const ThermalParticle &>(p).getTemperature();
    }

private:

    Mdouble initialTemperature_;
    Mdouble finalTemperature_;
        
    std::ofstream TemperatureFile_;
};



int main(int argc UNUSED, char *argv[] UNUSED)
{
    //define domain size
    Mdouble domainLength = 0.25e-3;
    Mdouble domainWidth = 0.25e-3;
    Mdouble domainDepth = 1.0e-3;


    //define Temperatures
    Mdouble initialTemperature = 363; //K
    Mdouble finalTemperature = 1200; //K
    
    //define particles in column
    unsigned int numberParticlesInColumn=20;

    //define properties of particle-particle contacts
    Mdouble density = 1000;
    Mdouble collisionTime = 3e-4;
    Mdouble restitution = 0.2;
    Mdouble friction = 0.5;

    //put above properties into a contact law (do not modify)
    ThermalLinearViscoelasticSlidingFrictionSpecies species;
    Mdouble mass = density*4./3.*pi*cubic(25e-6);
    species.setDensity(density);
    species.setCollisionTimeAndRestitutionCoefficient(collisionTime, restitution, mass);
    species.setSlidingFrictionCoefficient(friction);
    species.setSlidingStiffness(2.0/7.0*species.getStiffness());
    species.setSlidingDissipation(2.0/7.0*species.getDissipation());
    //c specific heat cap
    species.setHeatCapacity(500); //J/kg/K
    //K therm cond
    species.setThermalConductivity(20); //W/m/K
    logger(INFO,"Gravitational compression %\% of radius per layer",100*mass*9.8/species.getStiffness()/25e-6);

    //Create a solver and run the commands in the constructor PowderBed::PowderBed
    PowderBed pb (domainLength, domainWidth, domainDepth, species, initialTemperature, finalTemperature);
    pb.setTimeStep(0.05*collisionTime);
    pb.setSaveCount(5.*collisionTime/pb.getTimeStep()); //save every 20 collisions (good for visual output)
    //pb.setSaveCount(10); //save every 20 collisions (good for visual output)

    //Create particles of a certain distribution
    Mdouble meanRadius = 25e-6;
    Mdouble stdRadius = 0.1 * meanRadius;
    //pb.setGaussianDistribution(meanRadius,stdRadius);
    //pb.setUniformDistribution(minRadius,maxRadius);
    //pb.createThreeParticleStack (1.0*meanRadius,1.0*meanRadius,1.0*meanRadius);
    //pb.createColumn(numberParticlesInColumn,meanRadius,finalTemperature,initialTemperature);
    pb.createColumn(numberParticlesInColumn,meanRadius,initialTemperature,initialTemperature);

    //run initial conditions
   // pb.setParticlesWriteVTK(true);
    //pb.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    pb.solve();

    return 0;
}
