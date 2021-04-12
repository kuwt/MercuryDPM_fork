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

#include "Mercury3D.h"
#include "Species/SinterFrictionSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
using constants::pi;
using helpers::readFromFile;

/// One particle, sintering slowly to a wall
template <class SpeciesType>
class InitialConditions : public Mercury3D
{
    typedef typename SpeciesType::InteractionType InteractionType;


public:
    /**
     * Default constructor; creates (empty) species, particle, and geometry sets output details
     * @return
     */
    InitialConditions(Mdouble radius, Mdouble domainRadius, Mdouble domainHeight, unsigned dimensions) {
        //create empty species
        species = speciesHandler.copyAndAddObject(SpeciesType());

        //output settings
        setName("SinterBed0");
        //dataFile.setFileType(FileType::NO_FILE);
        fStatFile.setFileType(FileType::NO_FILE);
        //restartFile.setFileType(FileType::NO_FILE);
        setXBallsAdditionalArguments(" -v0 -solidf ");
        setSaveCount(1000);
        setTimeMax(0);
        restartFile.setSaveCount(restartFile.getSaveCount()*10);

        //set gravity and density
        setGravity(Vec3D(0, 0, -9.8));
        species->setDensity(1000); //overly heavy

        //define domain
        setXMin(-domainRadius);
        setYMin(-domainRadius);
        setZMin(0);
        setXMax(domainRadius);
        setYMax(domainRadius);
        setZMax(domainHeight);

        // Inserting wall
        InfiniteWall baseWall;
        baseWall.set(Vec3D(0.0,0.0,-1.0), Vec3D(0.0,0.0,0.0));
        baseWall.setSpecies(species);//!
        wallHandler.copyAndAddObject(baseWall);

//        AxisymmetricIntersectionOfWalls sideWall;
//        sideWall.setPosition(Vec3D(0, 0, 0));
//        sideWall.setOrientation(Vec3D(0, 0, 1));
//        sideWall.addObject(Vec3D(1,0,0),Vec3D(domainRadius,0,0));
//        sideWall.setSpecies(species);//!
//        wallHandler.copyAndAddObject(sideWall);

        PeriodicBoundary b;
        b.set(Vec3D(1,0,0),getXMin(),getXMax());
        boundaryHandler.copyAndAddObject(b);
        b.set(Vec3D(0,1,0),getYMin(),getYMax());
        boundaryHandler.copyAndAddObject(b);

        logger(INFO,"dimensions %",dimensions);

        // Inserting particles with 10% standard deviation
        hGridRebuild();
        SphericalParticle p;
        p.setRadius(1.1*radius);
        p.setSpecies(species);
        Mdouble volumeOfParticles = 0;
        Mdouble volumeOfParticlesMax;
        if (dimensions==3)
            volumeOfParticlesMax = 0.58*pi*domainRadius*domainRadius*domainHeight;
        else if (dimensions==2)
            volumeOfParticlesMax = 0.6*2.0*domainRadius*domainHeight;
        else
            volumeOfParticlesMax = domainHeight;

        Mdouble maxHeight=getZMin()+p.getRadius();
        while (volumeOfParticles < volumeOfParticlesMax)
        {
            Vec3D pos;
            pos.Z = random.getRandomNumber(getZMin()+p.getRadius(), maxHeight);
            if (dimensions>1)
                pos.X = random.getRandomNumber(-domainRadius+p.getRadius(), domainRadius-p.getRadius());
            if (dimensions>2)
                pos.Y = random.getRandomNumber(-domainRadius+p.getRadius(), domainRadius-p.getRadius());
            p.setPosition(pos);

            if (checkParticleForInteraction(p))
            {
                particleHandler.copyAndAddObject(p);
                if (dimensions==3)
                    volumeOfParticles += p.getVolume();
                else if (dimensions==2)
                    volumeOfParticles += pi*p.getRadius()*p.getRadius();
                else if (dimensions==1)
                    volumeOfParticles += p.getRadius()+p.getRadius();
                else
                    break;
                p.setRadius(random.getRandomNumber(0.9,1.1)*radius);
                if (particleHandler.getNumberOfObjects()%100==0) std::cout << "." << std::flush;
            } else {
                maxHeight+=.002*radius;
            }
        }
        std::cout << std::endl;
        logger(INFO,"Inserted % particles",particleHandler.getNumberOfObjects());

    }

    /** set time step */
    void setupInitialConditions() override
    {
//        //set species of outer wall to non-sintering
//        auto wallSpecies = speciesHandler.copyAndAddObject(species);
//        wallHandler.getLastObject()->setSpecies(wallSpecies);

        //set time step
        Mdouble mass = species->getSmallestParticleMass();
        setTimeStep(4.0*species->computeTimeStep(mass));
        Mdouble oldGravity = getGravity().Z;
        setGravity(getGravity()*1e5);//set to 5e5/(H/d)
        logger(INFO,"Time step: %", getTimeStep());
        logger(INFO,"Fudge factors: Gravity scaled by factor %, timestep scaled by %",getGravity().Z/oldGravity,getTimeStep()/species->computeTimeStep(mass));
    }

    //override continueSolve function such that the code stops when the packing is relaxed (Ekin<1e-5*Eela)
    bool continueSolve() const override
    {
        static unsigned int counter = 0;
        if (++counter>100)
        {
            counter=0;
            if (getKineticEnergy()<3e-3*getElasticEnergy())
                return false;
        }
        return true;
    }

    Mdouble getMeanRelativeContactRadius() const
    {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanRadius = particleHandler.getMeanRadius();
        //logger(INFO,"d % r %",meanPlasticOverlap,meanRadius);
        return sqrt(2.0*meanOverlap/meanRadius);
    }

    //override printTime function such that console output shows the state of relaxation (Ekin/Eela)
    void printTime() const override
    {
        std::cout << "t=" << std::setw(12) << getTime()
                  << " Ene " << std::setw(12) << getKineticEnergy()/getElasticEnergy()
                  << " X/a " << std::setw(12) << getMeanRelativeContactRadius()
                  << std::endl;
    }

//    //extra gravity force to make the explosion effect less
//    void computeExternalForces(BaseParticle* p) override
//    {
//        DPMBase::computeExternalForces(p);
//        if (getTime()>10.0)
//        {
//            p->setVelocity(0.99*p->getVelocity());
//        }
//    }

    void scaleMass(Mdouble scale)
    {
        logger(INFO,"Scaling density, inverse gravity, square of timestep, and square of diss. coeff. up by a factor %",scale);
        species->setDensity(scale * species->getDensity());
        setGravity(getGravity()/scale);
        //setTimeStep(sqrt(scale) * getTimeStep());
        species->setDissipation(sqrt(scale) * species->getDissipation());
    }

    SpeciesType* species;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    const Mdouble radius = readFromFile("in","radius",6e-6);
    const Mdouble domainRadius = readFromFile("in","domainDiameterOverParticleDiameter",5.0) * radius;
    const Mdouble domainHeight = readFromFile("in","domainHeightOverParticleDiameter",5.0) * 2.0 * radius;
    const unsigned dimensions = static_cast<const unsigned int>(readFromFile("in", "dimensions", 2));

    InitialConditions<SinterFrictionSpecies> ic(radius, domainRadius, domainHeight,dimensions);
    const unsigned long int random = readFromFile("in","random",0);
    ic.random.setRandomSeed(random);
    ic.scaleMass(1e12);
    //set species properties
    const Mdouble kMax=readFromFile("in","kMax",20000);
    const Mdouble kL = readFromFile("in","kL",0.1*kMax);
    ic.species->setPlasticParameters(kL, kMax, kL, 0.16);
    const Mdouble mass = ic.species->getSmallestParticleMass();
    ic.species->setStiffnessAndRestitutionCoefficient(ic.species->getLoadingStiffness(),0.5,mass);
    ic.species->setSinterType(SINTERTYPE::CONSTANT_RATE);
    ic.species->setSinterRate(0.0);
    ic.species->setSinterAdhesion(1200e-9/radius);
    ic.species->setSlidingDissipation(2./7.*ic.species->getDissipation());
    ic.species->setSlidingStiffness(2./7.*ic.species->getLoadingStiffness());
    ic.species->setRollingDissipation(2./7.*ic.species->getDissipation());
    ic.species->setRollingStiffness(2./7.*ic.species->getLoadingStiffness());
    ic.species->setTorsionDissipation(2./7.*ic.species->getDissipation());
    ic.species->setTorsionStiffness(2./7.*ic.species->getLoadingStiffness());

    ic.setTimeMax(1e20);
//    ic.interactionHandler.setWriteVTK(FileType::MULTIPLE_FILES);
//    ic.setParticlesWriteVTK(true);
//    ic.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    ic.solve();

    return 0;
    ///\todo TW automatically measure simulation time
    //r=3 -> 220 particles
}
