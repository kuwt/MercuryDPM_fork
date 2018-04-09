//
// Created by irana on 3/2/18.
//

#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "Boundaries/PeriodicBoundary.h"
#include "PeriodicChute.h"
#include "Walls/InfiniteWall.h"

PeriodicChute::PeriodicChute(std::string roughBottomFile, Mdouble height, bool isMonodisperse)
{
    this->isMonoDisperse = isMonodisperse;
    if (isMonodisperse)
    {
        numberOfParticles = static_cast<unsigned int> (300 * height);
    }
    else
    {
        numberOfLarge = static_cast<unsigned int> (height / 10 * 633);
        numberOfSmall = 8 * numberOfLarge;
    }
    //read bottom particles from restart file
    readRestartFile(roughBottomFile);
    setRestarted(false);
    setTime(0);
    resetFileCounter();
    
    //clear all handlers except particle handler
    boundaryHandler.clear();
    wallHandler.clear();
    speciesHandler.clear();
    periodicBoundaryHandler.clear();
    domainHandler.clear();
    interactionHandler.clear();
    cgHandler.clear();
    
    //set all properties
    setName("PeriodicChute" + std::to_string((int) height) + "isMonoDisperse" + std::to_string(isMonodisperse));
    setSpeciesProperties();
    
    //species of fixed particles
    for (BaseParticle* p : particleHandler)
    {
        p->setSpecies(speciesHandler.getObject(0));
    }
    setChuteProperties();
    setTimeAndSaveCount();
    if (!isMonodisperse)
    {
        std::string filename = "COM_H" + std::to_string(height) + "_bidisperse";
        comFile.open(filename, std::ofstream::out);
        comFile << std::setw(12);
        comFile << "time" << std::setw(12) << "com_large" << std::setw(12)
                << "com_small" << std::endl;
    }
}

///compute the centre of mass of a) all flow particles b) all large particles c) all small particles
///also check if the flow is arrested.
///note, that we do not check every time step, but only every 1t.
void PeriodicChute::actionsAfterTimeStep()
{
    if (!isMonoDisperse)
    {
        static Mdouble nextCheckedTime = 1;
        if (getTime() > nextCheckedTime)
        {
            //compute COM of species 1
            Mdouble com1 = 0;
            int n1 = 0;
            //compute COM of species 2
            Mdouble com2 = 0;
            int n2 = 0;
            for (const BaseParticle* const p : particleHandler)
            {
                if (p->getSpecies()->getIndex() == 1)
                {
                    com1 += p->getPosition().Z;
                    n1++;
                }
                else if (p->getSpecies()->getIndex() == 2)
                {
                    com2 += p->getPosition().Z;
                    n2++;
                }
            }
            com1 /= n1;
            if (n2 > 0)
                com2 /= n2;
            nextCheckedTime += 1;
            logger(INFO, "COM species 1: %, COM species 2: %", com1, com2);
            comFile << std::setw(12);
            comFile << getTime() << std::setw(12) << com1 << std::setw(12) << com2 << std::endl;
            if (getKineticEnergy() < 1e-5)
                logger(ERROR, "The flow has arrested");
        }
    }
}

void PeriodicChute::setupInitialConditions()
{
    //Make sure that there are no gaps that particles can fall through
    InfiniteWall w0;
    w0.setSpecies(speciesHandler.getObject(0));
    w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, -1));
    wallHandler.copyAndAddObject(w0);
    
    
    //after the bottom is made, the boundaries can be put in place and the particles inserted
    setBoundaries();
    
    insertParticles();
    logger(INFO, "initial conditions are set");
}

void PeriodicChute::insertParticles()
{
    if (isMonoDisperse)
    {
        BaseParticle p0;
        p0.setSpecies(speciesHandler.getObject(1));
        for (unsigned int numberInserted = 0; numberInserted < numberOfParticles; ++numberInserted)
        {
            const Mdouble radius = random.getRandomNumber(0.45, 0.55);
            p0.setRadius(radius);
            insertOneParticle(p0);
        }
    }
    else
    {
        insertParticles(numberOfLarge, numberOfSmall);
    }
}

void PeriodicChute::insertParticles(unsigned int numberOfLargeToGo, unsigned int numberOfSmallToGo)
{
    logger(DEBUG, "need % small and % large particles", numberOfLargeToGo, numberOfSmallToGo);
    
    while ((numberOfLargeToGo > 0 || numberOfSmallToGo > 0))
    {
        if (numberOfSmallToGo > 0)
        {
            insertSmallParticle(numberOfSmallToGo);
        }
        else
        {
            insertLargeParticle(numberOfLargeToGo);
        }
    }
}

void PeriodicChute::insertLargeParticle(unsigned int& numberOfLargeToGo)
{
    BaseParticle p0;
    p0.setSpecies(speciesHandler.getObject(1));
    const Mdouble radius = random.getRandomNumber((2.0 / 3.0 * 0.9), (2.0 / 3.0 * 1.1));
    p0.setRadius(radius);
    insertOneParticle(p0);
    numberOfLargeToGo--;
    logger(DEBUG, "Inserted large particle");
}

void PeriodicChute::insertSmallParticle(unsigned int& numberOfSmallToGo)
{
    BaseParticle p0;
    p0.setSpecies(speciesHandler.getObject(2));
    const Mdouble radius = random.getRandomNumber((1.0 / 3.0 * 0.9), (1.0 / 3.0 * 1.1));
    p0.setRadius(radius);
    insertOneParticle(p0);
    numberOfSmallToGo--;
    logger(DEBUG, "Inserted small particle");
}

//assume: species and radius of p0 are set
void PeriodicChute::insertOneParticle(BaseParticle& p0)
{
    Vec3D pos;
    p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
    unsigned int toBeFailed = 2000;
    do
    {
        //if we failed to insert the particle 1000 times, increase domain in z direction
        if (toBeFailed == 0)
        {
            setZMax(getZMax() + 2);
            toBeFailed = 1000;
            logger(DEBUG, "resetting failed and increasing zmax, new zmax: %", getZMax());
        }
        pos.X = random.getRandomNumber(getXMin() + p0.getRadius(), getXMax() - p0.getRadius());
        pos.Y = random.getRandomNumber(getYMin() + p0.getRadius(), getYMax() - p0.getRadius());
        pos.Z = random.getRandomNumber(getZMin() + p0.getRadius(), getZMax());
        p0.setPosition(pos);
        toBeFailed--;
    } while (!checkParticleForInteraction(p0));
    particleHandler.copyAndAddObject(p0);
}

void PeriodicChute::setBoundaries()
{
    
    PeriodicBoundary b0;
    b0.set(Vec3D(0.0, 1.0, 0.0), getYMin(), getYMax());
    boundaryHandler.copyAndAddObject(b0);
    b0.set(Vec3D(1.0, 0.0, 0.0), getXMin(), getXMax());
    boundaryHandler.copyAndAddObject(b0);
}

void PeriodicChute::setChuteProperties()
{
    setGravity({std::sin(26.0 * constants::pi / 180), 0, -std::cos(26.0 * constants::pi / 180)});
    setXMax(30);
    setYMax(10);
    setZMax(4);
}

void PeriodicChute::setTimeAndSaveCount()
{
    setTimeStep(1e-4);
    if (isMonoDisperse)
    {
        setTimeMax(1000);
    }
    else
    {
        setTimeMax(2000);
    }
    setSaveCount(1e4);
}

void PeriodicChute::setSpeciesProperties()
{
    const Mdouble density = 6.0 / constants::pi;
    
    auto sReference = LinearViscoelasticSlidingFrictionSpecies();
    //for the reference particles (d=1, m=1, see silbert):
    sReference.setDensity(density);
    sReference.setCollisionTimeAndRestitutionCoefficient(0.005, 0.88, 1);
    sReference.setSlidingDissipation(2.0 / 7 * sReference.getDissipation());
    sReference.setSlidingStiffness(2.0 / 7 * sReference.getStiffness());
    sReference.setSlidingFrictionCoefficient(0.5);
    //add species to the handler twice: once for bottom, one for bottom, purely for visualisation purposes
    speciesHandler.copyAndAddObject(sReference);
    if (isMonoDisperse)
    {
        speciesHandler.copyAndAddObject(sReference);
    }
    else
    {
        const Mdouble massLarge = 4.0 / 3 * constants::pi * pow(2.0 / 3.0, 3.0) * density;
        auto sLarge = sReference;
        sLarge.setCollisionTimeAndRestitutionCoefficient(0.005, 0.88, massLarge);
        sLarge.setSlidingDissipation(2.0 / 7 * sLarge.getDissipation());
        sLarge.setSlidingStiffness(2.0 / 7 * sLarge.getStiffness());
        sLarge.setSlidingFrictionCoefficient(0.5);
        
        const Mdouble massSmall = 4.0 / 3 * constants::pi * pow(1.0 / 3.0, 3.0) * density;
        auto sSmall = sReference;
        sSmall.setCollisionTimeAndRestitutionCoefficient(0.005, 0.88, massSmall);
        sSmall.setSlidingDissipation(2.0 / 7 * sSmall.getDissipation());
        sSmall.setSlidingStiffness(2.0 / 7 * sSmall.getStiffness());
        sSmall.setSlidingFrictionCoefficient(0.5);
        
        speciesHandler.copyAndAddObject(sLarge);
        speciesHandler.copyAndAddObject(sSmall);
    }
}
