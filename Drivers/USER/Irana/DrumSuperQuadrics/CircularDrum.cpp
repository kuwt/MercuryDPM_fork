
/// \todo Make the scale factor of moving the surface of the drum particles closer together a setable paramter
/// \todo Make it four stages with a counter
/// Stage 1 : Create the drum
/// Stage 2 : Get the particles in : At the moment due to how to this checked only low particles numbers can be placed in.
/// Stage 3 : Settle the particles
/// Repeat stage 2 until all required particles are in
/// Stage 4 : Settle to a very low KE
/// Stage 5 : Rotate the drum
/// \todo write a restarter for this code
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include <CG/Functions/Lucy.h>
#include <CG/TimeAveragedCG.h>
#include "Chute.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"
#include "MercuryTime.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"

enum class Stage
{
    SETTLE_PARTICLES, RELAX_AFTER_INSERTION, ROTATE_DRUM
};

/// \brief Class for simulation rotating drums whose outer walls are made of particles
class RotatingDrum : public Chute
{
public:
    /// \brief Default constructor call. Sets the diameter of the small particles to be one.
    RotatingDrum() : radius_s(0.5)
    {
        volumeRatioEllipsoidOverSphere = 1;
        aspectRatioEllipsoid = 2;
        setMin(0, 0, 0);
        scaleFactor = 1;
        checkTime = 1;
        numSphereToBeInserted = 1000;
        numEllipsoidToBeInserted = numSphereToBeInserted / volumeRatioEllipsoidOverSphere;
        setName("drumEllipsoidsV" + std::to_string((int) (10 * volumeRatioEllipsoidOverSphere)) + "AR" +
                std::to_string((int) aspectRatioEllipsoid));
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////
    /// /brief setupInitial conditions. Basics does step 1 only; creating the walls of the drum
    ////////////////////////////////////////////////////////////////////////////////////////
    void setupInitialConditions()
    {
        
        setTimeMax(15 * 2 * constants::pi / rotationSpeed); //simulate 15 rotations of the drum.
        logger(INFO, "time max: %", getTimeMax());
        setTimeStep(1. / (200.0 * 5.0));
        
        setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(301, getTimeMax(), getTimeStep()));
        
        addSpecies();
        
        constructFixedParticleDrum();
        
        // Remove the two existing boundaries insertion and periodic and put the periodic back if periodic else put
        boundaryHandler.clear();
        wallHandler.clear();
        
        if (getIsPeriodic())
        {
            PeriodicBoundary b0;
            b0.set(Vec3D(0.0, 1.0, 0.0), getYMin(), getYMax());
            boundaryHandler.copyAndAddObject(b0);
        }
        else
        {
            //Now add two extra solid walls at the end
            InfiniteWall w0;
            w0.setSpecies(speciesDrumWall);
            w0.set(Vec3D(0.0, -1.0, 0.0), getMin() + Vec3D(0, 2 * getFixedParticleRadius(), 0));
            wallHandler.copyAndAddObject(w0);
            w0.set(Vec3D(0.0, 1.0, 0.0), getMax() - Vec3D(0, 2 * getFixedParticleRadius(), 0));
            wallHandler.copyAndAddObject(w0);
        }
        insertParticles();
        
        auto cg0 = cgHandler.copyAndAddObject(TimeAveragedCG<CGCoordinates::XZ, CGFunctions::Lucy>());
        cg0->setWidth(1);
        cg0->statFile.setSaveCount(10. / getTimeStep());
        cg0->statFile.setName(getName() + "." + std::to_string(getRunNumber()) + ".LucyXZ.stat");
        cg0->setNZ(40);
        cg0->setNX(40);
        cg0->setX(getXMin(), getXMax());
        cg0->setZ(getZMin(), getZMax());
        cg0->setTimeMin(getTimeMax() - 100);
        auto cg1 = cgHandler.copyAndAddObject(cg0);
        cg1->setSelectedParticle([](const BaseInteractable* p)
                                 { return p->getSpecies()->getIndex() == 1; });
        cg1->statFile.setName(getName() + "." + std::to_string(getRunNumber()) + ".species1.LucyXZ.stat");
        auto cg2 = cgHandler.copyAndAddObject(cg0);
        cg2->setSelectedParticle([](const BaseInteractable* p)
                                 { return p->getSpecies()->getIndex() == 2; });
        cg2->statFile.setName(getName() + "." + std::to_string(getRunNumber()) + ".species2.LucyXZ.stat");
    }
    
    void constructFixedParticleDrum()
    {
        Chute::setupInitialConditions();
        convertBaseParticlesToSuperquadric();
        
        setMin(-drumRadius, getYMin(), -drumRadius);
        setMax(drumRadius, getYMax(), drumRadius);
        
        //Now map the bed particles to the drum
        logger(INFO, "\n \n \n STEP 1 : Creating the Drum\n---------------------------\n \n \n");
        
        for (unsigned int i = 0; i < particleHandler.getNumberOfObjects(); i++)
        {
            BaseParticle* p = particleHandler.getObject(i);
            Vec3D position = p->getPosition();
            Mdouble theta = (position.X) / (drumRadius * scaleFactor);
            Mdouble y = position.Y;
            //simple, cylindrical drum
            Mdouble r = drumRadius * (1 + position.Z / (drumRadius * constants::pi * 2.0));
            
            position.X = r * cos(theta);
            position.Y = y;
            position.Z = r * sin(theta);
            
            p->setPosition(position);
        }
        
        logger(INFO, "rotation speed %", rotationSpeed);
    }
    
    void convertBaseParticlesToSuperquadric()
    {
        const unsigned numberBaseParticles = particleHandler.getSize();
        for (unsigned i = 0; i < numberBaseParticles; ++i)
        {
            BaseParticle* particle = particleHandler.getObject(i);
            SuperQuadricParticle particleSuperQuadric = SuperQuadricParticle(*particle);
            particleHandler.copyAndAddObject(particleSuperQuadric);
        }
        for (int i = numberBaseParticles - 1; i >= 0; --i)
        {
            particleHandler.removeObject(i);
        }
    }
    
    void addSpecies()
    {
        Mdouble mu = 0.35;
        Mdouble muRolling = 0;
        Mdouble muTorsion = 0;
        // Set the non-dim density
        double rho1 = 6.0 / constants::pi;
        
        // Set the collision time and coefficient of resitution.
        double tc = 1.0 / 20.0;
        
        //Computer the mass of the partices
        double massSphere = 4.0 / 3.0 * constants::pi * std::pow(radius_s, 3.0) * rho1;
        double massEllipsoid = massSphere * volumeRatioEllipsoidOverSphere;
        double massWall = 4.0 / 3.0 * constants::pi * std::pow(0.25, 3.0) * rho1;
        logger(INFO, "mass small %, radius_s %", massSphere, radius_s);
        
        speciesDrumWall = new LinearViscoelasticFrictionSpecies;
        speciesParticles1 = new LinearViscoelasticFrictionSpecies;
        speciesParticles2 = new LinearViscoelasticFrictionSpecies;
        speciesDrumWall->setIsSuperquadricSpecies(true);
        speciesParticles1->setIsSuperquadricSpecies(true);
        speciesParticles2->setIsSuperquadricSpecies(true);
        
        speciesDrumWall->setDensity(rho1);
        speciesDrumWall->setCollisionTimeAndRestitutionCoefficient(tc, wallCOR, massWall);
        speciesDrumWall->setSlidingDissipation(
                speciesDrumWall->getDissipation()); //  Set the tangential dissipation equal to the normal disipation for small-small collsions
        speciesDrumWall->setSlidingStiffness(speciesDrumWall->getStiffness() * 2.0 / 7.0);
        speciesDrumWall->setSlidingFrictionCoefficient(mu);
        //////
        //
        speciesParticles1->setDensity(rho1);
        speciesParticles1->setCollisionTimeAndRestitutionCoefficient(tc, smallCOR, massSphere);
        speciesParticles1->setSlidingDissipation(
                speciesParticles1->getDissipation()); //  Set the tangential dissipationequal to the normal disipation for large-large collision
        speciesParticles1->setSlidingStiffness(speciesParticles1->getStiffness() * 2.0 / 7.0);
        speciesParticles1->setSlidingFrictionCoefficient(mu);
        
        speciesParticles1->setRollingStiffness(speciesParticles1->getStiffness() * 2.0 / 5.0);
        speciesParticles1->setRollingFrictionCoefficient(muRolling);
        speciesParticles1->setRollingDissipation(speciesParticles1->getDissipation());
        
        speciesParticles1->setTorsionStiffness(speciesParticles1->getStiffness() * 2.0 / 5.0);
        speciesParticles1->setTorsionFrictionCoefficient(muTorsion);
        speciesParticles1->setTorsionDissipation(speciesParticles1->getDissipation());
        
        //////////
        speciesParticles2->setDensity(rho1);
        speciesParticles2->setCollisionTimeAndRestitutionCoefficient(tc, largeCOR, massEllipsoid);
        speciesParticles2->setSlidingDissipation(speciesParticles2->getDissipation());
        speciesParticles2->setSlidingStiffness(speciesParticles2->getStiffness() * 2.0 / 7.0);
        speciesParticles2->setSlidingFrictionCoefficient(mu);
        
        speciesParticles2->setRollingStiffness(speciesParticles1->getStiffness() * 2.0 / 5.0);
        speciesParticles2->setRollingFrictionCoefficient(muRolling);
        speciesParticles2->setRollingDissipation(speciesParticles1->getDissipation());
        
        speciesParticles2->setTorsionStiffness(speciesParticles1->getStiffness() * 2.0 / 5.0);
        speciesParticles2->setTorsionFrictionCoefficient(muTorsion);
        speciesParticles2->setTorsionDissipation(speciesParticles1->getDissipation());
        
        speciesHandler.addObject(speciesDrumWall);
        speciesHandler.addObject(speciesParticles1);
        speciesHandler.addObject(speciesParticles2);
        speciesMixedDrumAnd1 = speciesHandler.getMixedObject(speciesDrumWall, speciesParticles1);
        speciesMixedDrumAnd2 = speciesHandler.getMixedObject(speciesDrumWall, speciesParticles2);
        speciesMixed1And2 = speciesHandler.getMixedObject(speciesParticles1, speciesParticles2);
        speciesMixedDrumAnd1->setIsSuperquadricSpecies(true);
        speciesMixedDrumAnd2->setIsSuperquadricSpecies(true);
        speciesMixed1And2->setIsSuperquadricSpecies(true);
        
        
        //  Set the contact time (tc), resitution coefficeient (r) and density (rho) for small for all particles
        
        
        speciesMixedDrumAnd1->setCollisionTimeAndRestitutionCoefficient(tc, (0.5 * (wallCOR + smallCOR)), massSphere,
                                                                        massSphere);
        speciesMixedDrumAnd1->setSlidingDissipation(
                speciesMixedDrumAnd1->getDissipation()); //  Set the tangential dissipation equal to the normal disipation for mixed collision
        speciesMixedDrumAnd1->setSlidingFrictionCoefficient(2 * mu);
        speciesMixedDrumAnd1->setSlidingStiffness(speciesMixedDrumAnd1->getStiffness() * 2.0 / 7.0);
        speciesMixedDrumAnd2->setCollisionTimeAndRestitutionCoefficient(tc, (0.5 * (wallCOR + largeCOR)), massSphere,
                                                                        massEllipsoid);
        
        speciesMixedDrumAnd1->setRollingStiffness(speciesParticles1->getStiffness() * 2.0 / 5.0);
        speciesMixedDrumAnd1->setRollingFrictionCoefficient(muRolling);
        speciesMixedDrumAnd1->setRollingDissipation(speciesParticles1->getDissipation());
        
        speciesMixedDrumAnd1->setTorsionStiffness(speciesParticles1->getStiffness() * 2.0 / 5.0);
        speciesMixedDrumAnd1->setTorsionFrictionCoefficient(muTorsion);
        speciesMixedDrumAnd1->setTorsionDissipation(speciesParticles1->getDissipation());
        
        speciesMixedDrumAnd2->setSlidingDissipation(
                speciesMixedDrumAnd2->getDissipation()); //  Set the tangential dissipation equal to the normal disipation
        speciesMixedDrumAnd2->setSlidingFrictionCoefficient(2 * mu);
        speciesMixedDrumAnd2->setSlidingStiffness(speciesMixedDrumAnd2->getStiffness() * 2.0 / 7.0);
        
        speciesMixedDrumAnd2->setRollingStiffness(speciesParticles1->getStiffness() * 2.0 / 5.0);
        speciesMixedDrumAnd2->setRollingFrictionCoefficient(muRolling);
        speciesMixedDrumAnd2->setRollingDissipation(speciesParticles1->getDissipation());
        
        speciesMixedDrumAnd2->setTorsionStiffness(speciesParticles1->getStiffness() * 2.0 / 5.0);
        speciesMixedDrumAnd2->setTorsionFrictionCoefficient(muTorsion);
        speciesMixedDrumAnd2->setTorsionDissipation(speciesParticles1->getDissipation());
        
        speciesMixed1And2->setCollisionTimeAndRestitutionCoefficient(tc, (0.5 * (largeCOR + smallCOR)), massSphere,
                                                                     massEllipsoid);
        speciesMixed1And2->setSlidingDissipation(
                speciesMixed1And2->getDissipation()); //  Set the tangential dissipation equal to the normal disipation
        speciesMixed1And2->setSlidingFrictionCoefficient(mu);
        speciesMixed1And2->setSlidingStiffness(speciesMixed1And2->getStiffness() * 2.0 / 7.0);
        
        speciesMixed1And2->setRollingStiffness(speciesParticles1->getStiffness() * 2.0 / 5.0);
        speciesMixed1And2->setRollingFrictionCoefficient(muRolling);
        speciesMixed1And2->setRollingDissipation(speciesParticles1->getDissipation());
        
        speciesMixed1And2->setTorsionStiffness(speciesParticles1->getStiffness() * 2.0 / 5.0);
        speciesMixed1And2->setTorsionFrictionCoefficient(muTorsion);
        speciesMixed1And2->setTorsionDissipation(speciesParticles1->getDissipation());
        
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Create and inserts in the drum new partciles random at non-overlapping locations in the drum.
    /// It also checks the locations is not already filled by a particle
    ////////////////////////////////////////////////////////////////////////////////////////////
    void createParticles()
    {
        // CREATE THE PARTICLES
        while ((numSphereToBeInserted > 0) || (numEllipsoidToBeInserted > 0))
        {
            bool isLargeParticle = false;
            SuperQuadricParticle P0;
            P0.setExponents(1, 1);
            //random to see if want to generate a large or small particles, helps makes the initial conditions homogenious
            if (random.getRandomNumber(1.0, numEllipsoidToBeInserted + numSphereToBeInserted) >
                (numEllipsoidToBeInserted))
            {
                P0.setAxes(0.5, 0.5, 0.5);
                P0.setAxes(1.8, std::sqrt(0.2), std::sqrt(0.2));
                P0.setSpecies(speciesParticles1);
            }
            else
            {
                Mdouble b = std::cbrt(1.0 / 8 / aspectRatioEllipsoid) * std::cbrt(volumeRatioEllipsoidOverSphere);
                Mdouble a = aspectRatioEllipsoid * b;
                P0.setAxes(a, b, b);
                P0.setAxes(1, 0.6, 0.6);
                P0.setSpecies(speciesParticles2);
                isLargeParticle = true;
            }
            //randomise particle position, zero initial velocity
            unsigned failCounter = 0;
            do
            {
                Mdouble r = random.getRandomNumber(-0.9, 0.9) * drumRadius;
                Mdouble theta = random.getRandomNumber(0, constants::pi * 2.0);
                Mdouble y = random.getRandomNumber(getYMin(), getYMax());
                Vec3D position;
                position.X = r * std::cos(theta);
                position.Y = y;
                position.Z = r * std::sin(theta);
                
                P0.setPosition(position);
                P0.setVelocity(Vec3D(0.0, 0.0, 0.0));
                
                failCounter++;
                
                if (failCounter == 10000)
                {
                    return;
                }
                
            } while (!checkParticleForInteraction(P0));
            
            particleHandler.copyAndAddObject(P0);
            if (isLargeParticle)
            {
                numEllipsoidToBeInserted--;
            }
            else
            {
                numSphereToBeInserted--;
            }
        }
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// /brief actionsBeforeTimeStep. This does stage 3: setlle particles, stage 4: relax particles and stage 5: start the drum rotating.
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    void actionsBeforeTimeStep()
    {
        if (step == Stage::ROTATE_DRUM || getTime() < checkTime)
        {
            return;
        }
        checkTime = getTime() + 1;
        std::cout << "Current KE " << getKineticEnergy() << std::endl;
        if (step == Stage::RELAX_AFTER_INSERTION)
        {
            if (getKineticEnergy() < (particleHandler.getNumberOfObjects() / 1000.0))
            {
                initialiseDrumRotation();
                step = Stage::ROTATE_DRUM;
                logger(INFO, "\n \n \nSTEP 5: Starting the drum rotation\n---------------------------\n \n \n");
                logger(INFO, "Time is %", getTime());
            }
        }
        
        if (step == Stage::SETTLE_PARTICLES && getKineticEnergy() < (particleHandler.getNumberOfObjects() / 100.0))
        {
            insertParticles();
        }
    }
    
    void actionsOnRestart() override
    {
        //offset is the time that the drum starts rotating
        Mdouble offset = 13;
        setName("restarted" + getName());
        drumRadius = 15;
        setFroude(0.01);
        logger(INFO, "rotation speed %", rotationSpeed);
        //cg from 10 rotations for 1/2 rotation
        setTimeMax( offset + 20 * constants::pi  / rotationSpeed + constants::pi / rotationSpeed);
        logger(INFO, "time max: %", getTimeMax());
        initialiseDrumRotation();
        step = Stage::ROTATE_DRUM;
        auto cg0 = cgHandler.copyAndAddObject(TimeAveragedCG<CGCoordinates::XZ, CGFunctions::Lucy>());
        cg0->setWidth(1);
        cg0->statFile.setSaveCount(10. / getTimeStep());
        cg0->statFile.setName(getName() + "." + ".LucyXZ.stat");
        cg0->setNZ(60);
        cg0->setNX(60);
        cg0->setX(getXMin(), getXMax());
        cg0->setZ(getZMin(), getZMax());
        cg0->setTimeMin(offset + 20 * constants::pi  / rotationSpeed);
        auto cg1 = cgHandler.copyAndAddObject(cg0);
        cg1->setSelectedParticle([](const BaseInteractable* p)
                                 { return p->getSpecies()->getIndex() == 1; });
        cg1->statFile.setName(getName() + ".species1.LucyXZ.stat");
        auto cg2 = cgHandler.copyAndAddObject(cg0);
        cg2->setSelectedParticle([](const BaseInteractable* p)
                                 { return p->getSpecies()->getIndex() == 2; });
        cg2->statFile.setName(getName() +  ".species2.LucyXZ.stat");
        
        auto cg3 = cgHandler.copyAndAddObject(
                TimeAveragedCG<CGCoordinates::XZ, CGFunctions::Lucy, CGFields::OrientationField>());
        cg3->setWidth(1);
        cg3->statFile.setName(getName() + ".orientation.stat");
        cg3->setNZ(60);
        cg3->setNX(60);
        cg3->setX(getXMin(), getXMax());
        cg3->setZ(getZMin(), getZMax());
        cg3->setSelectedParticle([](const BaseInteractable* p)
                                 { return p->getSpecies()->getIndex() == 2; });
        cg3->setTimeMin(offset + 20 * constants::pi  / rotationSpeed);
    }
    
    void initialiseDrumRotation()
    {
        prescribeFixedParticlePositions();
    }
    
    void insertParticles()
    {
        logger(INFO, "\n \n \nSTEP 2: Inserting particles\n---------------------------\n \n \n");
        logger(INFO, "Number of large particles to be inserted: %", numEllipsoidToBeInserted);
        logger(INFO, "Number of small particles to be inserted: %", numSphereToBeInserted);
        createParticles();
        logger(INFO, "Finished creating particles");
        logger(INFO, "Number of large particles still to be inserted: %", numEllipsoidToBeInserted);
        logger(INFO, "Number of small particles still to be inserted: %", numSphereToBeInserted);
        if ((numSphereToBeInserted == 0) && (numEllipsoidToBeInserted == 0))
        {
            // If you are here are partices are inserted and you are moving to step 4; relax to very low KE
            step = Stage::RELAX_AFTER_INSERTION;
            logger(INFO, "\n \n \nSTEP 4: Relaxing particles\n---------------------------\n \n \n");
        }
        else
        {
            // If you are here you are the dum is not full and you are settling the inserted particles
            step = Stage::SETTLE_PARTICLES;
            logger(INFO, "\n \n \nSTEP 3: Settling the inserted particles\n---------------------------\n \n \n");
        }
    }
    
    void prescribeFixedParticlePositions() const
    {
        for (BaseWall* w : wallHandler)
        {
            w->setPrescribedAngularVelocity([this](Mdouble time)
                                            { return Vec3D(0, rotationSpeed, 0); });
        }
        for (BaseParticle* p : particleHandler)
        {
            if (p->getInvMass() == 0)
            {
                p->setPrescribedAngularVelocity([this](Mdouble time)
                                                { return Vec3D(0, rotationSpeed, 0); });
                p->setPrescribedPosition([this, p](Mdouble time)
                                         {
                                             const Vec3D position = p->getPosition();
                                             const Mdouble theta =
                                                     std::atan2(position.Z, position.X) + getTimeStep() * rotationSpeed;
                    
                                             const Mdouble r = std::sqrt(
                                                     position.X * position.X + position.Z * position.Z);
                    
                                             return (Vec3D(
                                                     (r * std::cos(theta)),
                                                     position.Y,
                                                     (r * std::sin(theta))));
                                         });
                
                p->setPrescribedVelocity([this, p](Mdouble time)
                                         {
                                             const Vec3D position = p->getPosition();
                                             Mdouble theta = std::atan2(position.Z, position.X);
                                             const Mdouble r = std::sqrt(
                                                     position.X * position.X + position.Z * position.Z);
                    
                                             return (rotationSpeed * r * Vec3D(
                                                     -std::sin(theta),
                                                     0,
                                                     std::cos(theta)));
                                         });
            }
        }
    }
    
    void setRPM(double new_speed)
    {
        rotationSpeed = new_speed;
    }
    
    void setFroude(Mdouble froude)
    {
        logger.assert_always(drumRadius > 1e-10, "Set drum radius before froude number");
        rotationSpeed = std::sqrt(froude * getGravity().getLength() / drumRadius);
    }
    
    void set_particle_numbers(int new_num_small, int new_num_large)
    {
        logger.assert_always(new_num_small >= 0, "Please give a non-negative number of small particles");
        logger.assert_always(getTime() < getTimeStep(), "Please insert particle-number before starting simulation");
        numSphereToBeInserted = new_num_small;
        
        logger.assert_always(new_num_large >= 0, "Please give a non-negative number of large particles");
        numEllipsoidToBeInserted = new_num_large;
    }
    
    
    void setDrumRadius(double radius)
    {
        drumRadius = radius;
        
        setXMax(drumRadius * constants::pi * 2.0 * scaleFactor);
        setZMax(drumRadius * constants::pi * 2.0);
    }
    
    /// Set the particle size of the wall particles in dimeter and overlap between
    void setWallParameters(double particlesSize, double new_overlap)
    {
        setFixedParticleRadius(particlesSize);
        scaleFactor = new_overlap;
        
        setXMax(drumRadius * constants::pi * 2.0 * scaleFactor);
        setZMax(drumRadius * constants::pi * 2.0);
    }
    
    void setCoefficientOfRestitutionLarge(double newCOR)
    {
        largeCOR = newCOR;
    }
    
    void setCoefficientOfRestitutionSmall(double newCOR)
    {
        smallCOR = newCOR;
    }
    
    void setCoefficientOfRestitutionWall(double newCOR)
    {
        wallCOR = newCOR;
    }

private:
    double rotationSpeed;
    double volumeRatioEllipsoidOverSphere;
    double aspectRatioEllipsoid;
    double smallCOR;
    double largeCOR;
    double wallCOR;
    const double radius_s;
    int numSphereToBeInserted;
    int numEllipsoidToBeInserted;
    Stage step;
    
    // This is the factor the particles in the drum are overlapped by 1.0 means no overlap; 2.0 means 50% overlap etc.
    double scaleFactor;
    
    double drumRadius;
    Mdouble checkTime;
    LinearViscoelasticFrictionSpecies* speciesDrumWall;
    LinearViscoelasticFrictionSpecies* speciesParticles1;
    LinearViscoelasticFrictionSpecies* speciesParticles2;
    LinearViscoelasticFrictionMixedSpecies* speciesMixedDrumAnd1;
    LinearViscoelasticFrictionMixedSpecies* speciesMixedDrumAnd2;
    LinearViscoelasticFrictionMixedSpecies* speciesMixed1And2;
};

int main(int argc, char* argv[])
{
    
    // Problem parameters
    RotatingDrum problem;
    problem.autoNumber();
    
    ///\todo set axes of ellipsoids here somehow
    //problem.setCoefficientOfRestitution(0.9);
    problem.setCoefficientOfRestitutionLarge(0.1);
    problem.setCoefficientOfRestitutionSmall(0.1);
    problem.setCoefficientOfRestitutionWall(0.1);
    
    problem.setChuteAngleAndMagnitudeOfGravity(0.0, 1.0);
    
    // Chute properties : Simply remove the first line to add side walls.
    //problem.makeChutePeriodic();
    //problem.setYMax(40.5);
    problem.setYMax(6);
    problem.makeChutePeriodic();
    //problem.setDrumRadius(35.5);
    problem.setDrumRadius(27);
    problem.setFroude(0.01);
    
    
    //set fixed-particle radius and scale factor (number of times the bottom is wrapped around the drum)
    problem.setWallParameters(0.25, 1);
    
    //Swap the next two lines to swap between the different type of rought bottoms. Note MONOLAYER_ORDERED does not give enough friction for the circular case.
    //problem.setRoughBottomType(MULTILAYER);
    //problem.setRoughBottomType(MONOLAYER_DISORDERED);
    problem.setRoughBottomType(MONOLAYER_ORDERED);
    
    problem.restartFile.setFileType(FileType::MULTIPLE_FILES);
    problem.setXBallsAdditionalArguments("-cmode 8");
    problem.setSuperquadricParticlesWriteVTK(true);
    problem.readArguments(argc, argv);
    
    Time time1;
    time1.tic();
    
    problem.solve();
    std::cout << "Time to solve: " << time1.toc() << std::endl;
    
}
