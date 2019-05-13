//
// Created by irana on 10/3/17.
//

#include <iomanip>
#include <random>
#include "Chute.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "Boundaries/DeletionBoundary.h"
#include "Boundaries/SubcriticalMaserBoundary.h"
#include <random>

class ChuteWithMaserAndSideWalls : public Chute
{
public:
    /*! Constructor: construct the ChuteWithMaserAndSideWalls simulation, and write one file for all types except restart files;
     * write a new restart file every time.
     */
    ChuteWithMaserAndSideWalls() :
            Chute(),
            maserHasOpened(false)
    {
        setFileType(FileType::ONE_FILE);
        restartFile.setFileType(FileType::MULTIPLE_FILES);
        setChuteProperties();
        setMin(0, 0, 0);
        setMax(getChuteLength(), getChuteWidth(), 4);
        setSpeciesProperties();
        setTimeStep(1e-4);
    }
    
    void setupInitialConditions() override
    {
        logger.assert(!getIsPeriodic(), "Chute should not be periodic in y-direction");
        setupSideWalls();
        PeriodicBoundary b0;
        b0.set(Vec3D(1.0,0.0,0.0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(b0);
        
        createBottom();
        insertParticles();
    
        setSaveCount(.1 / getTimeStep());
        restartFile.setSaveCount(10./getTimeStep());
    }


/*!
 * A subfunction of setupInitialConditions; determines the initial content of the particleHandler.
 * It inserts the particles randomly into the domain in order to achieve a homogeneous distribution.
 * Note, large particles are assigned to species 1, small particles are assigned species 0,
 * fixed particles are species 0 (by default).
 */
    void insertParticles()
    {
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(1));
        while (particleHandler.getSize() - particleHandler.getNumberOfFixedObjects()
               < getChuteLength() * getChuteWidth() * getInflowHeight())
        {
            p0.setRadius(0.5 * distribution(generator));
            Vec3D pos;
            p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
            unsigned int toBeFailed = 1000;
            do
            {
                //if we failed to insert the particle 1000 times, increase domain in z direction
                if (toBeFailed == 0)
                {
                    setZMax(getZMax() + 1);
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
    }
    
    /*!
     * A function that overrides DPMBase::actionsBeforeTimeStep; opens the maser at a certain time by extending the base
     * and replacing the periodic boundary with a maser boundary. This also sets a flag that the maser is activated
     */
    void actionsBeforeTimeStep() override
    {
        if (!maserHasOpened && getTime() > 5)
        {
            populateBottomParticles();
            copyParticlesForMaser();
            openMaser();
        }
        if (maserHasOpened)
        {
            //dynamically extend the bottom every 1t if the particles moved too close to the edge.
            static Mdouble nextTimeToCheck = getTime() + getTimeStep();
            if (getTime() > nextTimeToCheck && bottomLength < 1e4)
            {
                Mdouble furthestX = 0;
                for (const BaseParticle* const p : particleHandler)
                {
                    if (!p->isFixed() && p->getPosition().X > furthestX)
                    {
                        furthestX = p->getPosition().X;
                    }
                }
                while (furthestX + 200 > bottomLength)
                {
                    //copy BaseParticles up to x = nextLength
                    for (BaseParticle* p : bottomParticles)
                    {
                        p->setPosition(p->getPosition() + Vec3D(20, 0, 0));
                        particleHandler.copyAndAddObject(p);
                    }
                    bottomLength += 20;
                }
                
                nextTimeToCheck += 5;
                logger(INFO, "New bottom length: %, t = %", bottomLength, getTime());
            }
        }
    }
    
    //this function replaces setupInitialConditions in case of a restart
    //if you have problem-specific variables, set them here
    void actionsOnRestart() override
    {
        setName(getName());
        if (getTime() < 2000.001)
        {
            maserHasOpened = false;
        }
        else
        {
            maserHasOpened = true;
            bottomLength = particleHandler.getHighestPositionX();
            populateBottomParticles();
        }
        setName("maser" + getName());
        restartFile.setFileType(FileType::MULTIPLE_FILES);
        restartFile.setSaveCount(10./getTimeStep());
        dataFile.setSaveCount(10./getTimeStep());
    }

private:
    
    /*!
 * A subfunction of the constructor; sets species properties such that all particles have the same density (6/pi),
 * and all collisions have the same contact time tc=sqrt(<d>/g)/200
 */
    void setSpeciesProperties()
    {
        const Mdouble density = 6.0 / constants::pi;
    
        auto sReference = LinearViscoelasticSlidingFrictionSpecies();
        //for the reference particles (d=1, m=1, see silbert):
        sReference.setDensity(density);
        sReference.setDissipation(25); //gamma^n
        sReference.setSlidingDissipation(2.0 / 7 * sReference.getDissipation()); //  gamma^t
        sReference.setStiffness(2e5); // k^n
        sReference.setSlidingStiffness(2.0 / 7 * sReference.getStiffness()); // k^t
        sReference.setSlidingFrictionCoefficient(0.5); //mu
    
        //add species to handler twice: once for bottom, once for flow (nicer for visualisation)
        speciesHandler.copyAndAddObject(sReference);
        speciesHandler.copyAndAddObject(sReference);
    }
    
    /*!
 * A subfunction of the constructor; sets gravity vector, domain size, and the rough bottom,
 * which is made of particles to simulate a rough surface.
 */
    void setChuteProperties()
    {
        setChuteAngleAndMagnitudeOfGravity(24, 1);
        setChuteLength(20);
        setChuteWidth(10);
        setInflowHeight(5);
        setZMax(4);
        setRoughBottomType(RoughBottomType::MULTILAYER);
        setFixedParticleRadius(0.5);
    }
    
    void populateBottomParticles()
    {
        for (const BaseParticle* const p : particleHandler)
        {
            if (p->isFixed() && p->getPosition().X < 20)
            {
                BaseParticle* pCopy = p->copy();
                logger(DEBUG, "species of bottom particle: %", p->getSpecies()->getIndex());
                bottomParticles.push_back(pCopy);
            }
        }
    }
    
    ///copy the initial set of particles, here the bottom. You can add flow particles here as well (optional).
    void copyParticlesForMaser()
    {
        //copy BaseParticles up to x = 200
        for (BaseParticle* p : bottomParticles)
        {
            while (p->getPosition().X < 200)
            {
                p->setPosition(p->getPosition() + Vec3D(20, 0, 0));
                particleHandler.copyAndAddObject(p);
            }
        }
        bottomLength = 200;
    }
    
    /*!
     * Subfunction of actionsBeforeTimeStep; replaces the periodic boundary in x-direction with a maser boundary
     */
    void openMaser()
    {
        logger.assert_always((dynamic_cast<PeriodicBoundary*> (boundaryHandler.getLastObject())) != nullptr,
                             "Cannot cast % to maser boundary", boundaryHandler.getLastObject()->getName());
        SubcriticalMaserBoundary maserBoundaryActual(*(dynamic_cast<PeriodicBoundary*> (boundaryHandler.getLastObject())));
        boundaryHandler.removeLastObject();
        boundaryHandler.copyAndAddObject(maserBoundaryActual);
        
        //make a deletion boundary for the particles falling off the end, here at x = 10'000
        DeletionBoundary deletionBoundary;
        deletionBoundary.set(Vec3D(1, 0, 0), 1e4);
        
        maserHasOpened = true;
        logger(INFO, "opened the maser boundary");
        
        setXMax(1e4);
    }
    
    Mdouble bottomLength = 20;
    bool maserHasOpened;
    std::vector<BaseParticle*> bottomParticles;
    
    std::default_random_engine generator;
    std::normal_distribution<Mdouble> distribution = std::normal_distribution<Mdouble>(1.0, 0.025);;
};

int main(int argc, char* argv[])
{
    ChuteWithMaserAndSideWalls problem;
    problem.setName("ChuteWithMaserAndSideWalls");
    problem.setTimeMax(100);
    problem.solve(argc,argv);
    return 0;
}
