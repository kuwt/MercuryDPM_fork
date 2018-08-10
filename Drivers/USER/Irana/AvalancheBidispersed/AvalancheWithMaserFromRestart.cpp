//
// Created by irana on 8/4/17.
//

#include <iomanip>
#include <sstream>
#include <Boundaries/DeletionBoundary.h>
#include <Boundaries/SubcriticalMaserBoundary.h>
#include <CG/Functions/Lucy.h>
#include <CG/TimeAveragedCG.h>
#include <CG/TimeSmoothedCG.h>

#include "Boundaries/ConstantMassFlowMaserBoundary.h"
#include "../BidispersedChute/BidispersedChute.h"

class AvalancheWithMaserFromRestart : public BidispersedChute
{
public:
    /*! Constructor: construct the MaserForAvalanche simulation, and write one file for all types except restart files;
     * write a new restart file every time. Files are written every 500 time steps.
     */
    AvalancheWithMaserFromRestart(BidispersedChuteParameters parameters) :
            BidispersedChute(parameters),
            maserHasOpened(false)
    {
        setFileType(FileType::ONE_FILE);
        setSaveCount(50. / getTimeStep());
        setMin(0, 0, 0);
        setMax(20, 10, 4);
    }
    
    void setupInitialConditions() override
    {
        BidispersedChute::setupInitialConditions();
    }
    
    /*!
     * A function that overrides DPMBase::actionsBeforeTimeStep; opens the maser at a certain time by extending the base
     * and replacing the periodic boundary with a maser boundary. This also sets a flag that the maser is activated
     */
    void actionsBeforeTimeStep() override
    {
        if (!maserHasOpened && getTime() > 500) //should be 2000
        {
            populateBottomParticles();
            copyParticlesForMaser();
            openMaser();
        }
        if (false)
        {
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
                
                nextTimeToCheck += 1;
                logger(INFO, "New bottom length: %, t = %", bottomLength, getTime());
            }
        }
    }
    
    //this function replaces setupInitialConditions in case of a restart
    //if you have problem-specific variables, set them here
    void actionsOnRestart() override
    {
        setName(getName());
        if (getTime() < 1)
        {
            maserHasOpened = false;
        }
        else
        {
            maserHasOpened = true;
            bottomLength = particleHandler.getHighestPositionX();
            populateBottomParticles();
        }
        setName("postProcessing");
        setTimeMax(getTime() + getTimeStep());
        restartFile.setFileType(FileType::MULTIPLE_FILES);
        restartFile.setSaveCount(1.0/getTimeStep());
        dataFile.setFileType(FileType::ONE_FILE);
        dataFile.setSaveCount(1.0/getTimeStep());
        //fStatFile.setFileType(FileType::ONE_FILE);
        /*auto cg0 = cgHandler.copyAndAddObject(CG<CGCoordinates::XZ,CGFunctions::Lucy>());
        cg0->setWidth(1);
        cg0->statFile.setSaveCount(100);
        cg0->statFile.setName(getName() + ".LucyXZ.stat");
        cg0->setNX(1000);
        cg0->setX(0, 2000);
        cg0->setNZ(20);
        cg0->setZ(-5, 15);*/
        setParticlesWriteVTK(true);
        writePythonFileForVTKVisualisation();
    }

private:
    
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
        //copy BaseParticles up to x = 20
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
        DeletionBoundary deletionBoundary;
        deletionBoundary.set(Vec3D(1, 0, 0), 1e4);
        
        maserHasOpened = true;
        logger(INFO, "opened the maser boundary");
    }
    
    Mdouble bottomLength = 20;
    bool maserHasOpened;
    std::vector<BaseParticle*> bottomParticles;
};

int main(int argc, char* argv[])
{
    //In principle, we restart from a small chute simulation, so all lines except solve can be ignored.
    BidispersedChuteParameters parameters(6, 23, 0.5);
    //set the large particle diameter to 1, and the small particle diameter sizeRatio times smaller
    Mdouble sizeRatio = 1.025;
    parameters.setLargeParticleRadius(0.5);
    parameters.setSmallParticleRadius(0.5/sizeRatio);
    AvalancheWithMaserFromRestart maserForAvalanche(parameters);
    std::stringstream name;
    name << 'H' << parameters.getInflowHeight()
         << "A" << parameters.getAngleInDegrees()
         << "Phi" << (int) 100 * parameters.getConcentrationSmall()
         << "S" << sizeRatio;
    maserForAvalanche.setName(name.str());
    maserForAvalanche.setTimeMax(3000);
    maserForAvalanche.solve(argc,argv);
    return 0;
}
