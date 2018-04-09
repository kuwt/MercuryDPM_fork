//
// Created by irana on 6/28/17.
//


#include <Boundaries/DeletionBoundary.h>
#include "Boundaries/ConstantMassFlowMaserBoundary.h"
#include "../BidispersedChute/BidispersedChute.h"

class RestartForBScStudents : public BidispersedChute
{
public:
    RestartForBScStudents()
    {
        readRestartFile("MonodispersedAvalanche");
        setTimeMax(1000);
        setXMax(1500);
        setZMax(100);
        setName("New" + getName());
        dataFile.setFileType(FileType::ONE_FILE);
        restartFile.setFileType(FileType::ONE_FILE);
        fStatFile.setFileType(FileType::ONE_FILE);
        eneFile.setFileType(FileType::ONE_FILE);
        for (unsigned int i = 0; i < boundaryHandler.getSize(); ++i)
        {
            if (boundaryHandler.getObject(i)->getName() == "ConstantMassFlowMaserBoundary")
            {
                maserBoundary = static_cast<ConstantMassFlowMaserBoundary*>(boundaryHandler.getObject(i));
                logger(INFO, "Found maser boundary");
            }
        }
        for (BaseParticle* p : particleHandler)
        {
            if (p->getPosition().X > maserBoundary->getDistanceLeft() &&
                p->getPosition().X < maserBoundary->getDistanceRight() &&
                p->isFixed())
            {
                BaseParticle* pCopy = p->copy();
                logger(DEBUG, "species of bottom particle: %", p->getSpecies()->getIndex());
                maserBoundary->removeParticleFromMaser(pCopy);
                baseParticles.push_back(pCopy);
            }
        }
        deletionBoundaryX = boundaryHandler.copyAndAddObject(DeletionBoundary());
        deletionBoundaryX->set({1,0,0}, getXMax());
        deletionBoundaryZ = boundaryHandler.copyAndAddObject(DeletionBoundary());
        deletionBoundaryZ->set({0,0,1}, getZMax());
    }
    
    void actionsAfterTimeStep() override
    {
        BidispersedChute::actionsAfterTimeStep();
        static Mdouble timeToExtend = 380;
        static Mdouble distanceToExtend = 1800;
        if (getTime() > timeToExtend)
        {
            for (BaseParticle* b : baseParticles)
            {
                while (b->getPosition().X < distanceToExtend)
                {
                    particleHandler.copyAndAddObject(b);
                    b->setPosition(b->getPosition() + Vec3D(20,0,0));
                }
            }
            deletionBoundaryX->set({1,0,0}, distanceToExtend);
            setXMax(distanceToExtend);
            timeToExtend += 20;
            distanceToExtend += 200;
        }
    }

private:
    void deleteHighXParticles()
    {
        for (auto it = particleHandler.end() - 1; it >= particleHandler.begin(); it--)
        {
            if ((*it)->getPosition().X > getXMax() || (*it)->getPosition().Z > getZMax())
            {
                particleHandler.removeObject((*it)->getIndex());
            }
        }
    }
    
    ConstantMassFlowMaserBoundary* maserBoundary;
    std::vector<BaseParticle*> baseParticles;
    DeletionBoundary* deletionBoundaryX;
    DeletionBoundary* deletionBoundaryZ;
};

int main()
{
    RestartForBScStudents problem;
    problem.setSaveCount(1/problem.getTimeStep());
    problem.setLastSavedTimeStep(NEVER);
    problem.solve();
}
