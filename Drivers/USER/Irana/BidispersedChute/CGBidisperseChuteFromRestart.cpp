#include "StatisticsVector.h"
#include "CG/TimeAveragedCG.h"
#include "BidispersedChute.h"
#include "MercuryTime.h"

/// Only use with restart file! ./CGBidisperseChuteFromRestart -r 'MySerialRestartFile.restart'
/// Restart the chute with live CG
class CGBidisperseChuteFromRestart : public BidispersedChute
{
public:
    void actionsOnRestart() override
    {
        setName("restarted" + getName());
        const Mdouble scaleFactor = 25;
        for(BaseParticle* p : particleHandler)
        {
            if (!p->isFixed())
            {
                const Vec3D position = p->getPosition();
                p->setPosition(Vec3D(position.X, position.Y, scaleFactor * position.Z));
            }
        }
        setTimeMax(getTime() + getTimeStep());
        dataFile.setFileType(FileType::ONE_FILE);
        //setParticlesWriteVTK(true); //uncomment for paraview-files
        //writePythonFileForVTKVisualisation(); //uncomment for paraview-files
        //setName("postprocessing");//paraview-files cannot have a . in the filename, so either rename + sed or uncomment
        auto cg0 = cgHandler.copyAndAddObject(TimeAveragedCG<CGCoordinates::XZ,CGFunctions::Lucy>());
        cg0->setWidth(scaleFactor / 2);
        cg0->statFile.setSaveCount(10000);
        cg0->statFile.setName(getName() + ".LucyXZ.stat");
        cg0->setNZ(50);
        cg0->setNX(50);
        cg0->setX(0, 1000);
        cg0->setZ(-1 * scaleFactor, 9 * scaleFactor);
        auto cg1 = cgHandler.copyAndAddObject(cg0);
        cg1->setSelectedParticle([](const BaseInteractable* p){ return  p->getSpecies()->getIndex() == 1;});
        cg1->statFile.setName(getName() + ".large.LucyXZ.stat");
        auto cg2 = cgHandler.copyAndAddObject(cg0);
        cg2->setSelectedParticle([](const BaseInteractable* p){ return p->getSpecies()->getIndex() == 2;});
        cg2->statFile.setName(getName() + ".small.LucyXZ.stat");
        setParticlesWriteVTK(false);
    }
};

int main(int argc, char* argv[])
{
    CGBidisperseChuteFromRestart problem;
    Time time;
    time.tic();
    problem.solve(argc, argv);
    std::cout << "Time needed for CG: " << time.toc() << std::endl;
    return 0;
}
