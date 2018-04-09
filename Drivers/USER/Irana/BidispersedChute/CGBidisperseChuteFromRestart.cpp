#include "StatisticsVector.h"
#include "CG/TimeAveragedCG.h"
#include "BidispersedChute.h"
#include "MercuryTime.h"

///Only use with restart file!
/// Restart the chute with live CG, run for 1 time units and write the files
class CGBidisperseChuteFromRestart : public BidispersedChute
{
public:
    void actionsOnRestart() override
    {
        setName("restarted" + getName());
        setTimeMax(getTime() + 10);
        setTimeStep(1e-5);
        setSaveCount(1./getTimeStep());
        dataFile.setFileType(FileType::ONE_FILE);
        writeRestartFile();
        auto cg0 = cgHandler.copyAndAddObject(TimeAveragedCG<CGCoordinates::Z,CGFunctions::Lucy>());
        cg0->setWidth(0.5);
        cg0->statFile.setSaveCount(100);
        cg0->statFile.setName(getName() + ".LucyZ.stat");
        cg0->setNZ(400);
        cg0->setZ(-5, 35);
        auto cg1 = cgHandler.copyAndAddObject(cg0);
        cg1->setSelectedParticle([](const BaseInteractable* p){ return  p->getSpecies()->getIndex() == 1;});
        cg1->statFile.setName(getName() + ".large.LucyZ.stat");
        auto cg2 = cgHandler.copyAndAddObject(cg0);
        cg2->setSelectedParticle([](const BaseInteractable* p){ return p->getSpecies()->getIndex() == 2;});
        cg2->statFile.setName(getName() + ".small.LucyZ.stat");
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
