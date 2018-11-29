#include "GCGPeriodic.h"
#include "Species/LinearViscoelasticSlidingFrictionReversibleAdhesiveSpecies.h"

/**
 * Simulates the mixing
 */
int main(int argc, char *argv[])
{
    //set all necessary parameters
    Mdouble particleRadius = helpers::readFromCommandLine(argc,argv,"-particleRadius",2.5e-3); // m
    Mdouble particleRadius2 = helpers::readFromCommandLine(argc,argv,"-particleRadius2",1.5e-3); /// ADDED by M
    Mdouble volFracSpec1 = helpers::readFromCommandLine(argc,argv,"volFracSpec1",1.0); /// ADDED by M
    Mdouble polydispersity = 0.05;
    Mdouble polydispersity2 = 0.05; /// ADDED by M
	Mdouble solidFraction = helpers::readFromCommandLine(argc,argv,"-solidFraction",0.1);
    Mdouble rpm = 120;
    auto particleType = SpeciesType::Glass;
    auto wallType = SpeciesType::Glass;

    //setup simulation
    GCGPeriodic gcg(particleType, particleType, wallType, solidFraction, rpm, particleRadius, particleRadius2, volFracSpec1, polydispersity, polydispersity2,15e-3);
    
    gcg.removeOldFiles();
    gcg.setTimeMax(1);
    if (NUMBER_OF_PROCESSORS < 8)
    {
        gcg.setNumberOfDomains({NUMBER_OF_PROCESSORS,1,1});
    }
    else
    {
        unsigned int halfDomains = NUMBER_OF_PROCESSORS / 2.0;
        gcg.setNumberOfDomains({halfDomains,2,1});
    }
    gcg.setParticlesWriteVTK(true);
    gcg.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    gcg.setSaveCount(0.005/gcg.getTimeStep());//save every 0.02s
    gcg.writeParaviewSnapshotScript();
    gcg.solve();
}
