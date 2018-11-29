#include "GCG.h"
#include "Species/LinearViscoelasticSlidingFrictionReversibleAdhesiveSpecies.h"

/**
 * Simulates mixing in a GCG-70
 */
int main(int argc, char *argv[])
{
    // PARAMETERS
    // default values are based on one of the experimental settings we received from Fette (case 11):
    //  - volumetric fill rate set to 7e-6 m^3/s, based on throughput of 40kg/h and MCC density of 1.582 kg/l
    //  - initial filling volume set to 0.14 l, based on a mean residence time of 20s and the volumetric fill rate
    //  - rpm of 460, long screw

    // allows passing in a restart file
    const std::string restart = helpers::readFromCommandLine(argc,argv,"-r",std::string(""));
    // simulation time
    const Mdouble timeMax = helpers::readFromCommandLine(argc,argv,"-timeMax",50.0e-4);
    // restart if requested, else start new simulation
    if (!restart.empty()) { GCG(restart, timeMax); return 0; }

    const Mdouble volFracSpec1 = helpers::readFromCommandLine(argc,argv,"-volFracSpec1",0.30); /// ADDED by M
    // relative standard deviation of particle radius
    const Mdouble maxFillRate = helpers::readFromCommandLine(argc,argv,"-fillRate",20e-6); //50 kg/h = 14 g/s = 10e-6 m^3/s
    // blade rotation speed [revolutions per minute]
    const Mdouble rpm = helpers::readFromCommandLine(argc,argv,"-rpm",460.0);
    // belt velocity [m/s]; 0 means no belt is created
    const Mdouble beltVelocity = helpers::readFromCommandLine(argc,argv,"-beltVelocity",1.0);
    // timestep between the output files
    const Mdouble timeStepOutput = helpers::readFromCommandLine(argc,argv,"-timeStepOutput",0.005);
    // boolean (0 or 1) to determine if vtk files get written
    const bool writeVTK = helpers::readFromCommandLine(argc,argv,"-writeVTK",true);
    // boolean (0 or 1) to determine if mixer gets initially filled
    const Mdouble initiallyFilled = helpers::readFromCommandLine(argc,argv,"-initiallyFilled",false);
    // boolean (0 or 1) to determine if mixer gets initially filled
    const bool longScrew = helpers::readFromCommandLine(argc,argv,"-longScrew",false);
    // mean residence time [s]; determines how much particles are initially in the blender
    const Mdouble meanResidenceTime = helpers::readFromCommandLine(argc,argv,"-meanResidenceTime",30.0);
    // casingRadius [m], or scaleFactor
    const Mdouble casingRadius = helpers::readFromCommandLine(argc,argv,"-casingRadius",70e-3);
    const Mdouble scaleFactor = helpers::readFromCommandLine(argc,argv,"-scaleFactor",casingRadius/35e-3);
    // blade inclination [deg]
    const Mdouble bladeInclination = helpers::readFromCommandLine(argc,argv,"-bladeInclination",20.0);
    // type of material (0=PH101--APAP, 1=SD100--MPT)
    const unsigned mixType = helpers::readFromCommandLine(argc,argv,"-mixType",(unsigned)0);
    const auto fillerType = (mixType==0)?SpeciesType::PH101:SpeciesType::SD100;
    const auto apiType = (mixType==0)?SpeciesType::APAPP:SpeciesType::MPT;
    std::cout << std::flush;
    
    // IMPLEMENTATION
    GCG gcg(fillerType, apiType, maxFillRate, rpm, beltVelocity, volFracSpec1, longScrew, scaleFactor, bladeInclination);
    if (initiallyFilled) gcg.initialFillBasedOnMeanResidenceTime(meanResidenceTime);
    gcg.setName("GCGMixtureSelfTest");
    gcg.removeOldFiles();
    gcg.setTimeMax(timeMax);
    gcg.setNumberOfDomains({NUMBER_OF_PROCESSORS,1,1});
    gcg.setParticlesWriteVTK(writeVTK);
    gcg.setWallsWriteVTK(writeVTK);
    gcg.setSaveCount(static_cast<unsigned>(timeStepOutput/gcg.getTimeStep()));

    // SOLVE
    gcg.writeParaviewSnapshotScript();
    helpers::writeCommandLineToFile(gcg.getName()+".cmd", argc, argv);
    gcg.addCG();
    gcg.solve();
    
//    // RESTART
//    if (NUMBER_OF_PROCESSORS==1)
//    {
//        logger(INFO, "Running 20% longer to test the restart capability");
//        std::string cmd = "./GCGSelfTest -r GCGSelfTest -timeMax " + helpers::to_string(1.2 * timeMax, 8);
//        system(cmd.c_str());
//    }
    
    return 0;
}
