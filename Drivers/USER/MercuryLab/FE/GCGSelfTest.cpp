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
    const Mdouble timeMax = helpers::readFromCommandLine(argc,argv,"-timeMax",0.05);
    // restart if requested, else start new simulation
    if (!restart.empty()) { GCG(restart, timeMax); return 0; }

    // set particle radius (in meters) from command line, e.g. "-particleRadius 0.001" for 1mm particles
    const Mdouble particleRadius = helpers::readFromCommandLine(argc,argv,"-particleRadius",4e-3);
    // set particle radius (in meters) from command line, e.g. "-particleRadius 0.001" for 1mm particles
    const Mdouble particleRadius2 = helpers::readFromCommandLine(argc,argv,"-particleRadius2",2.5e-3); /// ADDED by M
    // set volume fraction of species 1
    const Mdouble volFracSpec1 = helpers::readFromCommandLine(argc,argv,"-volFracSpec1",1.0); /// ADDED by M
    // relative standard deviation of particle radius
    const Mdouble polydispersity = helpers::readFromCommandLine(argc,argv,"-polydispersity",0.1);
    // relative standard deviation of particle radius
    const Mdouble polydispersity2 = helpers::readFromCommandLine(argc,argv,"-polydispersity2",0.1); /// ADDED by M
    // maximum volumetric fill rate [m^3/s] (max, b/c it cannot insert particles if the inflow hopper is full)
    const Mdouble maxFillRate = helpers::readFromCommandLine(argc,argv,"-fillRate",7e-6);
    // blade rotation speed [revolutions per minute]
    const Mdouble rpm = helpers::readFromCommandLine(argc,argv,"-rpm",460);
    // belt velocity [m/s]; 0 means no belt is created
    const Mdouble beltVelocity = helpers::readFromCommandLine(argc,argv,"-beltVelocity",1);
    // timestep between the output files
    const Mdouble timeStepOutput = helpers::readFromCommandLine(argc,argv,"-timeStepOutput",0.005);
    // boolean (0 or 1) to determine if vtk files get written
    const Mdouble writeVTK = helpers::readFromCommandLine(argc,argv,"-writeVTK",true);
    // boolean (0 or 1) to determine if mixer gets initially filled
    const Mdouble initiallyFilled = helpers::readFromCommandLine(argc,argv,"-initiallyFilled",true);
    // boolean (0 or 1) to determine if mixer gets initially filled
    const Mdouble longScrew = helpers::readFromCommandLine(argc,argv,"-longScrew",false);
    // mean residence time [s]; determines how much particles are initially in the blender
    const Mdouble meanResidenceTime = helpers::readFromCommandLine(argc,argv,"-meanResidenceTime",1);
    // casingRadius [m]
    const Mdouble casingRadius = helpers::readFromCommandLine(argc,argv,"-casingRadius",35e-3);
    // blade inclination [deg]
    const Mdouble bladeInclination = helpers::readFromCommandLine(argc,argv,"-bladeInclination",20);
    // type of material
    const auto particleType = SpeciesType::Cocoa;

    // IMPLEMENTATION
    GCG gcg(particleType, maxFillRate, rpm, beltVelocity, particleRadius, particleRadius2, volFracSpec1, polydispersity, polydispersity2, longScrew, casingRadius/35e-3, bladeInclination);
    if (initiallyFilled) gcg.initialFillBasedOnMeanResidenceTime(meanResidenceTime);
    gcg.setName("GCGSelfTest");
    gcg.removeOldFiles();
    gcg.setTimeMax(timeMax);
    gcg.setNumberOfDomains({NUMBER_OF_PROCESSORS,1,1});
    gcg.setParticlesWriteVTK(writeVTK);
    gcg.setWallsWriteVTK(writeVTK);
    gcg.setSaveCount(timeStepOutput/gcg.getTimeStep());

    // SOLVE
    gcg.writeParaviewSnapshotScript();
    helpers::writeCommandLineToFile(gcg.getName()+".cmd", argc, argv);
    gcg.addCG();
    gcg.solve();
    
//    // RESTART
//    if (NUMBER_OF_PROCESSORS==1)
//    {
//        logger(INFO, "Running 20% longer to test the restart capability");
//        std::string cmd = "./GCGSelfTest -r GCGSelfTest -timeMax " + helpers::to_string(1.2 * timeMax, 3);
//        system(cmd.c_str());
//    }
    
    return 0;
}
