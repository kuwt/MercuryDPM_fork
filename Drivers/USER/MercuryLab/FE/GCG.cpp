#include "GCG.h"
#include "Species/LinearViscoelasticSlidingFrictionReversibleAdhesiveSpecies.h"

using helpers::to_string;

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

    // set particle radius (in meters) from command line, e.g. "-particleRadius 0.001" for 1mm particles
    const Mdouble particleRadius = helpers::readFromCommandLine(argc,argv,"-particleRadius",2e-3);
    // set particle radius (in meters) from command line, e.g. "-particleRadius 0.001" for 1mm particles
    const Mdouble particleRadius2 = helpers::readFromCommandLine(argc,argv,"-particleRadius2",1.5e-3); /// ADDED by M
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
    const Mdouble beltVelocity = helpers::readFromCommandLine(argc,argv,"-beltVelocity",0.1);
    // simulation time
    const Mdouble timeMax = helpers::readFromCommandLine(argc,argv,"-timeMax",1.0);
    // timestep between the output files
    const Mdouble timeStepOutput = helpers::readFromCommandLine(argc,argv,"-timeStepOutput",0.005);
    // boolean (0 or 1) to determine if vtk files get written
    const bool writeVTK = helpers::readFromCommandLine(argc,argv,"-writeVTK",true);
    // boolean (0 or 1) to determine if mixer gets initially filled
    const bool initiallyFilled = helpers::readFromCommandLine(argc,argv,"-initiallyFilled",false);
    // boolean (0 or 1) to determine if mixer gets initially filled
    const bool longScrew = helpers::readFromCommandLine(argc,argv,"-longScrew",true);
    // mean residence time [s]; determines how much particles are initially in the blender
    const Mdouble meanResidenceTime = helpers::readFromCommandLine(argc,argv,"-meanResidenceTime",20);
    // casingRadius [m]
    const Mdouble casingRadius = helpers::readFromCommandLine(argc,argv,"-casingRadius",35e-3);
    // blade inclination [deg]
    const Mdouble bladeInclination = helpers::readFromCommandLine(argc,argv,"-bladeInclination",20);
    // type of material
    const auto particleType = SpeciesType::Cocoa;

    // how to name a simulation (change as needed)
    int r = particleRadius*1e6;
    int r2 = particleRadius2*1e6;
    std::string name = "GCG_R" + to_string<int>(particleRadius*1e6)
                       + (volFracSpec1==1.0?"":
                         ("_" + to_string<int>(particleRadius2*1e6) + "_VF" + to_string<int>(volFracSpec1*100)))
                       + (initiallyFilled?"_Filled":"")
                       + "_NP" + to_string<int>(NUMBER_OF_PROCESSORS)
                       + "_RPM" + to_string<int>(rpm)
                       + "_FR" + to_string<int>(maxFillRate*1e9)
                       + (longScrew?"":"_Short")
                       + (casingRadius==35e-3?"":"_CR"+ to_string<int>(casingRadius*1e3))
                       + (bladeInclination==20?"":"_BI"+ to_string<int>(bladeInclination))
                       ;
    logger(INFO,"Simulation data will be outputted to %.* files",name);

    // IMPLEMENTATION
    GCG gcg(particleType, maxFillRate, rpm, beltVelocity, particleRadius, particleRadius2, volFracSpec1, polydispersity, polydispersity2, longScrew, casingRadius/35e-3, bladeInclination);
    if (initiallyFilled) gcg.initialFillBasedOnMeanResidenceTime(meanResidenceTime);
    gcg.setName(name);
    gcg.removeOldFiles();
    gcg.setTimeMax(timeMax);
    if (NUMBER_OF_PROCESSORS < 100)
    {
        gcg.setNumberOfDomains({NUMBER_OF_PROCESSORS,1,1});
    }
    else
    {
        unsigned int halfDomains = NUMBER_OF_PROCESSORS / 2.0;
        gcg.setNumberOfDomains({halfDomains,2,1});
    }
    gcg.setParticlesWriteVTK(writeVTK);
    gcg.setWallsWriteVTK(writeVTK);
    gcg.setSaveCount(timeStepOutput/gcg.getTimeStep());

    // SOLVE
    gcg.addCG();
    gcg.writeParaviewSnapshotScript();
    helpers::writeCommandLineToFile(gcg.getName()+".cmd", argc, argv);
    gcg.solve();

    //realistic settings (case 11)
    // fill rate: 40/1.427/3600
    //./GCG -rpm 460 -fillRate 0.0000078 -meanResidenceTime 20 -initiallyFilled 1
}
