#include "GCG.h"
#include "inflowData.h"

enum RunID {
    _, _1, _2, _2b, _3, _4, _5, _6, _8, _11, _14, _17, _19, _19b, _21, _21b, _22, _22b, _30, _82
};

/*!
 * write to file
 */
std::ostream &operator<<(std::ostream &os, RunID runID) {
    switch (runID) {
        case _:
            os << "";
            break;
        case _1:
            os << "1";
            break;
        case _2:
            os << "2";
            break;
        case _2b:
            os << "2b";
            break;
        case _3:
            os << "3";
            break;
        case _4:
            os << "4";
            break;
        case _5:
            os << "5";
            break;
        case _6:
            os << "6";
            break;
        case _8:
            os << "8";
            break;
        case _11:
            os << "11";
            break;
        case _14:
            os << "14";
            break;
        case _17:
            os << "17";
            break;
        case _19:
            os << "19";
            break;
        case _19b:
            os << "19b";
            break;
        case _21:
            os << "21";
            break;
        case _21b:
            os << "21b";
            break;
        case _22:
            os << "22";
            break;
        case _22b:
            os << "22b";
            break;
        case _30:
            os << "30";
            break;
        case _82:
            os << "82";
            break;
        default:
            logger(ERROR, "runID % not found", static_cast<unsigned>(runID));
    }
    return os;
}

/*!
 * write from file
 */
RunID getRunID(std::string s) {
    RunID runID;
    if (s == "") runID = _;
    else if (s == "1") runID = _1;
    else if (s == "2") runID = _2;
    else if (s == "2b") runID = _2b;
    else if (s == "3") runID = _3;
    else if (s == "4") runID = _4;
    else if (s == "5") runID = _5;
    else if (s == "6") runID = _6;
    else if (s == "8") runID = _8;
    else if (s == "11") runID = _11;
    else if (s == "14") runID = _14;
    else if (s == "17") runID = _17;
    else if (s == "19") runID = _19;
    else if (s == "19b") runID = _19b;
    else if (s == "21") runID = _21;
    else if (s == "21b") runID = _21b;
    else if (s == "22") runID = _22;
    else if (s == "22b") runID = _22b;
    else if (s == "30") runID = _30;
    else if (s == "82") runID = _82;
    else logger(ERROR, "SpeciesType % not found", s);
    return runID;
}


/**
 * Simulates mixing in a GCG-70
 */
int main(int argc, char *argv[]) {
    // GENERAL PARAMETERS

    // read simulation time (default 60)
    Mdouble timeMax = helpers::readFromCommandLine(argc, argv, "-timeMax", 60.0);

    // restart from restart file if requested; else start new simulation
    const std::string restart = helpers::readFromCommandLine(argc, argv, "-r", std::string(""));
    if (!restart.empty()) {
        GCG(restart, timeMax);
        return 0;
    }

    // belt velocity [m/s]; 0 means no belt is created
    const Mdouble beltVelocity = helpers::readFromCommandLine(argc, argv, "-beltVelocity", 1.0);

    // timestep between the output files
    // use 0.005 for video output
    const Mdouble timeStepOutput = helpers::readFromCommandLine(argc, argv, "-timeStepOutput", 0.005);

    // Determines if vtk files get written
    const bool writeVTK = helpers::readFromCommandLine(argc, argv, "-writeVTK", false);

    // Determines if vtk files get written
    Mdouble shaftRadius = 15e-3;

    // Determines if mixer gets initially filled
    const Mdouble initiallyFilled = helpers::readFromCommandLine(argc, argv, "-initiallyFilled", false);
    // Mean residence time [s]; determines how much particles are initially in the blender
    const Mdouble meanResidenceTime = initiallyFilled ?
            helpers::readFromCommandLine(argc, argv, "-meanResidenceTime",30.0) : 0.0;

    // FIXED SIMULATION PARAMETERS

    // Scale of casing radius vs default value of 35mm
    Mdouble scaleFactor = 1;

    // blade inclination [deg]
    Mdouble bladeInclination = 20.0;

    // VARIABLE SIMULATION PARAMETERS

    // type of filler and api (determines density, contact properties and psd)
    SpeciesType fillerType, apiType;
    // solid volume fraction of filler (solid volume fraction of api is 1-fillerVolumeFraction)
    Mdouble fillerVolumeFraction;
    // relative standard deviation of particle radius
    Mdouble maxFillRate;
    // blade rotation speed [revolutions per minute]
    Mdouble rpm;
    // blade pattern (default: OldScrew)
    BladePattern bladePattern;
    // variable flow rate and sampling interval for filler and api
    std::vector<Mdouble> fillerFlowRate;
    Mdouble fillerSamplingInterval;
    std::vector<Mdouble> apiFlowRate;
    Mdouble apiSamplingInterval;



    // runID
    RunID runID = getRunID(helpers::readFromCommandLine(argc, argv, "-runID", std::string("14")));

    if (runID == _) { //if no runID is selected use values from command line
        // type of material (0=PH101--APAP, 1=SD100--MPT)
        const unsigned mixType = helpers::readFromCommandLine(argc, argv, "-mixType", (unsigned) 0);
        fillerType = (mixType == 0) ? SpeciesType::PH101 : SpeciesType::SD100;
        apiType = (mixType == 0) ? SpeciesType::APAPP : SpeciesType::MPT;
        fillerVolumeFraction = helpers::readFromCommandLine(argc, argv, "-fillerVolumeFraction", 0.99);
        maxFillRate = helpers::readFromCommandLine(argc, argv, "-fillRate", 20e-6);
        rpm = helpers::readFromCommandLine(argc, argv, "-rpm", 460.0);
        bladePattern = (BladePattern) helpers::readFromCommandLine(argc, argv, "-bladePattern", 1);
    } else {
        logger(INFO, "runID %", runID);

        //select fillerType, apiType, fillerVolumeFraction, and check if runID exists
        //also set variable flow rate and sampling interval for filler and api
        switch (runID) {
            case _1:
            case _2:
            case _2b:
            case _3:
            case _8:
            case _14:
            case _19:
            case _21:
            case _22:
                fillerType = SpeciesType::SD100;
                apiType = SpeciesType::MPT;
                fillerVolumeFraction = 0.99;
                if (runID == _2b) {
                    fillerFlowRate = pearlitol24bWeight;
                    fillerSamplingInterval = pearlitol24bSamplingInterval;
                    apiFlowRate = MPT0_24bWeight;
                    apiSamplingInterval = MPT0_24bSamplingInterval;
                } else {
                    fillerFlowRate = pearlitol12Weight;
                    fillerSamplingInterval = pearlitol12SamplingInterval;
                    apiFlowRate = MPT0_12Weight;
                    apiSamplingInterval = MPT0_12SamplingInterval;
                }
                break;
            case _4:
            case _5:
            case _6:
            case _11:
            case _17:
            case _19b:
            case _21b:
            case _22b:
            case _30:
                if (runID==_30) shaftRadius = 10e-3;
                fillerType = SpeciesType::PH101;
                apiType = SpeciesType::APAPP;
                fillerVolumeFraction = 0.3;
                fillerFlowRate = avicel3_6Weight;
                fillerSamplingInterval = avicel3_6SamplingInterval;
                apiFlowRate = apap8_4Weight;
                apiSamplingInterval = apap8_4SamplingInterval;
                break;
            case _82:
                fillerType = SpeciesType::MPT;
                apiType = fillerType;
                fillerVolumeFraction = 0.5;
                fillerFlowRate = MPT2Weight;
                fillerSamplingInterval = MPT2SamplingInterval;
                apiFlowRate = fillerFlowRate;
                apiSamplingInterval = fillerSamplingInterval;
                break;
            default:
                logger(ERROR, "RunID % does not exist", runID);
        }
        logger(INFO, "fillerType %, apiType %, fillerVolumeFraction %", fillerType, apiType, fillerVolumeFraction);

        //turn mass flow rate into volume flow rate
        Mdouble fillerDensity = getDensity(fillerType);
        for (auto& v : fillerFlowRate) {
            v /= fillerDensity;
        }
        Mdouble apiDensity = getDensity(apiType);
        for (auto& v : apiFlowRate) {
            v /= apiDensity;
        }

        //select maxFillRate
        const Mdouble massFlow = (runID == _2b) ? 24 : 12;
        const Mdouble meanDensity = fillerVolumeFraction * fillerDensity
                + (1.0 - fillerVolumeFraction) * apiDensity;
        maxFillRate = massFlow / 3600. / meanDensity;
        logger(INFO, "massFlow % kg/h, meanDensity % kg/m^3, maxFillRate % ml/s", massFlow, meanDensity, maxFillRate*1e6);


        //select rpm
        switch (runID) {
            case _1:
            case _4:
            case _82:
                rpm = 112;
                break;
            case _3:
            case _6:
                rpm = 460;
                break;
            default:
                rpm = 286;
        }
        logger(INFO, "rpm %", rpm);

        //select timeMax
        switch (runID) {
            case _1:
            case _2:
            case _3:
            case _22:
                timeMax = 120;
                break;
            default:
                timeMax = 60;
        }
        logger(INFO, "timeMax %", timeMax);

        //select bladePattern
        switch (runID) {
            case _8:
            case _11:
                bladePattern = BladePattern::FullMixingWithForward; //3F
                break;
            case _14:
            case _17:
                bladePattern = BladePattern::FullMixingWithForwardWithPins; //3F*
                break;
            case _19:
            case _19b:
                bladePattern = BladePattern::FullAltShifted; //11
                break;
            case _21:
            case _21b:
                bladePattern = BladePattern::FullAlt; //6
                break;
            case _22:
            case _22b:
                bladePattern = BladePattern::FullMixingShifted; //8
                break;
            default:
                bladePattern = BladePattern::FullMixing; //3
        }
        logger(INFO, "bladePattern %", bladePattern);

        //allows the user to reset certain values
        //fillerVolumeFraction = helpers::readFromCommandLine(argc, argv, "-fillerVolumeFraction", fillerVolumeFraction);
        //maxFillRate = helpers::readFromCommandLine(argc, argv, "-fillRate", maxFillRate*1e6)*1e-6;
        rpm = helpers::readFromCommandLine(argc, argv, "-rpm", rpm);
    }

    // IMPLEMENTATION
    logger(INFO, "\n\nCalling constructor ...");
    GCG gcg(fillerType, apiType, maxFillRate, rpm, beltVelocity, fillerVolumeFraction, bladePattern,
            scaleFactor, bladeInclination, shaftRadius);
    if (fillerFlowRate.size()!=0) gcg.useVariableFillRates(fillerFlowRate, fillerSamplingInterval, apiFlowRate, apiSamplingInterval);
    if (initiallyFilled) gcg.initialFillBasedOnMeanResidenceTime(meanResidenceTime);
    //logger(INFO, "\n\nSystem length %",gcg.getXMax()-gcg.getXMin());
    std::stringstream name;
    name << "GCG" << runID;
    gcg.setName(name.str());
    logger(INFO, "\n\nRemoving old files ...");
    gcg.removeOldFiles();
    gcg.setTimeMax(timeMax);
    gcg.setNumberOfDomains({NUMBER_OF_PROCESSORS, 1, 1});
    gcg.setParticlesWriteVTK(writeVTK);
    gcg.setWallsWriteVTK(true);
    gcg.fStatFile.writeFirstAndLastTimeStep();
    gcg.setSaveCount(static_cast<unsigned>(timeStepOutput / gcg.getTimeStep()));
    gcg.restartFile.setSaveCount(static_cast<unsigned>(1.0 / gcg.getTimeStep()));

    // SOLVE
    gcg.writeParaviewSnapshotScript();
    helpers::writeCommandLineToFile(gcg.getName() + ".cmd", argc, argv);
    //gcg.addCG();

    auto filler = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(gcg.speciesHandler.getObject(0));
    auto api = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(gcg.speciesHandler.getObject(1));
    auto wall = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(gcg.speciesHandler.getObject(2));
    auto fillerApi = dynamic_cast<LinearPlasticViscoelasticFrictionMixedSpecies*>(gcg.speciesHandler.getMixedObject(0,1));
    auto fillerWall = dynamic_cast<LinearPlasticViscoelasticFrictionMixedSpecies*>(gcg.speciesHandler.getMixedObject(0,2));
    auto apiWall = dynamic_cast<LinearPlasticViscoelasticFrictionMixedSpecies*>(gcg.speciesHandler.getMixedObject(1,2));

    logger(INFO, "Filler (%)\tmus % \tkc/k* % ",fillerType,filler->getSlidingFrictionCoefficient(),
            filler->getCohesionStiffness()/filler->getUnloadingStiffnessMax());
    logger(INFO, "API (%)\tmus % \tkc/k* % ",apiType,api->getSlidingFrictionCoefficient(),
            api->getCohesionStiffness()/api->getUnloadingStiffnessMax());
    logger(INFO, "Filler-API  \tmus % \tkc/k* % ",fillerApi->getSlidingFrictionCoefficient(),
            fillerApi->getCohesionStiffness()/fillerApi->getUnloadingStiffnessMax());
    logger(INFO, "Filler-Wall \tmus % \tkc/k* % ",fillerWall->getSlidingFrictionCoefficient(),
            fillerWall->getCohesionStiffness()/fillerWall->getUnloadingStiffnessMax());
    logger(INFO, "API-Wall    \tmus % \tkc/k* % ",apiWall->getSlidingFrictionCoefficient(),
            apiWall->getCohesionStiffness()/apiWall->getUnloadingStiffnessMax());

    logger(INFO, "\n\nCalling solve ...");
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
