//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <Mercury3D.h>
#include <CG/CG.h>
#include <CG/TimeSmoothedCG.h>
#include <CG/TimeAveragedCG.h>
#include <CG/Functions/Linear.h>
#include <CG/Fields/GradVelocityField.h>
#include <CG/Fields/LiquidMigrationFields.h>
#include <algorithm>
#include <Interactions/AdhesiveForceInteractions/LiquidMigrationWilletInteraction.h>

template<template<class, template<class> class, class> class CGType,
        template<class> class BaseFunction, class Fields>
void addObject(CGHandler &cg, std::string coordinate) {
    if (coordinate == "O") {
        cg.copyAndAddObject(CGType<CGCoordinates::O, BaseFunction, Fields>());
    } else if (coordinate == "X") {
        cg.copyAndAddObject(CGType<CGCoordinates::X, BaseFunction, Fields>());
    } else if (coordinate == "Y") {
        cg.copyAndAddObject(CGType<CGCoordinates::Y, BaseFunction, Fields>());
    } else if (coordinate == "Z") {
        cg.copyAndAddObject(CGType<CGCoordinates::Z, BaseFunction, Fields>());
    } else if (coordinate == "XY") {
        cg.copyAndAddObject(CGType<CGCoordinates::XY, BaseFunction, Fields>());
    } else if (coordinate == "XZ") {
        cg.copyAndAddObject(CGType<CGCoordinates::XZ, BaseFunction, Fields>());
    } else if (coordinate == "YZ") {
        cg.copyAndAddObject(CGType<CGCoordinates::YZ, BaseFunction, Fields>());
    } else if (coordinate == "XYZ") {
        cg.copyAndAddObject(CGType<CGCoordinates::XYZ, BaseFunction, Fields>());
    } else {
        logger(ERROR, "CGCoordinate % not understood; options are O,X,Y,Z,XY,XZ,YZ,XYZ", coordinate);
    }
}

template<template<class> class BaseFunction, class Fields>
void addObject(CGHandler &cg, std::string typeOrig, std::string coordinate) {
    std::string type = typeOrig;
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    if (type == "cg") {
        addObject<CG, BaseFunction, Fields>(cg, coordinate);
    } else if (type == "timesmoothedcg") {
        addObject<TimeSmoothedCG, BaseFunction, Fields>(cg, coordinate);
    } else if (type == "timeaveragedcg") {
        addObject<TimeAveragedCG, BaseFunction, Fields>(cg, coordinate);
    } else {
        logger(ERROR, "CGType % not understood; options are CG, TimeSmoothedCG, TimeAveragedCG", typeOrig);
    }
}

template<class Fields>
void addObject(CGHandler &cg, std::string type, std::string coordinate, std::string function) {
    if (function == "Gauss") {
        addObject<CGFunctions::Gauss, Fields>(cg, type, coordinate);
    } else if (function == "Lucy") {
        addObject<CGFunctions::Lucy, Fields>(cg, type, coordinate);
    } else if (function == "Linear") {
        addObject<CGFunctions::Linear, Fields>(cg, type, coordinate);
    } else if (function == "Heaviside") {
        addObject<CGFunctions::Heaviside, Fields>(cg, type, coordinate);
    } else {
        logger(ERROR, "CGFunction % not understood; options are Gauss, Lucy, Linear, Heaviside", function);
    }
}

/**
 * Creates a CG-object with the template parameters and parent class specified by string arguments.
 * The object is added to the cgHandler and a BaseCG pointer to it is returned.
 * To obtain the right template parameters from the string arguments, a cascade of 4 functions is used (see above).
 * Each function expresses one of the strings to a tempate parameter.
 */
BaseCG *addObject(CGHandler &cg, std::string type, std::string coordinate, std::string function, std::string fieldsOrig) {
    //change strings to lower case
    std::string fields = fieldsOrig;
    std::transform(fields.begin(), fields.end(), fields.begin(), ::tolower);
    if (fields == "standardfields" || fields == "standard") {
        addObject<CGFields::StandardFields>(cg, type, coordinate, function);
    } else if (fields == "gradvelocityfield" || fields == "gradvelocity") {
        addObject<CGFields::GradVelocityField>(cg, type, coordinate, function);
    } else if (fields == "liquidmigrationfields" || fields == "liquidmigration") {
        addObject<CGFields::LiquidMigrationFields>(cg, type, coordinate, function);
    } else {
        logger(ERROR, "CGFields % not understood; options are standard, gradVelocity, liquidMigration ", fieldsOrig);
    }
    return cg.getLastObject();
}

/**
 * Creates an object in the cgHandler of the dpm class that satisfies the options specified by the command line arguments.
 *
 * First, it checks if -help was specified. In this case, the help file is displayed, and the program exits.
 *
 * Second, the correct type of CG object is determined, based on the following command line arguments:
 *   name, -coordinate, -function, -fields, -timeaverage, -timesmooth
 *
 * Third, all other command line arguments are read and passed into the CG object.
 */
void commandLineCG(Mercury3D &dpm, int argc, char **argv)
{
    //change options to lower case
    for (unsigned i = 2; i < argc; i++) {
        if (argv[i][0] == '-')
            for (unsigned j = 1; j < strlen(argv[i]); ++j) {
                argv[i][j] = tolower(argv[i][j]);
            }
    }

    //checks if if -help option was requested
    for (unsigned i = 1; i < argc; i++)
        if (!strcmp(argv[i], "-help")) {
            helpers::more(std::string(argv[0]) + ".txt");
            exit(-1);
        }

    //coordinate
    //Stores the CGCoordinate type.
    //This variable is a template parameter of CG and needs to be specified by the argument -stattype
    //If unspecified, it is set to default value O.
    std::string coordinate = "O";
    for (unsigned i = 1; i < argc; i++)
        if (!strcmp(argv[i], "-coordinates") || !strcmp(argv[i], "-stattype")) {
            coordinate = argv[i + 1];
            logger(INFO, "Set CGCoordinates to %", coordinate);
        }

    //function
    //Stores the CGFunction type.
    //This variable is a template parameter of CG and needs to be specified by the argument -cgtype
    //If unspecified, it is set to default value Lucy.
    std::string function = "Lucy";
    for (unsigned i = 1; i < argc; i++)
        if (!strcmp(argv[i], "-function") || !strcmp(argv[i], "-cgtype")) {
            function = argv[i + 1];
            logger(INFO, "Set CGFunction to %", function);
        }

    //fields
    //Stores the CGField type.
    //This variable is a template parameter of CG and needs to be specified by the argument -function
    //If unspecified, it is set to default value Standard.
    std::string fields = "StandardFields";
    for (unsigned i = 1; i < argc; i++)
        if (!strcmp(argv[i], "-fields")) {
            fields = argv[i + 1];
            logger(INFO, "Set CGFields to %", fields);
        }

    //type
    //Determines the CG type (CG or TimeAveragedCG)
    //It needs to be specified by the argument -timeAverage
    //If unspecified, it is set to default value 'false'.
    // take the 'T' off the coordinate name if time averaged
    std::string type = "CG";
    if (coordinate[0] == 'T') {
        type = "TimeSmoothedCG";
        coordinate = coordinate.substr(1);
    }
    for (unsigned i = 1; i < argc; i++)
        if (!strcmp(argv[i], "-timeaverage")) {
            type = "TimeAveragedCG";
            logger(INFO, "Creating time-averaged CG");
        } else if (!strcmp(argv[i], "-timeaveraging")) {
            logger(ERROR, "% is not a valid argument; use -timeaverage instead",argv[i]);
        }
    for (unsigned i = 1; i < argc; i++)
        if (!strcmp(argv[i], "-timesmooth")) {
            type = "TimeSmoothedCG";
            logger(INFO, "Activating time smoothing");
        }

    //Now create the right CG object and a pointer to it
    BaseCG *cg = addObject(dpm.cgHandler, type, coordinate, function, fields);

    //Reads the name of the file name to be restarted.
    //This variable needs to be specified as the second argument ("./MercuryCG name ...")
    //If unspecified, it is set to default value "Chain".
    std::string name = "";
    if (argc > 1 && argv[1][0] != '-') {
        name = argv[1];
        //logger(INFO, "Evaluating files %.*", name);
    }
    if (name.empty()) {
        logger(ERROR, "Please enter a base name base name of the files to be analysed.\n"
                      "For more information, enter './MercuryCG -help'");
    }

    //restarting from file
    dpm.cgHandler.restart(name);

    //set a default stat file name
    cg->statFile.setName(dpm.getName() + ".stat");

    //set min/max, so h can be set
    cg->setMin(dpm.getMin());
    cg->setMax(dpm.getMax());

    // Now all other arguments are read
    // find first argument that begins with '-'
    unsigned i = 1;
    while (i < argc && argv[i][0] != '-') i++;
    //logger(INFO, "Reading argument % of %", i, argc);
    // interpret arguments
    for (; i < argc; i += 2) {
        if (!strcmp(argv[i], "-w") || !strcmp(argv[i], "-width")) {
            cg->setWidth(atof(argv[i + 1]));
            logger(INFO, "Set cg width to %", cg->getWidth());
        } else if (!strcmp(argv[i], "-std")) {
            cg->setStandardDeviation(atof(argv[i + 1]));
            logger(INFO, "Set cg width to % (std %)", cg->getWidth(),atof(argv[i + 1]));
        } else if (!strcmp(argv[i], "-n")) {
            cg->setN(atoi(argv[i + 1]));
            logger(INFO, "Set n to %", argv[i + 1]);
        } else if (!strcmp(argv[i], "-nx")) {
            cg->setNX(atoi(argv[i + 1]));
            logger(INFO, "Set nx to %", argv[i + 1]);
        } else if (!strcmp(argv[i], "-ny")) {
            cg->setNY(atoi(argv[i + 1]));
            logger(INFO, "Set ny to %", argv[i + 1]);
        } else if (!strcmp(argv[i], "-nz")) {
            cg->setNZ(atoi(argv[i + 1]));
            logger(INFO, "Set nz to %", argv[i + 1]);
        } else if (!strcmp(argv[i], "-h")) {
            cg->setH(atof(argv[i + 1]));
            logger(INFO, "Set n to %x%x% for h=%", cg->getNX(), cg->getNY(), cg->getNZ(), argv[i + 1]);
        } else if (!strcmp(argv[i], "-hx")) {
            cg->setHX(atof(argv[i + 1]));
            logger(INFO, "Set nx to % to satisfy hx=%", cg->getNX(), argv[i + 1]);
        } else if (!strcmp(argv[i], "-hy")) {
            cg->setHY(atof(argv[i + 1]));
            logger(INFO, "Set ny to % to satisfy hy=%", cg->getNY(), argv[i + 1]);
        } else if (!strcmp(argv[i], "-hz")) {
            cg->setHZ(atof(argv[i + 1]));
            logger(INFO, "Set nz to % to satisfy hz=%", cg->getNZ(), argv[i + 1]);
        } else if (!strcmp(argv[i], "-x")) {
            cg->setX(atof(argv[i + 1]), atof(argv[i + 2]));
            logger(INFO, "Set x to (%,%)", argv[i + 1], argv[i + 2]);
            ++i;
        } else if (!strcmp(argv[i], "-y")) {
            cg->setY(atof(argv[i + 1]), atof(argv[i + 2]));
            logger(INFO, "Set y to (%,%)", argv[i + 1], argv[i + 2]);
            ++i;
        } else if (!strcmp(argv[i], "-z")) {
            cg->setZ(atof(argv[i + 1]), atof(argv[i + 2]));
            logger(INFO, "Set z to (%,%)", argv[i + 1], argv[i + 2]);
            ++i;
        } else if (!strcmp(argv[i], "-t")) {
            cg->setTimeMin(atof(argv[i + 1]));
            cg->setTimeMax(atof(argv[i + 2]));
            logger(INFO, "Set t to (%,%)", argv[i + 1], argv[i + 2]);
            i++;
        } else if (!strcmp(argv[i], "-tmin")) {
            cg->setTimeMin(atof(argv[i + 1]));
            logger(INFO, "Set tMin to %", argv[i + 1]);
        } else if (!strcmp(argv[i], "-timemin")) {
            logger(ERROR, "% is not a valid argument; use -tMin instead",argv[i]);
        } else if (!strcmp(argv[i], "-tmax")) {
            cg->setTimeMax(atof(argv[i + 1]));
            logger(INFO, "Set tMax to %", argv[i + 1]);
        } else if (!strcmp(argv[i], "-timemax")) {
            logger(ERROR, "% is not a valid argument; use -tMax instead",argv[i]);
        } else if (!strcmp(argv[i], "-o")) {
            cg->statFile.setName(argv[i + 1]);
            logger(INFO, "Set output file name to %", argv[i + 1]);
        } else if (!strcmp(argv[i], "-coordinates") || !strcmp(argv[i], "-function") || !strcmp(argv[i], "-fields")) {
        } else if (!strcmp(argv[i], "-help") || !strcmp(argv[i], "-restart") || !strcmp(argv[i], "-timeaverage") || !strcmp(argv[i], "-timesmooth")) {
            --i;
        } else {
            logger(ERROR, "Could not read argument %", argv[i]);
        }
    }

    logger(INFO, "Created object of type %", cg->getName());

    bool useDataFiles = true;
    for (unsigned i = 1; i < argc; i++)
        if (!strcmp(argv[i], "-restart")) {
            useDataFiles = false;
            logger(INFO, "Using restart files instead of data files");
        }
    if (useDataFiles) {
        if (!dpm.cgHandler.evaluateDataFiles()) {
            logger(ERROR,"Evaluation of data files has failed. Check input files.");
        }
    } else {
        if (dpm.restartFile.getFileType()==FileType::ONE_FILE) {
            logger(ERROR,"Evaluation of restart files has failed. Check input files.");
        }
        dpm.cgHandler.evaluateRestartFiles();
    }
    logger(INFO, "\n"
                 "MercuryCG has finished.\n"
                 "Coarse-grained output is written to %.\n"
                 "To load output into Matlab, use data=readMercuryCG('%')",
           dpm.cgHandler.getLastObject()->statFile.getName(), dpm.cgHandler.getLastObject()->statFile.getName());
}


/**
 * \brief MercuryCG is the postprocessing tool for extracting coarse-grained fields from particle data.
 * See \ref MercuryCG for details.
 */
int main(int argc, char *argv[]) {
    Mercury3D dpm;
    commandLineCG(dpm, argc, argv);
    return 0;
}