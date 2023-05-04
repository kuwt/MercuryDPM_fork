//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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
#include <Logger.h>
#include <Math/Helpers.h>

int main(int argc, char* argv[])
{
    //read in params
    // material type
    std::string type = helpers::readFromCommandLine(argc,argv,"-speciesType",std::string("X"));
    // type of fit function from in- to output
    std::string fit = helpers::readFromCommandLine(argc,argv,"-fit",std::string("X"));
    // the string that is written into the output file
    std::string out;
    // the name of the output file
    std::string outFile = "CalibrationTest_" + type;

    // case specific translation from in- to output
    if (fit == "identity1") {
        // 1D input, 1D output: f(x)=x
        // read param0
        std::string string0 = helpers::readFromCommandLine(argc, argv, "-param0", std::string("0"));
        double param0 = std::stof(string0);
        // prepare output
        out = std::to_string(param0);
        outFile += "_" + string0 ;
    } else if (fit == "identity2") {
        // 2D input, 2D output: f(x)=x
        // read param0
        std::string string0 = helpers::readFromCommandLine(argc, argv, "-param0", std::string("0"));
        double param0 = std::stof(string0);
        std::string string1 = helpers::readFromCommandLine(argc, argv, "-param1", std::string("0"));
        double param1 = std::stof(string1);
        // prepare output
        out = std::to_string(param0) + ' ' + std::to_string(param1);
        outFile += "_" + string0 + "_" + string1;
    } else if (fit == "calibration43") {
        // 4D input, 3D output
        // read param0
        std::string string0 = helpers::readFromCommandLine(argc, argv, "-restitutionCoefficient", std::string("0"));
        double restitutionCoefficient = std::stof(string0);
        std::string string1 = helpers::readFromCommandLine(argc, argv, "-slidingFriction", std::string("0"));
        double slidingFriction = std::stof(string1);
        std::string string2 = helpers::readFromCommandLine(argc, argv, "-rollingFriction", std::string("0"));
        double rollingFriction = std::stof(string2);
        std::string string3 = helpers::readFromCommandLine(argc, argv, "-bondNumber", std::string("0"));
        double bondNumber = std::stof(string3);
        // fit
        //double angleOfRepose = slidingFriction;
        //double drum = slidingFriction + rollingFriction;
        //double ffc = 1.0/(slidingFriction+rollingFriction+bondNumber);
        double ffc = 1.0/bondNumber;
        double angleOfRepose = slidingFriction;
        double drum = rollingFriction;
        // prepare output
        out = std::to_string(ffc) + ' ' + std::to_string(angleOfRepose) + ' ' + std::to_string(drum);
        outFile += "_" + string0 + "_" + string1 + "_" + string2 + "_" + string3;
    } else if (fit == "calibration44") {
        // 4D input, 4D output
        // read param0
        std::string string0 = helpers::readFromCommandLine(argc, argv, "-restitutionCoefficient", std::string("0"));
        double restitutionCoefficient = std::stof(string0);
        std::string string1 = helpers::readFromCommandLine(argc, argv, "-slidingFriction", std::string("0"));
        double slidingFriction = std::stof(string1);
        std::string string2 = helpers::readFromCommandLine(argc, argv, "-rollingFriction", std::string("0"));
        double rollingFriction = std::stof(string2);
        std::string string3 = helpers::readFromCommandLine(argc, argv, "-bondNumber", std::string("0"));
        double bondNumber = std::stof(string3);
        // fit
        double dummy = restitutionCoefficient;
        double ffc = 1.0/bondNumber;
        double angleOfRepose = slidingFriction;
        double drum = rollingFriction;
        // prepare output
        out = std::to_string(dummy) + ' ' + std::to_string(ffc) + ' ' + std::to_string(angleOfRepose) + ' ' + std::to_string(drum);
        outFile += "_" + string0 + "_" + string1 + "_" + string2 + "_" + string3;
    } else {
        logger(ERROR,"fit function % unknown", fit);
    }

    //write output

    helpers::writeToFile(outFile + ".txt", out);
    logger(INFO,"Output to %: %",outFile, out);
    return 0;
}