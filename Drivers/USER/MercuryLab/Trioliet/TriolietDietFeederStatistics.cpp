//Copyright (c) 2013-2014, The MercuryDPM Developers Team. All rights reserved.
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

#include "Statistics.h"
#include "Mercury3D.h"
#include <iostream>  // ostream
#include <iomanip>   // setw
#include <string.h>
#include <stdlib.h>
#include <sstream>

void computeTorque() {
    logger(INFO, "Computing torques");
    //Energy consumption (sum of torque on moving parts)
    Mercury3D dpm;
    dpm.setName("TriolietDietFeeder");
    dpm.readRestartFile();
    File &fStatFile = dpm.fStatFile;
    WallHandler &walls = dpm.wallHandler;
    fStatFile.setCounter(0);
    std::string dummy;

    struct Data {
        Mdouble time = 0;
        Vec3D torqueLeft = {0,0,0};
        Vec3D torqueRight = {0,0,0};
    };

    std::vector<Data> data;
    //read all fstat files
    while (fStatFile.open(std::fstream::in)) {
        getline(fStatFile.getFstream(), dummy);
        getline(fStatFile.getFstream(), dummy);
        getline(fStatFile.getFstream(), dummy);

        Mdouble time, delta, ctheta, fdotn, fdott;
        int index1 = 0, index2 = 0, counter = 0;
        Vec3D Contact, P1_P2_normal, P1_P2_tangential, force, torque;
        //go through each line in the fstat file; break when eof or '#' or newline
        Data step;
        while ((fStatFile.getFstream().peek() != -1) && (fStatFile.getFstream().peek() != '#')) { //(!fStatFile.eof())
            counter++;
            fStatFile.getFstream()
            >> time
            >> index1
            >> index2
            >> Contact
            >> delta
            >> ctheta
            >> fdotn
            >> fdott
            >> P1_P2_normal
            >> P1_P2_tangential;
            fStatFile.getFstream().ignore(256, '\n');
            force = fdotn * P1_P2_normal + 0*fdott * P1_P2_tangential;
            //Finished reading fstat file

            switch (index2) {
                case -1: //TriolietScrew Left
                case -3: //TriolietScrew Left
                case -5: //TriolietScrew Left
                case -7: //TriolietScrew Left
                    step.torqueLeft += Vec3D::cross(Contact - walls.getObject(4)->getPosition(), force);
                    break;
                case -2: //TriolietScrew Right
                case -4: //TriolietScrew Right
                case -6: //TriolietScrew Right
                case -8: //TriolietScrew Right
                    step.torqueRight += Vec3D::cross(Contact - walls.getObject(5)->getPosition(), force);
                    break;
            }
        }
        //logger(INFO, "finished timestep %", time);
        step.time = time;
        data.push_back(step);
        logger(INFO, "timeStep % (%): torque left %, right %", step.time, fStatFile.getCounter(), step.torqueLeft, step.torqueRight);
    }
    logger(INFO, "writing torque data to TriolietDietFeeder.torques");
    std::stringstream ss;
    ss << "time torqueLeft torqueRight\n";
    for (auto step : data) {
        ss << step.time << ' '
        << step.torqueLeft << ' '
        << step.torqueRight << '\n';
    }
    helpers::writeToFile("TriolietDietFeeder.torques",ss.str());
    //p 'TriolietDietFeeder.torques' u 1:7 w l title 'torque on left screw', '' u 1:4 w l title 'torque on right screw'
    fStatFile.close();
}

//common settings
void commonSettings(StatisticsVector<XY>& stats) {
	Mdouble h=0.2;
	Mdouble tMin = 0;
	auto s0 = stats.speciesHandler.getLastObject();
	auto s1 = stats.speciesHandler.copyAndAddObject(s0);
	stats.dataFile.setCounter(0);
	stats.readDataFile("TriolietDietFeeder.data");
	for (auto p : stats.particleHandler) {
		if (p->getPosition().Y > 0) {
			p->setSpecies(s1);
		}
	}
	stats.setCGTimeMin(tMin);
	stats.setCGWidth(0.5*h);
	stats.setTimeMaxStat(1e20);
	stats.setDoTimeAverage(false);
	stats.set_h(h);
}

int main(int argc, char *argv[]) {
    ///\todo Consider running the three statistics in parallel:
    if (argc <= 1) {
        logger(INFO,
               "\nRunning TriolietDietFeederStatistics without arguments only creates the torque data file.\n"
                "Postprocessing the TDF data takes time; therefore, the postprocessing steps are run in parallel.\n"
                "To do so, run the following 3 commands (for full and partial statistics):\n"
                "  ./TriolietDietFeederStatistics -o TriolietDietFeeder.timeAveraged.stat -h 0.05  -w 0.05 -tmin 5 -timeaverage 1 &\n"
                "  ./TriolietDietFeederStatistics -o TriolietDietFeeder.stat &\n"
                "  ./TriolietDietFeederStatistics -o TriolietDietFeeder.0.stat -indSpecies 0 &\n"
                "  ./TriolietDietFeederStatistics -o TriolietDietFeeder.1.stat -indSpecies 1 &\n"
                "  ./TriolietDietFeederStatistics -XSplit -o TriolietDietFeeder.2.stat -indSpecies 0 &\n"
                "  ./TriolietDietFeederStatistics -XSplit -o TriolietDietFeeder.3.stat -indSpecies 1 &\n");
        computeTorque();

        exit(-1);
    }

    logger(INFO, "Creating statistics:");
	StatisticsVector<XY> stats("TriolietDietFeeder");
    Mdouble h=0.2;
    Mdouble tMin = 0;
    stats.setCGTimeMin(tMin);
    stats.setCGWidth(0.5*h);
    stats.setTimeMaxStat(1e20);
    stats.setDoTimeAverage(false);
    stats.set_h(h);
    auto s0 = stats.speciesHandler.getLastObject();
    auto s1 = stats.speciesHandler.copyAndAddObject(s0);
    stats.dataFile.setCounter(0);
    stats.readDataFile("TriolietDietFeeder.data");

	if (!strcmp(argv[1],"-XSplit"))
    {
        logger(INFO,"Splitting data into two species at x=0");
        for (auto p : stats.particleHandler) {
            if (p->getPosition().X > 0) {
                p->setSpecies(s1);
            }
        }
        stats.readStatArguments(argc, argv);
    } else {
        logger(INFO,"Splitting data into two species at y=0");
        for (auto p : stats.particleHandler) {
            if (p->getPosition().Y > 0) {
                p->setSpecies(s1);
            }
        }
        stats.readStatArguments(argc + 1, argv - 1);
    }
	stats.statistics_from_fstat_and_data();

    logger(INFO,"To view output in MATLAB, copy loadstatistics.m and TriolietDietFeederStatistics.m into this directory and run");
}
