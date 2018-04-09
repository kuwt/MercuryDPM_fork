#include "CSK.h"
#include<sstream>
#include<chrono>
#include<string>
#include<iostream>
/**
 * Simulates the mixing stage
 */

int main()
{
    //start measuring time
    std::chrono::time_point<std::chrono::system_clock> startClock, endClock;
    startClock = std::chrono::system_clock::now();

    Mdouble rpm = 36.0;
    unsigned nDoms = 2; //number of domains for parallel decomposition
    //std::string name = "longRateTest_Short-r0.035-rpm36-tc0.005";
    std::string name = "new_para_phase3_fill-r0.03-rpm6-nDomains2";
    //setup simulation
    CSK csk(name,false);
    csk.setRPM(rpm);
    csk.setTimeMax(600);
    csk.setSaveCount(5000);
    std::stringstream ss;
    ss << name << "-mixingStage" << "-rpm" << rpm;
    csk.setName(ss.str());
    //csk.write(std::cout,false);
    //csk.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    //Set the number of domains for parallel decomposition
    csk.setNumberOfDomains({nDoms,1,1});
    csk.solve();

    endClock = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = endClock - startClock;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << std::endl;
    logger(INFO, "Elapsed time for solving the PDE: % s", elapsed_seconds.count());


    return 0;
}
