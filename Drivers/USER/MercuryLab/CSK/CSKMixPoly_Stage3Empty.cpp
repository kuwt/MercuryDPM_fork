#include "CSK.h"
#include<sstream>
#include<chrono>
#include<string>
#include<iostream>
/**
 * Simulates the emptying stage
 */

int main()
{

    //start measuring time
    std::chrono::time_point<std::chrono::system_clock> startClock, endClock;
    startClock = std::chrono::system_clock::now();
    Mdouble timeMax = 600.0;
    Mdouble rpm = 36.0; //Set at 6 rpm to correspond to 'real' system
    std::string name = "phaseThreeInitialTest-r0.04-rpm6-mixingStage-rpm36"; 
    //setting the outflow rate (in m^3 per second - note initial volume should be 0.6)
    double outflowRate = 0.6 / timeMax;    
    //allowing the collision time to be reset if desirable
    //Mdouble collisionTime = 0.005 / 2.5; //Currently at Thomas-suggested value :-)
    //setup simulation for outflow
    CSK csk(name, outflowRate);
    csk.setRPM(rpm);
    csk.setTimeMax(timeMax);
    //csk.setTimeStep(0.02 * collisionTime);
    csk.setSaveCount(5000);
    std::stringstream ss;
    ss << name << "-emptyingStage" << "-rpm" << rpm << "-tMax" << timeMax;
    csk.setName(ss.str());
    //csk.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    csk.solve();

    //end measuring time
    endClock = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = endClock - startClock;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << std::endl;
    logger(INFO, "Elapsed time for solving the PDE: % s", elapsed_seconds.count());
    return 0;
}
