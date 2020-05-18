//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
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

#include<iostream>
#include "Mercury3DRestart.h"
using namespace std;

class Mercury3DRestarter : public Mercury3DRestart
{
//	void actionsBeforeTimeStep(){};
//	void setupInitialConditions(){setTime(time);}
//public:
//	double time;
    
//    virtual void actionsAfterTimeStep()
//    {
//        std::cout 
//            << "Time" << getTime() 
//            << " sc"<< restartFile.getSaveCount() 
//            << " nst"<< restartFile.getNextSavedTimeStep() 
//            << " t"<< getNumberOfTimeSteps() << std::endl;
//    }

    
    void writeOutputFiles() override
    {
        Mercury3D::writeOutputFiles();

        if (getWallTime() - getInitialWallTime() > getMaxWallTime())
        {
    		actionsAfterSolve();
    		finishStatistics();
    		closeFiles();
            std::cout << "Exiting for restarting after "
                << getWallTime() - getInitialWallTime() << "s" << std::endl;
            
            //set the restart command
            std::stringstream com("");
            restartFile.setCounter(restartFile.getCounter()-1);
            com << getClusterCommand() << " ./Mercury3DRestart -r " << restartFile.getFullName();
            
            //std::cout << com << std::endl;        
            std::cout << "system output:" << system(com.str().c_str()) << std::endl;        
	        exit(0);
    	}
    }
};

int main(int argc, char *argv[])
{
	if (argc<2) {
		cout << "Please enter the name of the simulation you want to restart and, optionally, the name of the simulation after restart" << endl;
		exit(-1);
	}
 	Mercury3DRestarter problem;
    problem.setSaveCount(2);
    problem.setMaxWallTime(10);
    problem.setClusterCommand("");
    problem.solve(argc, argv);
}
