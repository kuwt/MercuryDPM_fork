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

//based on /storage2/usr/people/sluding/MDCC/C3DshearXL30/MU0_LONG2
#include "Mercury3D.h"
#include <sys/time.h>
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

class CSCRestart : public Mercury3D {
public:
    CSCRestart ()
    {
        initialTime = get_wall_time();
        maxWallTime = 20*3600; //kill after 20 hours
        setName("CSCRun");
        std::cout << "Reading file " << restartFile.getName() << std::endl;
        readRestartFile();
        restartFile.getFstream().precision(18);
        setAppend(true);
        printTime();
    }

    void writeOutputFiles()
    {
        if (get_wall_time()-initialTime<maxWallTime)
        {
            DPMBase::writeOutputFiles();
        }
        else
        {
            setTimeMax(getTime());
        }
    }

    void printTime() const
    {
        std::cout << "t=" << getTime() 
            << " Ene " << getKineticEnergy()/getElasticEnergy()
            << " wallTime " << get_wall_time()-initialTime << std::endl;

    }
    
    void restart() 
    {
        std::stringstream com("");
        com << "/home/weinhartt/bin/sclusterscriptexecute ./CSCRestart";
		std::cout << system(com.str().c_str()) << std::endl;	    
    }

    double maxWallTime;
    double initialTime;
};

int main(int argc, char *argv[]) {
    CSCRestart SC;
    SC.setTimeMax(2000);
    SC.solve(argc, argv);
    if (SC.getTimeMax()<2000)
        SC.restart();
    return 0;
}
