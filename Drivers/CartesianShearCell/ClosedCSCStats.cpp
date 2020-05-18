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
#include "StatisticsVector.h"

class ClosedCSCStats : public StatisticsVector<XZ> {
public:
    ClosedCSCStats ()
    {
        setName("ClosedCSCRun");
        //restartFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
        std::cout << "Reading file " << restartFile.getFullName() << std::endl;
        readRestartFile();
        setRestarted(false);
        setName("ClosedCSCStats");
        setFileType(FileType::NO_FILE);
        restartFile.setFileType(FileType::ONE_FILE);
        statFile.setFileType(FileType::ONE_FILE);
        
        setNTimeAverageReset(100);
        setDoPeriodicWalls(false);
		setN(50);
		//setCGShape("Lucy");
        setCGWidth(1.0);
		setSaveCount(50);
		setCGTimeMin(getTime());
		setTimeMax(20); //40'000 time steps
    }

    void printTime() const
    {
        std::cout << "t=" << getTime() 
            << " Ene " << getKineticEnergy()/getElasticEnergy() << std::endl;
    }

    //add flow particles
    void setupInitialConditions() 
    {
        std::string filename;
        filename = getName() 
            + "SC" + std::to_string(restartFile.getSaveCount())
            + "N" + std::to_string(getNX())
            + "T" + std::to_string(static_cast<int>(round(getTimeMax())));
        setName(filename);
        setTime(getCGTimeMin());
        setTimeMax(getCGTimeMin()+getTimeMax()); //40'000 time steps
    }
};

int main(int argc, char *argv[]) {
    ClosedCSCStats SC;
    argv--; argc++;
    SC.readStatArguments(argc, argv);
    SC.solve();
    return 0;
}

// ./fstatistics CSCRun -timeaverage 0 -o CSCRun.O.stat 
// ./fstatistics CSCRun -tmin 50 -n 50 -stattype XZ
// ./fstatistics CSCRun -w 2 -tmin 70 -n 50 -stattype XZ -o CSCRun.W2.stat
