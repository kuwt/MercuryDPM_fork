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
#include "Mercury3DRestart.h"

class CSCRun : public Mercury3DRestart
{
public:

    CSCRun(Mdouble shearVelocity): shearVelocity_(shearVelocity){};
    
    void setupInitialConditions()
    {
        Mdouble timeMax = getTimeMax();
        setName("CSCInit");
        std::cout << "Reading file " << restartFile.getName() << std::endl;
        readRestartFile();
        //setRunNumber(0); //restart doesn't work with autonumbered init files
        setTimeMax(timeMax);
        setTime(0);
        setRestarted(false);
        setName("CSCRun");
        writeXBallsScript();
        setXBallsAdditionalArguments("-v0 -solidf -3dturn 1");
        setFileType(FileType::ONE_FILE);
        //restartFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
        std::cout << "loaded " << particleHandler.getNumberOfObjects() <<
            " fixed particles" << std::endl;

        //setWall velocity
        for (BaseParticle* p : particleHandler) {
            if (p->isFixed()) {
                if (p->getPosition().X > 0.0)
                    p->setVelocity(Vec3D(0.0, 0.5 * shearVelocity_, 0.0));
                else
                    p->setVelocity(Vec3D(0.0, -0.5 * shearVelocity_, 0.0));
            }
        }
        //std::cout << "Shear velocity " << shearVelocity_ << std::endl;
        
        //set save count such that wall particles move 1/4 part. diam. per saved time step
        setSaveCount(0.25 / shearVelocity_ / getTimeStep());
        //std::cout << "Saving every " << dataFile.getSaveCount() * getTimeStep()
        //    << " time units" << std::endl;
    }

    void printTime() const
    {
        std::cout << "t=" << getTime()
            << " ene " << getKineticEnergy() / getElasticEnergy()
            << " wallTime " << getWallTime() - getInitialWallTime() << std::endl;
    }
    
    Mdouble shearVelocity_;
};

int main(int argc, char *argv[])
{
    Mdouble shearVelocity = 1.0/40.0; //divide by 20 to get max inertial number
    CSCRun SC(shearVelocity);
    if (false) {
        //set to true to test the restarting
        SC.setSaveCount(2);
        SC.setMaxWallTime(10);
        SC.setClusterCommand("");
        SC.setTimeMax(1);
    } else {
        SC.setSaveCount(100);
        SC.setMaxWallTime(20*3600);//*3600);//kill after 20 hours
        SC.setClusterCommand("~/bin/sclusterscriptexecute");
        SC.setTimeMax(4000.0);
    }
    SC.solve(argc, argv);
    return 0;
}
