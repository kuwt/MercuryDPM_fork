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
#include "Walls/InfiniteWall.h"

class ClosedCSCRun : public Mercury3D {
public:
    ClosedCSCRun (Mdouble shearVelocity, Mdouble pressure)
    {
        std::cout << "Reading file ClosedCSCWall.restart" << std::endl;
        setName("ClosedCSCWalls");
        readRestartFile();
        setRestarted(false);
        setName("ClosedCSCRun");
        writeXBallsScript();
        setXBallsAdditionalArguments("-v0 -solidf -3dturn 0");
        setFileType(FileType::ONE_FILE);
        //restartFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
        std::cout << "loaded " << particleHandler.getNumberOfObjects() << 
            " fixed particles" << std::endl;

        //calculate required lidForce and lidMass
        lidForce = pressure*(getXMax()-getXMin())*(getYMax()-getYMin());
        lidMass = pressure*(getXMax()-getXMin())*(getYMax()-getYMin());
        lid = static_cast<InfiniteWall*>(wallHandler.getObject(wallHandler.getNumberOfObjects()-1));
        
        //setWall velocity
        for (BaseParticle* p: particleHandler)
        {
            if ( p->isFixed() ) 
            {
                if (p->getPosition().X>0.0)
                    p->setVelocity(Vec3D(0.0,0.5*shearVelocity,0.0));
                else
                    p->setVelocity(Vec3D(0.0,-0.5*shearVelocity,0.0));
            }
        }
        std::cout << "shear velocity " << shearVelocity << std::endl;
        setSaveCount(0.25/shearVelocity/getTimeStep());
        std::cout << "saving every " << dataFile.getSaveCount()*getTimeStep()
            << " time units" << std::endl;
    }

    void actionsAfterTimeStep()
    {
        Mdouble zMax = lid->getPosition().Z-1.0;
        // calculate the force applied to all lid particles
        Mdouble force = 0;
        for (BaseParticle* p: particleHandler)
            if ( p->isFixed() && p->getPosition().Z>zMax)
                force += p->getForce().Z;
        // calculate the change of lid velocity due to the lid force
        Mdouble dv = getTimeStep()*(force-lidForce)/lidMass;
        lid->setVelocity(lid->getVelocity()+Vec3D(0.0,0.0,dv));
        //apply the lid velocity to all lid particles
        for (BaseParticle* p: particleHandler)
            if ( p->isFixed() && p->getPosition().Z>zMax)
                p->setVelocity(p->getVelocity()+Vec3D(0.0,0.0,dv));            
    }

    
    void printTime() const
    {
        std::cout << "t=" << getTime() 
            << " Ene " << getKineticEnergy()/getElasticEnergy() 
            << " lid velocity " << lid->getVelocity().Z 
            << std::endl;
    }

    //add flow particles
    void setupInitialConditions() 
    {
    }
    
    Mdouble lidForce;
    Mdouble lidMass;
    InfiniteWall* lid;
};

int main(int argc, char *argv[]) {
    Mdouble shearVelocity = 1.0/40.0;
    Mdouble pressure = 30;
    ClosedCSCRun SC(shearVelocity,pressure);
    SC.setTimeMax(0);
    SC.solve(argc, argv);
    return 0;
}
