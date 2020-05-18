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

class CSCInit : public Mercury3D {
public:
    CSCInit ()
    {
        eneTolerance = 1;
        sizeDistribution = cbrt(2.0);

        std::cout << "Reading file CSCWalls.restart" << std::endl;
        setName("CSCWalls");
        readRestartFile();
        setRestarted(false);
        setName("CSCInit");
        setXBallsAdditionalArguments("-v0 -solidf -3dturn 0");
        writeXBallsScript();
        species = speciesHandler.getObject(0);
        std::cout << "loaded " << particleHandler.getNumberOfObjects() 
            << " fixed particles" << std::endl;
    }

    bool continueSolve() const
    {
        static unsigned int counter = 0;
        if (++counter>100)
        {
            counter=0;
            if (getKineticEnergy()<eneTolerance*getElasticEnergy())
                return false;
        }
        return true;
    }

    void printTime() const
    {
        std::cout << "t=" << getTime() 
            << " Ene " << getKineticEnergy()/getElasticEnergy() << std::endl;
    }

    //add flow particles
    void setupInitialConditions() 
    {
        //hGridRebuild();
        //exit(-1);
        //number of particles to be inserted
        unsigned int n=(getXMax()-getXMin())
            *(getYMax()-getYMin())*(getZMax()-getZMin());
        std::cout << "Inserting " << n << " particles" << std::endl;
        //try to find new insertable particles
        unsigned int i = 0;
        SphericalParticle p;
        p.setSpecies(species);
        Mdouble s = sizeDistribution;
        Mdouble rMin=cbrt(0.5/(s*s + 1.0)/(s + 1.0));
        p.setRadius(s*rMin);
        std::cout << "Particle sizes from " << rMin << " to " << s*rMin 
            << " (sizeDistribution " << sizeDistribution << ")" << std::endl;
        Vec3D position;
        Mdouble a=0.0;
        while (i<n){
            position.X = random.getRandomNumber(getXMin()+p.getRadius(), 
                                                getXMax()-p.getRadius());
            position.Y = random.getRandomNumber(getYMin()+p.getRadius(), 
                                                getYMax()-p.getRadius());
            position.Z = random.getRandomNumber(getZMin()+p.getRadius(), 
                                                a*getZMax()-p.getRadius());
            p.setPosition(position);
            if (checkParticleForInteraction(p)) {
                //std::cout << i << std::endl;
                particleHandler.copyAndAddObject(p);
                p.setRadius(random.getRandomNumber(rMin, s*rMin));
                i++;
            } else { 
                a += 0.01/n;
            }
        }
        std::cout << "Inserted " << n << " particles" << std::endl;
    }

    void save()
    {
        std::cout << "Save " << particleHandler.getNumberOfObjects() 
            << " particles" << std::endl;

        //save data and restart file of first time step
        dataFile.open();
        outputXBallsData(dataFile.getFstream());
        dataFile.close();
        
        restartFile.open();
        writeRestartFile();
        restartFile.close();
    }

    ParticleSpecies* species;
    Mdouble eneTolerance;
    Mdouble sizeDistribution;
};

int main(int argc, char *argv[]) {
    CSCInit SC;
    SC.autoNumber();
    SC.random.randomise();
    SC.solve(argc, argv);
    SC.save();
    return 0;
}
