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
#include "Boundaries/AngledPeriodicBoundary.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Particles/LiquidFilmParticle.h"
#include "Species/LinearViscoelasticFrictionLiquidMigrationWilletSpecies.h"

class ShearCell3D : public Mercury3D {
public:
    ShearCell3D (Mdouble splitRadius, Mdouble angularVelocity, std::string fileName)
    {
        setName(fileName);
        readRestartFile();
        setRestarted(false);

        unsigned int i = 0;
        for (BaseParticle* p : particleHandler)
        {
            if (p->isFixed() && mathsFunc::square(p->getPosition().X)+mathsFunc::square(p->getPosition().Y)>mathsFunc::square(splitRadius))
            {
                MPIndex.push_back(p->getIndex());
                Vec3D pos = p->getPosition();
                MPRadius.push_back(sqrt(pos.X * pos.X + pos.Y * pos.Y));
                MPAngle.push_back(atan(pos.Y / pos.X));
                MPHeight.push_back(atan(pos.Z));
                i++;
            }
        }
        AngularVelocity = angularVelocity;
        AngledPeriodicBoundary* b = dynamic_cast<AngledPeriodicBoundary*>(boundaryHandler.getObject(0));
        openingAngle = b->getOpeningAngle();
        //write(std::cout,false);
        std::cout << "openingAngle " << openingAngle << std::endl;
        setFileType(FileType::NO_FILE);
        dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
        restartFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
        setSaveCount(1);
        setTimeMax(0.0015);
        //setTimeStep(getTimeStep()/3.0);
        //restartFile.getFstream().precision(15);
    }

private:
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ///This is were the walls are implemented
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    void setupInitialConditions() override
    {
    }

    void printTime() const override
    {
        Mdouble volTotP = 0.0;
        unsigned int nLB = 0;
        for (BaseParticle * const p : particleHandler)
        {
            LiquidFilmParticle* l = dynamic_cast<LiquidFilmParticle*>(p);
            //             if  (l->getLiquidVolume())
            //                 std::cout << " P" << l->getIndex() << " " << l->getLiquidVolume() <<std::endl;
            if (l == nullptr)
                std::cerr << "Error" << std::endl;
            volTotP += l->getLiquidVolume();
        }
        Mdouble volTotI = 0.0;
        for (BaseInteraction* i : interactionHandler)
        {
            LiquidMigrationWilletInteraction* l = dynamic_cast<LiquidMigrationWilletInteraction*>(i);
            if (l == nullptr)
                std::cerr << "Error" << std::endl;
            //            if  (l->getLiquidBridgeVolume())
            //                std::cout << " I" << l->getIndex() << " " << l->getLiquidBridgeVolume();
            volTotI += l->getLiquidBridgeVolume();
            if (l->getLiquidBridgeVolume() != 0.0)
                nLB++;
        }
        std::cout
        << "t=" << std::setprecision(6) << std::left << std::setw(6) << getNumberOfTimeSteps()
        << ", Vol="  << std::setprecision(8) << std::left << std::setw(6) << volTotP+volTotI
        << "(" << std::setprecision(6) << std::left << std::setw(6) << volTotP
        << " in particle)," << "(" << std::setprecision(6) << std::left << std::setw(6) << volTotI
        << " in bridge)," <<"#LB="  << std::setprecision(6) << std::left << std::setw(12) << nLB
        << std::endl;
        std::cout.flush();

    }


    void actionsAfterTimeStep () {
        Mdouble dAlpha = getTime()*AngularVelocity;
        for (unsigned int i=0; i<MPIndex.size(); i++) {
            Mdouble a = MPAngle[i]+dAlpha;
            if (a> openingAngle) MPAngle[i] -= openingAngle;
            //if (a<=0.0) MPAngle[i] += openingAngle;
            Mdouble s = MPRadius[i]*sin(a);
            Mdouble c = MPRadius[i]*cos(a);
            BaseParticle* p = particleHandler.getObject(MPIndex[i]);
            p->setPosition
                (Vec3D(c, s, MPHeight[i] ));
            p->setVelocity
                (Vec3D(-AngularVelocity*s, AngularVelocity*c, 0.0));
            p->setAngularVelocity
                (Vec3D(0.0,0.0,AngularVelocity));
        }
    }

    Mdouble openingAngle;
    Mdouble AngularVelocity;
    std::vector<unsigned int> MPIndex;
    std::vector<Mdouble> MPRadius;
    std::vector<Mdouble> MPAngle;
    std::vector<Mdouble> MPHeight;
};


int main(int argc, char *argv[]) {
    
    std::string fileName ("ShearCell3DInitialConditions");
    if (argc>1) 
    {
        fileName = argv[1];
        std::cout << "restarting from " << fileName << std::endl;
    }

    Mdouble splitRadius = 85e-3;
    Mdouble angularVelocity = 0.01*2.0*constants::pi;
    ShearCell3D SC(splitRadius, angularVelocity,fileName);
    SC.setXBallsAdditionalArguments("-v0 -solidf");
    SC.solve();
    return 0;
}
