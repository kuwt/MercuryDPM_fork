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
#include <Species/LinearPlasticViscoelasticSlidingFrictionSpecies.h>
#include "Mercury3D.h"
#include "Boundaries/AngledPeriodicBoundary.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"

class ShearCell3D : public Mercury3D {
public:
    ShearCell3D (Mdouble splitRadius, Mdouble angularVelocity)
    {
        setName("ShearCell3DAbhiInitialConditions");
        readRestartFile();
        setRestarted(false);
        setName("ShearCell3DAbhi");
        std::cout << "openingAngle " << openingAngle << std::endl;
        setSaveCount(100000);
        species = dynamic_cast<LinearPlasticViscoelasticSlidingFrictionSpecies*>(speciesHandler.getObject(0));
        
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
        std::cout << "openingAngle " << openingAngle << std::endl;
    }

private:
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ///This is were the walls are implemented
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    void setupInitialConditions()
    {
    }

    void printTime() const
    {
        std::cout
            << "t=" << getTime()
            << " rotationAngleDeg " << getTime()*AngularVelocity*180.0/2.0/ constants::pi
            << " Ene " << getKineticEnergy()/getElasticEnergy()
            << std::endl;
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
public:
    LinearPlasticViscoelasticSlidingFrictionSpecies* species;
};


int main() {

    Mdouble splitRadius = 85/2.2; // Non-dimensional
    Mdouble angularVelocity = 0.001*2.0*constants::pi; // Slow rotation rate
    ShearCell3D SC(splitRadius, angularVelocity);

    SC.species->setSlidingFrictionCoefficient(0.5); //Change back to 0.01; for now simulation to 0.5
    SC.setXBallsAdditionalArguments("-v0 -solidf");
    SC.setHGridMaxLevels(2);
    SC.setTimeMax(1000.0);//200
    SC.dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    SC.fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    SC.restartFile.setFileType(FileType::ONE_FILE);
    SC.setSaveCount(10000);//1000
    SC.solve();
    return 0;
}
