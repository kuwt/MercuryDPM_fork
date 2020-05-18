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
    ShearCell3D (Mdouble splitRadius, Mdouble innerAngularVelocity, Mdouble outerAngularVelocity)
    {
        setName("ShearCell3DAbhiInitialConditions");
        readRestartFile();
        setRestarted(false);
        setName("ShearCell3DAbhi");
        std::cout << "openingAngle " << openingAngle << std::endl;
        setSaveCount(10000);
        species = dynamic_cast<LinearPlasticViscoelasticSlidingFrictionSpecies*>(speciesHandler.getObject(0));
        
        for (BaseParticle* p : particleHandler)
        {
			if (p->isFixed())
			{
				if (mathsFunc::square(p->getPosition().X)+mathsFunc::square(p->getPosition().Y)>mathsFunc::square(splitRadius))
				{
					OuterMPIndex.push_back(p->getIndex());
					Vec3D pos = p->getPosition();
					OuterMPRadius.push_back(sqrt(pos.X * pos.X + pos.Y * pos.Y));
					OuterMPAngle.push_back(atan(pos.Y / pos.X));
					OuterMPHeight.push_back(atan(pos.Z));
				} else {
					InnerMPIndex.push_back(p->getIndex());
					Vec3D pos = p->getPosition();
					InnerMPRadius.push_back(sqrt(pos.X * pos.X + pos.Y * pos.Y));
					InnerMPAngle.push_back(atan(pos.Y / pos.X));
					InnerMPHeight.push_back(atan(pos.Z));
				}
			}
		}
        InnerAngularVelocity = innerAngularVelocity;
        OuterAngularVelocity = outerAngularVelocity;
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
            << " rotationAngleDeg " << getTime()*(OuterAngularVelocity-InnerAngularVelocity)*180.0/2.0/ constants::pi
            << " Ene " << getKineticEnergy()/getElasticEnergy()
            << std::endl;
    }

    void actionsAfterTimeStep () {
        Mdouble dAlphaOuter = getTime()*OuterAngularVelocity;
        for (unsigned int i=0; i<OuterMPIndex.size(); i++) {
            Mdouble a = OuterMPAngle[i]+dAlphaOuter;
            if (a> openingAngle) OuterMPAngle[i] -= openingAngle;
            //if (a<=0.0) MPAngle[i] += openingAngle;
            Mdouble s = OuterMPRadius[i]*sin(a);
            Mdouble c = OuterMPRadius[i]*cos(a);
            BaseParticle* p = particleHandler.getObject(OuterMPIndex[i]);
            p->setPosition
                (Vec3D(c, s, OuterMPHeight[i] ));
            p->setVelocity
                (Vec3D(-OuterAngularVelocity*s, OuterAngularVelocity*c, 0.0));
            p->setAngularVelocity
                (Vec3D(0.0,0.0,OuterAngularVelocity));
        }
        Mdouble dAlphaInner = getTime()*InnerAngularVelocity;
        for (unsigned int i=0; i<InnerMPIndex.size(); i++) {
            Mdouble a = InnerMPAngle[i]+dAlphaInner;
            if (a> openingAngle) InnerMPAngle[i] -= openingAngle;
            //if (a<=0.0) MPAngle[i] += openingAngle;
            Mdouble s = InnerMPRadius[i]*sin(a);
            Mdouble c = InnerMPRadius[i]*cos(a);
            BaseParticle* p = particleHandler.getObject(InnerMPIndex[i]);
            p->setPosition
                (Vec3D(c, s, InnerMPHeight[i] ));
            p->setVelocity
                (Vec3D(-InnerAngularVelocity*s, InnerAngularVelocity*c, 0.0));
            p->setAngularVelocity
                (Vec3D(0.0,0.0,InnerAngularVelocity));
        }
    }

    Mdouble openingAngle;
    Mdouble InnerAngularVelocity;
    Mdouble OuterAngularVelocity;
    std::vector<unsigned int> OuterMPIndex;
    std::vector<Mdouble> OuterMPRadius;
    std::vector<Mdouble> OuterMPAngle;
    std::vector<Mdouble> OuterMPHeight;
    std::vector<unsigned int> InnerMPIndex;
    std::vector<Mdouble> InnerMPRadius;
    std::vector<Mdouble> InnerMPAngle;
    std::vector<Mdouble> InnerMPHeight;
public:
    LinearPlasticViscoelasticSlidingFrictionSpecies* species;
};


int main() {

    Mdouble splitRadius = 85e-3; // change back to 85
    Mdouble innerAngularVelocity = 0.01*2.0*constants::pi; // Slow rotation rate
    Mdouble outerAngularVelocity = 0.0*2.0*constants::pi; // Slow rotation rate
    ShearCell3D SC(splitRadius, innerAngularVelocity, outerAngularVelocity);

    SC.species->setSlidingFrictionCoefficient(0.5); //Change back to 0.01; for now simulation to 0.5
    SC.setXBallsAdditionalArguments("-v0 -solidf");
    SC.setHGridMaxLevels(2);
    SC.setTimeMax(200.0);
    SC.dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    SC.fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
    SC.restartFile.setFileType(FileType::ONE_FILE);// MULTIPLE_FILES_PADDED
    SC.setSaveCount(10000);
    SC.solve();
    return 0;
}
