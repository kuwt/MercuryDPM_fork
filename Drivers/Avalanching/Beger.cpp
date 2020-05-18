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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Chute.h"
#include "Walls/InfiniteWall.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"


//used to simulate Dirk Begers LawinenBox
class LawinenBox : public Chute {
public:

    LawinenBox() {
        readRestartFile("../ini/BegerThin.restart");
        setXBallsAdditionalArguments("-v0 -solidf -o 400");
        setAppend(false);
        setName("BegerMove");
    }

    ///rotates the chute counterclockwise about the origin (in degrees)
    void rotateChute(Mdouble Angle) {
        setChuteAngle(getChuteAngleDegrees() - Angle);
    };

    void writeEneTimeStep(std::ostream &os) const override {
        Mdouble m = particleHandler.getMass();
        Vec3D mom = particleHandler.getMomentum();
        static int width = (int) (os.precision() + 6);
        const InfiniteWall *iw = dynamic_cast <const InfiniteWall *> (wallHandler.getObject(2));
        Mdouble Angle = std::atan(iw->getNormal().Z / iw->getNormal().X) * 180. / constants::pi;
        os << std::setw(width) << getTime()
        << " " << std::setw(width) << -Vec3D::dot(getGravity(), mom)
        << " " << std::setw(width) << particleHandler.getKineticEnergy()
        << " " << std::setw(width) << particleHandler.getRotationalEnergy()
        << " " << std::setw(width) << getElasticEnergy()
        << " " << std::setw(width) << mom / m
        << " " << std::setw(width) << Angle
        << std::endl;
    }

    void actionsBeforeTimeStep() override {
        Mdouble RotationSpeed = 0.1; //in degree per second
        Mdouble CheckInterval = 0.5; //check every .. seconds

        static bool rotate = true;
        if (rotate) rotateChute(-RotationSpeed * getTimeStep()); //rotate 1 degree per second

        if (ceil(getTime() / CheckInterval) !=
            ceil((getTime() + getTimeStep()) / CheckInterval)) { //check every 0.1 seconds
            InfiniteWall *iw = dynamic_cast <InfiniteWall *> (wallHandler.getObject(2));
            Mdouble Angle = std::atan(iw->getNormal().Z / iw->getNormal().X) * 180. / constants::pi;
            if (Angle > 20) {
                RotationSpeed /= 10; //in degree per second
                Mdouble ene_kin = 0;
                for (unsigned int i = 0; i < particleHandler.getNumberOfObjects(); i++)
                    if (!particleHandler.getObject(i)->isFixed())
                        ene_kin += .5 * particleHandler.getObject(i)->getMass() *
                                   particleHandler.getObject(i)->getVelocity().getLengthSquared();
                rotate = ene_kin < 1e-5; //stop rotating if the kinetic energy increases
                std::cout << "Ene_kin=" << ene_kin
                << ", RotationSpeed=" << RotationSpeed
                << ", rotate=" << rotate
                << ", t=" << getTime()
                << ", theta=" << std::setprecision(6) << Angle
                << std::endl;
                if (Angle > 60) {
                    std::cout << "60 degree reached; exiting" << std::endl;
                    setTimeMax(getTime());
                } //end program
            }
        }
    }

    void create_inflow_particle() {
        Vec3D position;
        inflowParticle_.setRadius(random.getRandomNumber(getMinInflowParticleRadius(), getMaxInflowParticleRadius()));
        //inflowParticle_.computeMass();
        inflowParticle_.setVelocity(Vec3D(0, 0, 0));

        position.X = random.getRandomNumber(getXMin() + inflowParticle_.getRadius(),
                                            getXMax() - inflowParticle_.getRadius());
        position.Y = random.getRandomNumber(getYMin() + inflowParticle_.getRadius(),
                                            getYMax() - inflowParticle_.getRadius());
        position.Z = random.getRandomNumber(getZMin() + inflowParticle_.getRadius(),
                                            getZMax() - inflowParticle_.getRadius());
        inflowParticle_.setPosition(position);
    }

    void printTime() const override {
        Mdouble ene_kin = 0;
        for (unsigned int i = 0; i < particleHandler.getNumberOfObjects(); i++)
            if (!particleHandler.getObject(i)->isFixed())
                ene_kin += .5 * particleHandler.getObject(i)->getMass() *
                           particleHandler.getObject(i)->getVelocity().getLengthSquared();
        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
        << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
        << ", enekin=" << std::setprecision(3) << std::left << std::setw(6) << ene_kin
        << std::endl;
    }

    void setupInitialConditions() override {
    }

    SphericalParticle inflowParticle_;
};

int main(int argc, char *argv[]) {
    LawinenBox md;
    md.setTimeMax(1e20);
    md.setSaveCount((int) 1e5);
    md.restartFile.setFileType(FileType::MULTIPLE_FILES);
    md.eneFile.setSaveCount((int) 1e3);
    md.write(std::cout, false);
    md.solve(argc, argv);
}
