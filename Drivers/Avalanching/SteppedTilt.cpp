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
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include "Chute.h"
#include "Walls/InfiniteWall.h"


//used to simulate Dirk Begers LawinenBox
class LawinenBox : public Chute {
public:

    LawinenBox() {
        auto S = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        random.randomise();
        //load Beger's bottom
        readDataFile("../ini/Bottom.data");
        //exchange x and y position
        Vec3D P;
        for (unsigned int i = 0; i < particleHandler.getNumberOfObjects(); i++) {
            P = particleHandler.getObject(i)->getPosition();
            particleHandler.getObject(i)->setPosition(Vec3D(P.Y, P.X, 0));
        }
        double xmax = getXMax();
        setXMax(getYMax());
        setYMax(xmax);

        setName("SteppedTilt");
        setXBallsAdditionalArguments("-v0 -solidf");
        setFixedParticleRadius(particleHandler.getObject(0)->getRadius());
        setInflowParticleRadius(particleHandler.getObject(0)->getRadius());
        S->setDensity(2500);
        setGravity(Vec3D(0, 0, -9.8));
        //collision time 5 ms
        Mdouble tc = 50e-4;
        Mdouble eps = 0.97; //0.97;
        Mdouble beta = 0.44; //0.44;
        Mdouble mass = S->getDensity() * mathsFunc::cubic(10e-3) * constants::pi / 6.;
        S->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, eps, beta, mass);
        setTimeStep(tc / 50.);
        S->setSlidingFrictionCoefficient(0.1);

        S->setRollingStiffness(2. / 5. * S->getStiffness());
        S->setRollingDissipation(2. / 5. * S->getDissipation());
        S->setRollingFrictionCoefficient(S->getSlidingFrictionCoefficient() * 0.2);
        S->setSlidingFrictionCoefficient(S->getSlidingFrictionCoefficient() - S->getRollingFrictionCoefficient());

        nCreated_ = 0;

        for (unsigned int i = 0; i < particleHandler.getNumberOfObjects(); i++)
            particleHandler.getObject(i)->fixParticle();
    }

    void actionsBeforeTimeStep() override {
        Mdouble DeltaAngle = 0.025; //in degree
        Mdouble CheckInterval = 0.1; //in seconds
        static Mdouble NextCheck = CheckInterval;
        static Mdouble NextIncrease = 0;
        static Mdouble ChuteAngle = 0;

        //move chute while you are in the increase interval
        if (getTime() < NextIncrease) {
            //http://www.comp-fu.com/2012/01/nukes-smooth-ramp-functions/
            // smooth: traditional smoothstep
            double x = (getTime() - (NextIncrease - CheckInterval)) / CheckInterval;
            setChuteAngle(ChuteAngle + x * x * (3 - 2 * x) * DeltaAngle);
            if (getChuteAngleDegrees() > 80) setTimeMax(getTime());
            //~ cout
            //~ << " t=" << x << get_t()
            //~ << " theta=" << getChuteAngleDegrees()
            //~ << endl;
        }

        //move chute while you are in the check interval
        if (getTime() > NextCheck) {
            NextCheck = CheckInterval + getTime();

            Mdouble ene_kin = 0;
            for (unsigned int i = 0; i < particleHandler.getNumberOfObjects(); i++)
                if (!particleHandler.getObject(i)->isFixed())
                    ene_kin += .5 * particleHandler.getObject(i)->getMass() *
                               particleHandler.getObject(i)->getVelocity().getLengthSquared();
            if (ene_kin < 1e-11) {
                NextIncrease = CheckInterval + getTime();
                ChuteAngle = getChuteAngleDegrees();
            }
            std::cout
            << " t=" << getTime()
            << " theta=" << getChuteAngleDegrees()
            << " Ene_kin=" << ene_kin
            << std::endl;
        }
    }

    void printTime() const override {
    }

    void writeEneTimeStep(std::ostream &os) const override {
        Mdouble m = particleHandler.getMass();
        Vec3D mom = particleHandler.getMomentum();
        static int width = (int) (os.precision() + 6);
        ///todo{Why is there a +6 here?  TW: to ensure the numbers fit into a constant width column}
        os << std::setw(width) << getTime()
        << " " << std::setw(width) << -Vec3D::dot(getGravity(), mom)
        << " " << std::setw(width) << particleHandler.getKineticEnergy()
        << " " << std::setw(width) << particleHandler.getRotationalEnergy()
        << " " << std::setw(width) << getElasticEnergy()
        << " " << std::setw(width) << mom / m
        << " " << std::setw(width) << getChuteAngleDegrees()
        << std::endl;
    }

    void create_inflow_particle() {
        inflowParticle_.setRadius(random.getRandomNumber(getMinInflowParticleRadius(), getMaxInflowParticleRadius()));
        //inflowParticle_.computeMass();
        inflowParticle_.setVelocity(Vec3D(0, 0, 0));

        Vec3D position;
        position.X = random.getRandomNumber(getXMin() + inflowParticle_.getRadius(),
                                            getXMax() - inflowParticle_.getRadius());
        position.Y = random.getRandomNumber(getYMin() + inflowParticle_.getRadius(),
                                            getYMax() - inflowParticle_.getRadius());
        position.Z = random.getRandomNumber(getZMin() + inflowParticle_.getRadius(),
                                            getZMax() - inflowParticle_.getRadius());
        inflowParticle_.setPosition(position);
    }

    void setupInitialConditions() override {
        InfiniteWall w0;
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0, 0.0, getZMin()));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(-1.0, 0.0, 0.0), Vec3D(getXMin(), 0.0, 0.0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(1.0, 0.0, 0.0), Vec3D(getXMax(), 0.0, 0.0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0.0, getYMin(), 0.0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 1.0, 0.0), Vec3D(0.0, getYMax(), 0.0));
        wallHandler.copyAndAddObject(w0);
        wallHandler.copyAndAddObject(w0);

        //		WallsPeriodic.resize(0);
        //		Walls.resize(5);
        //		Walls[0].set(Vec3D( 0.0, 0.0,-1.0), -getZMin());
        //		Walls[1].set(Vec3D(-1.0, 0.0, 0.0), -getXMin());
        //		Walls[2].set(Vec3D( 1.0, 0.0, 0.0),  getXMax());
        //		Walls[3].set(Vec3D( 0.0,-1.0, 0.0), -getYMin());
        //		Walls[4].set(Vec3D( 0.0, 1.0, 0.0),  getYMax());

        //number of flowing particles
        ///
        int M = 10000;
        particleHandler.setStorageCapacity(particleHandler.getNumberOfObjects() + M);
        hGridActionsBeforeTimeLoop();
        hGridActionsBeforeTimeStep();

        while (getNCreated() < M) {
            create_inflow_particle();
            if (checkParticleForInteraction(inflowParticle_)) {
                increaseNCreated();
            } else setZMax(getZMax() + 0.00001);
        };
    }

    int getNCreated() const {
        return nCreated_;
    }

    void increaseNCreated() {
        nCreated_++;
    }

    int nCreated_;
    SphericalParticle inflowParticle_;

};

int main(int argc, char *argv[]) {
    LawinenBox md;
    md.autoNumber();
    md.setTimeMax(1000); //actually finishes at 400
    md.setSaveCount(5000); //every half second
    md.eneFile.setSaveCount(100); //to get good plotting resolution
    md.restartFile.setFileType(FileType::ONE_FILE);
    md.dataFile.setFileType(FileType::ONE_FILE);
    md.solve(argc, argv);
    md.writeRestartFile();
}
