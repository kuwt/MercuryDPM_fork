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
        //different particle positions every time the code is run
        //random.randomise();

        //define a new species (i.e. particle and contact properties); needs to be done before particles get inserted
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());

        //load Beger's data; the base plate is constructed from randomly placed particles, making it more frictional than a flat base plate
        logger.assert_always(readDataFile("Bottom.data"),"Input file could not be read");
        //exchange x and y position
        Vec3D pos;
        for (unsigned int i = 0; i < particleHandler.getNumberOfObjects(); i++) {
            pos = particleHandler.getObject(i)->getPosition();
            particleHandler.getObject(i)->setPosition(Vec3D(pos.Y, pos.X, 0));
        }
        setMax(Vec3D(getYMax(),getXMax(),getZMax()));
        // fix particles to the floor
        for (auto particle : particleHandler) {
            particle->fixParticle();
            particle->setSpecies(species);
        }
        // set fixed particle radius (not used, since we load the base plate from a data file)
        setFixedParticleRadius(particleHandler.getObject(0)->getRadius());

        //set name of output files
        setName("ContinuousTilt");
        // options for xballs output
        setXBallsAdditionalArguments("-v0 -solidf");

        // make size of inserted particles equal to the size of the particles on the base plate
        setInflowParticleRadius(particleHandler.getObject(0)->getRadius());
        inflowParticle_.setSpecies(species);

        // set particle density
        species->setDensity(2500);
        //set gravity (which is initially in downward direction)
        setGravity(Vec3D(0, 0, -9.8));
        //set stiffness, dissipation of the particle contacts such that the collision time is 5 ms and the restitution is 0.95/0.44 in normal/tangential direction (these values come from experimental data for glass particles, I think the Luger paper)
        Mdouble tc = 50e-4;
        Mdouble eps = 0.97; //0.97;
        Mdouble beta = 0.44; //0.44;
        Mdouble mass = species->getDensity() * mathsFunc::cubic(10e-3) * constants::pi / 6.;
        species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(tc, eps, beta, mass);
        // set an appropriate time step
        setTimeStep(tc / 50.);
        // set rolling/sliding friction coefficients and tangential stiffness/dissipation
        species->setSlidingFrictionCoefficient(0.1);
        species->setRollingStiffness(2. / 5. * species->getStiffness());
        species->setRollingDissipation(2. / 5. * species->getDissipation());
        species->setRollingFrictionCoefficient(species->getSlidingFrictionCoefficient() * 0.2);
        species->setSlidingFrictionCoefficient(species->getSlidingFrictionCoefficient() - species->getRollingFrictionCoefficient());
        logger(INFO,"Species properties: %",*species);
    }

    // this is outputted every time an output file is written
    void printTime() const override {
        std::cout
        << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
        << " a=" << std::setprecision(3) << std::left << std::setw(6) << getChuteAngleDegrees()
        << std::endl;
    }

    //this code is run before each time step of the simulation
    void actionsBeforeTimeStep() override {
        // after an initial settling time, the gravity vector is tilted at a continuous rate
        if (getTime() < 10) return;
        const Mdouble rate = 0.02; //in degree/second
        //initially increase the angle 10 times quicker (using a smooth function)
        const Mdouble da = 15;
        Mdouble dt = 0.02 * da / rate;
        Mdouble x = std::min(1.0, (getTime() - 10) / dt);
        Mdouble angle = da * x * x * (3 - 2 * x);
        setChuteAngle(rate * (getTime() - 10) + angle);
    }

    // determines what is written into the .ene output files
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

    // this defines the properties of the particles that are inserted
    void create_inflow_particle() {
        //particle radius
        inflowParticle_.setRadius(random.getRandomNumber(getMinInflowParticleRadius(), getMaxInflowParticleRadius()));
        //particle velocity
        inflowParticle_.setVelocity(Vec3D(0, 0, 0));
        //particle position
        Vec3D position;
        position.X = random.getRandomNumber(getXMin() + inflowParticle_.getRadius(),
                                            getXMax() - inflowParticle_.getRadius());
        position.Y = random.getRandomNumber(getYMin() + inflowParticle_.getRadius(),
                                            getYMax() - inflowParticle_.getRadius());
        position.Z = random.getRandomNumber(getZMin() + inflowParticle_.getRadius(),
                                            getZMax() - inflowParticle_.getRadius());
        inflowParticle_.setPosition(position);
    }

    // this code is run once at the beginning of the simulation
    void setupInitialConditions() override {
        // add flat walls in x, y, and negative z direction
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getLastObject());
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
        logger(INFO,"Added Nw=% walls",wallHandler.getNumberOfObjects());

        // use this code to make the simulation a periodic box, and remove the side walls
        //	WallsPeriodic.resize(0);
        //	Walls.resize(5);
        //	Walls[0].set(Vec3D( 0.0, 0.0,-1.0), -getZMin());
        //	Walls[1].set(Vec3D(-1.0, 0.0, 0.0), -getXMin());
        //	Walls[2].set(Vec3D( 1.0, 0.0, 0.0),  getXMax());
        //	Walls[3].set(Vec3D( 0.0,-1.0, 0.0), -getYMin());
        //	Walls[4].set(Vec3D( 0.0, 1.0, 0.0),  getYMax());

        //add particles above the base plate
        while (nCreated_ < numParticles) {
            create_inflow_particle();
            if (checkParticleForInteraction(inflowParticle_)) {
                particleHandler.copyAndAddObject(inflowParticle_);
                nCreated_++;
            } else setZMax(getZMax() + 0.00001);
        };
        logger(INFO,"Inserted N=% flowing and Nf=% fixed particles",particleHandler.getNumberOfObjects()-particleHandler.getNumberOfFixedObjects(),particleHandler.getNumberOfFixedObjects());
    }

    // how many particles to insert
    int numParticles = 10000;
    // track the number of inserted particles
    int nCreated_ = 0;
    // keep a default particle
    SphericalParticle inflowParticle_;
};

int main(int argc, char *argv[]) {
    LawinenBox md;
    // uncomment to add a unique number to the output files (so you can run the code repeatedly without overwriting the previous data)
    //md.autoNumber();
    // maximum duration of simulation
    md.setTimeMax(500); //actually finishes at 400
    // how often to write output
    md.setSaveCount(5000); //every half second
    md.eneFile.setSaveCount(100); //to get good plotting resolution
    // which output to write
    md.restartFile.setFileType(FileType::ONE_FILE);
    md.dataFile.setFileType(FileType::ONE_FILE);
    // this starts the simulation
    md.solve(argc, argv);
}
