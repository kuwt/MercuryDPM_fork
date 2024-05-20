//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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

#include <filesystem>
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Mercury3D.h"
#include "Math/ExtendedMath.h"

/*!
 * \class RotatingDrumBidisperse
 * \brief This driver file is used to set up a rotating drum simulation with a bidisperse particle distribution.
 *
 * \details The simulation is divided into setting up the initial particle packing and the actual rotation of the drum.
 *
 * The Froude number determines the rotation rate.
 * The relative drum radius and width are based on the mean particle diameter (d50) and determine how many particles of
 * the d50 will fit into the drum in a row.
 * The fillFraction determines the fill level of the drum in percentage.
 *
 * The RotatingDrumBidisperseInitialise is setting up the initial condition and is divided into four stages:
 * 1. Insertion of large particles until settled
 * 2. Insertion of small particles until settled
 * 3. Scaling up stiffness by 10% each timestep
 * 4. Settling of particles
 *
 * The RotatingDrumBidisperse is the actual simulation of the rotating drum. It is divided into two stages:
 * 1. Rotation of the drum for 2 rotations
 * 2. Writing the angle of repose to to a .txt file
 *
 * \param[in] diameterRatio     ratio of the large to small particle diameter
 * \param[out] angle_drum       angle of repose of the particle packing
 *
 */
class RotatingDrumBidisperseInitialise : public Mercury3D
{

public:
    RotatingDrumBidisperseInitialise(Mdouble diameterRatio)
    {
        // parameters
        const Mdouble dLarge = 3e-3; // particle diamter of the large particles [m]
        const Mdouble d50 = (dLarge+dLarge*diameterRatio)/2.0; // mean particle diameter (d50) [m]
        const Mdouble relDrumRadius = 10; // relative drum radius based on particle diameters [m]
        const Mdouble relDrumWidth = 5; // relative drum width based on particle diameters [m]
        const Mdouble fillFraction = 0.5; // fill level of the drum [-]

        const Mdouble stiffness = 100; // stiffness of particle contacts (N/m
        const Mdouble dissipation = 10e-3; // dissipation of particle contacts (Ns/m)
        const Mdouble density = 1500; // particle density (kg/m^3)
        const Mdouble slidingFriction = 0.5; // sliding friction coefficient (-)
        const Mdouble rollingFriction = 0.5; // rolling friction coefficient (-)

        // set name of output files
        std::string simulationName = "RotatingDrumBidisperseInitialise" + std::to_string(diameterRatio);
        setName(simulationName);

        // remove old output files
        removeOldFiles();

        // turn on paraview output
//        setParticlesWriteVTK(true);
//        wallHandler.setWriteVTK(true);

        //set drum radius, depth
        drumRadius_ = relDrumRadius*d50;
        logger(INFO, "drumRadius = %*d50", relDrumRadius);
        drumDepth_ = relDrumWidth*d50;
        logger(INFO, "drumDepth = %*d50", relDrumWidth);
        // domain
        setMax({drumRadius_,0.5*drumDepth_,drumRadius_});
        setMin({-drumRadius_,-0.5*drumDepth_,-drumRadius_});

        //set gravity
        setGravity({0,0,-9.8});
        //rotate twice
        setTimeMax(1e20);

        // set species and contact properties
        LinearViscoelasticFrictionSpecies speciesLarge;
        speciesLarge.setDensity(density);
        speciesLarge.setStiffness(stiffness);
        speciesLarge.setDissipation(dissipation);
        speciesLarge.setSlidingFrictionCoefficient(slidingFriction);
        speciesLarge.setSlidingStiffness(2./7.*stiffness);
        speciesLarge.setSlidingDissipation(2./7.*dissipation);
        speciesLarge.setRollingFrictionCoefficient(rollingFriction);
        speciesLarge.setRollingStiffness(2./7.*stiffness);
        speciesLarge.setRollingDissipation(2./7.*dissipation);
        sLarge_ = speciesHandler.copyAndAddObject(speciesLarge);

        LinearViscoelasticFrictionSpecies speciesSmall;
        speciesSmall.setDensity(density);
        speciesSmall.setStiffness(stiffness);
        speciesSmall.setDissipation(dissipation);
        speciesSmall.setSlidingFrictionCoefficient(slidingFriction);
        speciesSmall.setSlidingStiffness(2./7.*stiffness);
        speciesSmall.setSlidingDissipation(2./7.*dissipation);
        speciesSmall.setRollingFrictionCoefficient(rollingFriction);
        speciesSmall.setRollingStiffness(2./7.*stiffness);
        speciesSmall.setRollingDissipation(2./7.*dissipation);
        sSmall_ = speciesHandler.copyAndAddObject(speciesSmall);

        // get collisionTime and restitutionCoefficient
        double collisionTime = sSmall_->getCollisionTime(sSmall_->getMassFromRadius(dLarge*diameterRatio/2));
        double restitution = sSmall_->getRestitutionCoefficient(sSmall_->getMassFromRadius(dLarge*diameterRatio/2));
        logger(INFO,"collisionTime %, restitution coefficient %", collisionTime, restitution);

        setTimeStep(0.05 * collisionTime);
        setSaveCount(1000);

        // set wall species
        // define side-wall species (no friction/cohesion)
        auto frictionlessWallSpecies = speciesHandler.copyAndAddObject(speciesLarge);
        auto mixedSpecies = speciesHandler.getMixedObject(sLarge_, frictionlessWallSpecies);
        mixedSpecies->setRollingFrictionCoefficient(0.0);
        mixedSpecies->setSlidingFrictionCoefficient(0.0);
        // define drum-wall species (high friction)
        auto frictionalWallSpecies = speciesHandler.copyAndAddObject(speciesLarge);
        mixedSpecies = speciesHandler.getMixedObject(sLarge_, frictionalWallSpecies);
        mixedSpecies->setRollingFrictionCoefficient(std::max(1.0,
                                                             sLarge_->getSlidingFrictionCoefficient()));
        mixedSpecies->setSlidingFrictionCoefficient(std::max(1.0, sLarge_->getRollingFrictionCoefficient()));

        // adding cylindrical wall
        AxisymmetricIntersectionOfWalls cylinder;
        cylinder.setSpecies(frictionalWallSpecies);
        cylinder.setAxis({0.0,1.0,0.0});
        cylinder.addObject({1.0,0.0,0.0},{drumRadius_,0.0,0.0});
        wallHandler.copyAndAddObject(cylinder);

        //adding side walls
        InfiniteWall side;
        side.setSpecies(frictionlessWallSpecies);
        side.set(Vec3D(0., -1., 0.), Vec3D(0.0, -0.5 * drumDepth_, 0.0));
        wallHandler.copyAndAddObject(side);
        side.set(Vec3D(0., 1., 0.), Vec3D(0.0, 0.5 * drumDepth_, 0.0));
        wallHandler.copyAndAddObject(side);
        logger(INFO,"Walls have been placed");

        //add type of particles to insert
        SphericalParticle pLarge;
        pLarge.setSpecies(sLarge_);
        pLarge.setRadius(dLarge/2);
        SphericalParticle pSmall;
        pSmall.setSpecies(sSmall_);
        pSmall.setRadius(dLarge*diameterRatio/2);

        // insert particles
        CubeInsertionBoundary insertionLarge;
        insertionLarge.set(&pLarge, 1000, getMin(), getMax(), Vec3D(0, 0, 0), Vec3D(0, 0, 0));
        insertionLarge.setInitialVolume(0.5*fillFraction*constants::pi/2*pow(drumRadius_,2)*drumDepth_);
        CubeInsertionBoundary* insertionLargeP = boundaryHandler.copyAndAddObject(insertionLarge);
        CubeInsertionBoundary insertionSmall;
        insertionSmall.set(&pSmall, 1000, {getXMin(),getYMin(),0}, getMax(), Vec3D(0, 0, 0), Vec3D(0, 0, 0));
        insertionSmall.setInitialVolume(0);
        CubeInsertionBoundary* insertionSmallP = boundaryHandler.copyAndAddObject(insertionSmall);
        insertionSmall_ = insertionSmallP;
        insertionLarge_ = insertionLargeP;
        logger(INFO, "Starting stage %: insertion of large particles with an initial volume of %\n", stage_, insertionLargeP->getInitialVolume());

        // t_c/t_g < 200; t_g = \sqrt(d/g); A rotating drum should be driven by gravitational forces and ignore microscale
        // effects due to particle-particle collisions. Therefore, we try to keep the gravitational timescale higher than the
        // collision timescale.
        if (sqrt(dLarge/getGravity().Z) > 200)
        {
            logger(WARN, "Beware, collision timescale is much larger than gravitational timescale.\nReduce the collision time.");
        }
    }

    void actionsAfterTimeStep() override
    {
        switch(stage_)
        {
            // first stage is insertion of large particles until settled
            case 1:
                if (getKineticEnergy()/getElasticEnergy() < 1)
                {
                    // initiate second stage: insetion of small particles
                    if (insertionLarge_->getInitialVolume() <= insertionLarge_->getVolumeOfParticlesInserted())
                    {
                        insertionSmall_->setInitialVolume(insertionLarge_->getInitialVolume());
                        stage_ = 2;
                        logger(INFO, "Initiating stage %: Insertion and settling of small particles\n", stage_);
                        break;
                    }
                }
            // second stage is insertion of small particles until settled
            case 2:
                if (getKineticEnergy()/getElasticEnergy() < 1e-3)
                {
                    // initiate third stage: scaling up stiffness
                    stage_ = 3;
                    logger(INFO, "Initiating stage %: scaling up stiffness by 10% each timestep\n", stage_);
                    // activate backgroundDrag to reduce kinetic energy from overlaps while scaling up stiffness
                    setBackgroundDrag(insertionLarge_->getMassOfParticlesInserted());
                    logger(INFO, "Activating BackgroundDrag = %\n", getBackgroundDrag());
                    break;
                }
            case 3:
                if (getKineticEnergy()/getElasticEnergy() < 1e-3)
                {

                    // scale up stiffness by 10 percent of the actual stiffness each timestep
                    sSmall_->setStiffness(sSmall_->getStiffness()*1.1);
                    sLarge_->setStiffness(sLarge_->getStiffness()*1.1);
                    // get collisionTime and restitutionCoefficient
                    double collisionTime = sSmall_->getCollisionTime(sSmall_->getMassFromRadius(particleHandler.getSmallestParticle()->getRadius()));
                    // set new timestep
                    setTimeStep(0.05 * collisionTime);
//                    logger(INFO,"new timestep = %", getTimeStep());
                    // if final stiffness is reached
                    if (sLarge_->getStiffness() >= 1e5 && sSmall_->getStiffness() >= 1e5)
                    {
                        // initiate fourth stage: settling of particles
                        stage_ = 4;
                        setBackgroundDrag(0);
                        logger(INFO, "Initiating stage %: settling of particles\nDeactivating BackgroundDrag", stage_);
                        break;
                    }
                }
        }
    }

    bool continueSolve() const override
    {
        static unsigned int counter = 0;
        if (++counter>100)
        {
            counter=0;
            if (stage_ == 4 && getKineticEnergy()/getElasticEnergy() < 1e-3)
            {
                return false;
            }
        }
        return true;
    }

    void printTime() const override
    {
        logger(INFO, "t=% NLarge=% NSmall=% Ene=%", getTime(),insertionLarge_->getNumberOfParticlesInserted(),insertionSmall_->getNumberOfParticlesInserted() , getKineticEnergy() / getElasticEnergy());
    }

    void save()
    {
        logger(INFO, "Save % particles", particleHandler.getNumberOfObjects());

        //save data and restart file of first time step
        dataFile.open();
        outputXBallsData(dataFile.getFstream());
        dataFile.close();

        restartFile.open();
        writeRestartFile();
        restartFile.close();
    }


private:

    int stage_ = 1;
    Species<LinearViscoelasticNormalSpecies,FrictionSpecies>* sLarge_;
    Species<LinearViscoelasticNormalSpecies,FrictionSpecies>* sSmall_;
    CubeInsertionBoundary* insertionSmall_;
    CubeInsertionBoundary* insertionLarge_;
    Mdouble drumRadius_;
    Mdouble drumDepth_;
};

/*!
 *
 */

class RotatingDrumBidisperse : public Mercury3D
{
public:
    RotatingDrumBidisperse(int argc, char* argv[])
    {
        // add -diameterRatio xx to the execution command to change the diemater ratio of large to small particles
        diameterRatio_ = helpers::readFromCommandLine(argc, argv, "-diameterRatio", 1.0/3.0); // ratio to smaller diameter

        std::string restartFile = "RotatingDrumBidisperseInitialise" + std::to_string(diameterRatio_) + ".restart";
        setName("RotatingDrumBidisperseInitialise" + std::to_string(diameterRatio_));
        setXBallsAdditionalArguments("-v0 -solidf");
        writeXBallsScript();
        if (std::filesystem::exists(restartFile))
        {
            logger(INFO, "Reading file RotatingDrumBidisperseInitialise%.restart\n", diameterRatio_, Flusher::NO_FLUSH);
        }
        else
        {
            RotatingDrumBidisperseInitialise drumInit(diameterRatio_);
            drumInit.random.randomise();
            drumInit.solve();
            logger(INFO, "Initial conditions have been set.\nPlease check if everything went well and run me again to start the simulation\n");
            exit(0);
        }
        readRestartFile();
        setTime(0);
        setRestarted(false);
        setName("RotatingDrumBidisperse" + std::to_string(diameterRatio_));
        // remove old output files
        removeOldFiles();
        particleSpecies_ = speciesHandler.getObject(0);
        cylinderSpecies_ = wallHandler.getObject(0);
        logger(INFO, "loaded % fixed particles", particleHandler.getNumberOfObjects());

        // parameters
        const Mdouble froude = 0.01; // froude number
        const Mdouble revolutions = 2.0; // total number of revolutions

        // turn on paraview output
        setParticlesWriteVTK(true);
        wallHandler.setWriteVTK(true);

        //set gravity and rotation rate
        rotRate_ = froude*sqrt(getGravity().getLength()/getXMax());
        //rotate twice
        setTimeMax(2.0*constants::pi/rotRate_*revolutions);
        logger(INFO,"Simulating % revolutions at Froude number %",revolutions, froude);

        // set save count of files
        setSaveCount(getTimeMax()/getTimeStep()/1000);

        cylinderSpecies_->setAngularVelocity(Vec3D(0.0,rotRate_,0.0));
    }

    void printTime() const override
    {
        logger(INFO, "t % N % AoR % del %", getTime()/getTimeMax(), particleHandler.getNumberOfObjects(),
               atan(getCentreOfMass().X/getCentreOfMass().Z)*180.0/constants::pi,
               interactionHandler.getMeanOverlap()/particleHandler.getMeanRadius());
    }

    void writeResults() {
        // center of mass; both x and z-values  should be negative
        Vec3D com = getCentreOfMass();
        angle_drum_ = atan(com.X/com.Z)*180/constants::pi;
        logger(INFO,"Angle % COM (%, %)",angle_drum_,com.X,com.Z);
        double r = sqrt(com.X * com.X + com.Z * com.Z)/getXMax();
        if (r<0.2) {
            logger(WARN,"Center of mass is too close to the center, |COM|/R=%, setting angle to 90",r);
            angle_drum_=90;
        }
        if (com.X>0 || com.Z>0) {
            logger(WARN,"Center of mass (% %) is not in fourth quadrant, setting angle to 90",com.X,com.Z);
            angle_drum_=90;
        }
        helpers::writeToFile(getName()+".txt", helpers::toString(angle_drum_));
    }

private:

    ParticleSpecies* particleSpecies_;
    BaseWall *cylinderSpecies_;
    Mdouble diameterRatio_;
    Mdouble rotRate_;// this is in rad/s, so no need to convert again
    Mdouble angle_drum_;
};

//! Defines a rotating drum simulation
int main(int argc, char* argv[])
{
    RotatingDrumBidisperse drum(argc, argv);
    drum.solve();
    drum.writeResults();
    return 0;
}