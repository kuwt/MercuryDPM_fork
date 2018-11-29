//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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
//#include <Species/Species.h>
//#include <Species/LinearViscoelasticSpecies.h>
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Mercury3D.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Boundaries/StressStrainControlBoundary.h>
#include <Boundaries/CubeInsertionBoundary.h>
#include "Boundaries/LeesEdwardsBoundary.h"
#include "addSpecies.h"
#include <iomanip>
#include <sstream>
#include <utility>
#include <cmath>
#include <Walls/InfiniteWall.h>

/**
 * The simulation is in one of two stages:
 * - Compression: compress the samples until a fixed normal stress is reached
 * - Shear: shears the sample under fixed normal stress
 */
enum class Stage {
    Compression,
    Shear
};

/**
 * Simulation of a shear box, with particles between a fixed lower wall and an upper wall in z-direction, and periodic walls in x- and y-direction. The upper wall can move up and down to control the normal stress, and moves sideways during the shear stage.
 */
class ShearCellWall : public Mercury3D {
public:

    /**
     * Constructor.
     * @param type determines particle species
     * @param compressiveStress vector of goal pressures that have to be reached in successive stages
     * @param shearVelocity velocity of the upper wall during the shear stages
     * @param relativeCohesionStiffness determines cohesive behaviour of wall-particle interactions
     * @param slidingFriction sliding friction coefficient of wall-particle interactions
     * @param rollingFriction rolling friction coefficient of wall-particle interactions
     * @param shearRate the shear strain rate which determines the actual shearVelocity at sample height after compression
     */
    ShearCellWall(SpeciesType type, std::vector<double> compressiveStress, Vec3D shearVelocity,
                  double relativeCohesionStiffness, double slidingFriction, double rollingFriction, double shearRate)
            : compressiveStress_(std::move(compressiveStress)), shearVelocity_(shearVelocity), shearRate_(shearRate) {
        //some standard settings
        setXBallsAdditionalArguments("-3dturn 1 -v0 -solidf");
        setGravity({0, 0, 0});
        setTimeMax(1e20);

        //add species
        const auto &psd = getPSD(type);
        logger.assert(hasPSD(type), "Material type % not allowed", type);
        auto particleSpecies = addSingleSpecies(type, getMedianParticleRadius(psd), getMinParticleRadius(psd), *this,
                                                true);
        auto wallSpecies = addSingleSpecies(SpeciesType::SteelWall, getMedianParticleRadius(psd),
                                            getMinParticleRadius(psd), *this, true);

        //set wall-particle properties
        auto wallParticleSpecies = dynamic_cast<LinearPlasticViscoelasticFrictionMixedSpecies *>(speciesHandler.getMixedObject(
                particleSpecies, wallSpecies));
        logger.assert_always(wallParticleSpecies, "wallParticleSpecies is null pointer");
        double cohesionStiffness = relativeCohesionStiffness * wallParticleSpecies->getUnloadingStiffnessMax();
        wallParticleSpecies->setCohesionStiffness(cohesionStiffness);
        wallParticleSpecies->setSlidingFrictionCoefficient(slidingFriction);
        wallParticleSpecies->setRollingFrictionCoefficient(rollingFriction);
        logger(INFO, "Particle and wallParticle species\n%\n%", *particleSpecies, *wallParticleSpecies);

        // Create a cubic volume of 8.5 mean particle diameters side length.
        // It's not isotropically compressed, so I need to strech the z-dimension initially to get a cubic volume.
        setMax(17.0*getMedianParticleRadius(psd)*Vec3D(1,1,3));
        setMin(-getMax());

        // add periodic boundaries
        PeriodicBoundary boundary;
        boundary.set({1, 0, 0}, getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(boundary);
        boundary.set({0, 1, 0}, getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(boundary);

        // add walls
        InfiniteWall wall;
        wall.setSpecies(wallSpecies);
        wall.set({0, 0, -1}, getMin());
        wallHandler.copyAndAddObject(wall);
        wall.set({0, 0, 1}, getMax());
        //wall.setVelocity({0,0,-1E-3});
        shearedWall_ = wallHandler.copyAndAddObject(wall);

        //Here we setup the displacement control
        logger.assert_always(!compressiveStress_.empty(), "Need at least one compressive stress");
        auto wallArea = (getXMax() - getXMin()) * (getYMax() - getYMin());
        auto forceGoal = Vec3D(0, 0, compressiveStress_[nStage_] * wallArea);
        shearedWall_->setVelocityControl(forceGoal, gainFactor);
        logger(INFO, "Start with compression #% (p=%)", nStage_, compressiveStress_[nStage_]);

        //define initial volume fraction
        const Mdouble initialVolumeFraction = 0.25;

        //add particles
        CubeInsertionBoundary insertion;
        BaseParticle particle(speciesHandler.getObject(0));
        insertion.set(&particle, 100000, getMin(), getMax(), Vec3D(0, 0, 0), Vec3D(0, 0, 0), 0, 0);
        insertion.setPSD(getPSD(type));
        insertion.setInitialVolume(initialVolumeFraction * getTotalVolume());
        insertion.checkBoundaryBeforeTimeStep(this);
        logger(INFO, "Inserted % particles", particleHandler.getSize());
        //write(std::cout, false);
    }

    /// switch between shear and compression stage
    void actionsAfterTimeStep() override {
        //check every 50 timesteps
        static unsigned countN = 0;
        if (++countN >= 50) countN = 0; else return;

        printTime();

        const double minTimeCompression = 0.01;
        const double minTimeShear = 0.1;
        static double minEndTimeStage = minTimeCompression;

        auto wallArea = (getXMax() - getXMin()) * (getYMax() - getYMin());
        const double pressureGoal = compressiveStress_[nStage_] * wallArea;

        if (stage_ == Stage::Compression) {
            //change to shear stage if
            // - the minimum time value is reached,
            // - the system has relaxed, and
            // - the pressure is near its goal.
            if (getTime() > minEndTimeStage
                && getKineticEnergy() < 1e-3 * getElasticEnergy()
                && fabs(shearedWall_->getForce().Z / pressureGoal - 1.0) < relTolerance_) {
                stage_ = Stage::Shear;
                minEndTimeStage = getTime() + minTimeShear;
                logger(INFO, "Switch to shear stage #%", nStage_);
                shearVelocity_ = Vec3D(shearRate_ * (shearedWall_->getPosition().Z - getZMin()), 0, 0);
                shearedWall_->setVelocityControl({0, 0, pressureGoal}, gainFactor, shearVelocity_);
            }
        } else if (stage_ == Stage::Shear) { //now the shear stage is turned on and we shear until
            //here we store last 50 shear stress values from last 500 timesteps, and check the mean, STD
            static size_t count = 0;
            const size_t num = 50; //number of values stored
            static std::array<Mdouble, num> shearStress;
            static std::array<Mdouble, num> compressiveStress;

            shearStress[count] = fabs(shearedWall_->getForce().X / wallArea);
            compressiveStress[count] = fabs(shearedWall_->getForce().Z / wallArea);
            if (count == num - 1) count = 0; else count++;

            Mdouble meanShearStress = 0;
            for (auto s : shearStress) {
                meanShearStress += s;
            }
            meanShearStress /= shearStress.size();

            Mdouble var = 0;
            for (auto s : shearStress) {
                var += mathsFunc::square(s - meanShearStress);
            }
            var /= shearStress.size();

            Mdouble relStd = sqrt(var) / fabs(meanShearStress);

            double relToleranceVary_;
            if (getTime() <= minEndTimeStage + 0.5) {
                relToleranceVary_ = 0.005;
            } else if (getTime() <= minEndTimeStage + 1) {
                relToleranceVary_ = 0.01;
            } else if (getTime() <= minEndTimeStage + 1.5) {
                relToleranceVary_ = 0.02;
            } else if (getTime() <= minEndTimeStage + 2) {
                relToleranceVary_ = 0.04;
            } else if (getTime() <= minEndTimeStage + 3) {
                relToleranceVary_ = 0.1;
            } else {
                relToleranceVary_ = 10.0;
            }

            logger(INFO, "relStd %\t relTol %\t", relStd, relToleranceVary_);

            //Change to compression stage (or end simulation) if
            // - the minimum time value is reached,
            // - the pressure is near its goal
            // - the shear stress is near its average of the last 50 values
            if (getTime()>minEndTimeStage
                && fabs(shearedWall_->getForce().Z/pressureGoal-1.0)<relTolerance_
                && relStd < relToleranceVary_)
            {
                stage_ = Stage::Compression;
                minEndTimeStage = getTime() + minTimeCompression;
                nStage_++;
                shearStress_.push_back(meanShearStress);
                // If this was the last shear stage, end the simulation, else continue with next compression stage.
                if (nStage_ == compressiveStress_.size()) {
                    setTimeMax(getTime());
                    logger(INFO, "End of simulation after % stages", nStage_);
                    return;
                } else {
                    shearedWall_->setVelocityControl({0, 0, compressiveStress_[nStage_] * wallArea}, gainFactor);
                    logger(INFO, "Switch to compression #% (p=%, relStd % < relTol %)", nStage_, compressiveStress_[nStage_],relStd, relToleranceVary_);
                }
            }
        } else {
            logger(ERROR, "This should never be called");
        }
    }

    /// write p, tau, pGoal to ene file so it can later be plotted with gnuplot
    void writeEneTimeStep(std::ostream &os) const override {
        static auto wallArea = (getXMax() - getXMin()) * (getYMax() - getYMin());
        auto forceGoal = compressiveStress_[nStage_] * wallArea;
        os << "time " << getTime()
           << " pGoal " << forceGoal
           << " p " << shearedWall_->getForce().Z
           << " tau " << fabs(shearedWall_->getForce().X)
           << " shear " << (shearedWall_->getVelocity().X!=0)
           //current volume compared with initial volume
           << " relVolume " << (shearedWall_->getPosition().Z-getZMin())/(getZMax()-getZMin())
           << " eneRatio " << -shearedWall_->getForce().X
           << '\n';
    }

    /// write p/pGoal, tau/pGoal, eneKin/eneEla to screen
    void printTime() const override {
        static auto wallArea = (getXMax() - getXMin()) * (getYMax() - getYMin());
        auto forceGoal = compressiveStress_[nStage_] * wallArea;
        logger(INFO, "time % relCompression %\trelShearStress %\teneRatio %\trelVolume %",
                getTime(),
                shearedWall_->getForce().Z/forceGoal,
                fabs(shearedWall_->getForce().X)/forceGoal,
                getKineticEnergy()/getElasticEnergy(),
                (shearedWall_->getPosition().Z-getZMin())/(getZMax()-getZMin())
                );
    }

    /// pass the stored shear stresses to the main function
    std::vector<double> getShearStress() { return shearStress_; }

private:

    Vec3D shearVelocity_;
    const Mdouble shearRate_;
    const std::vector<double> compressiveStress_;
    unsigned nStage_ = 0;
    Stage stage_ = Stage::Compression;
    std::vector<double> shearStress_;
    Mdouble relTolerance_ = 0.01;
    //Mdouble relToleranceVary_ = 0.001;
    InfiniteWall *shearedWall_ = nullptr;
    Vec3D gainFactor = Vec3D(0, 0, 1e-1);
};


int main(int argc, char *argv[]) {
    //Reading in the five parameters/arguments, which can switch material and some of material properties
    SpeciesType speciesType;
    const std::string speciesTypeString = helpers::readFromCommandLine(argc, argv, "-speciesType", std::string("APAPP"));
    // wall properties
    const Mdouble relativeCohesionStiffness = helpers::readFromCommandLine(argc, argv, "-wallRelativeCohesionStiffness",
                                                                           0.000103663);
    const Mdouble slidingFriction = helpers::readFromCommandLine(argc, argv, "-wallSlidingFriction", 0.938272);
    const Mdouble rollingFriction = helpers::readFromCommandLine(argc, argv, "-wallRollingFriction", 0.928001);

    //set speciesType
    if (speciesTypeString == "APAPP") speciesType = SpeciesType::APAPP;
    else if (speciesTypeString == "MPT") speciesType = SpeciesType::MPT;
    else if (speciesTypeString == "PH101") speciesType = SpeciesType::PH101;
    else if (speciesTypeString == "SD100") speciesType = SpeciesType::SD100;
    else logger(ERROR, "SpeciesType % not found", speciesType);

    //set six compressive stress values
    //const std::vector<double> compressiveStress = {408, 1128, 1848, 2569, 3289, 4010}; // [Pa]
    const std::vector<double> compressiveStress = {4010, 1848, 408}; // [Pa]
    const Mdouble shearRate = 0.5; // [1/s] Then shear.
    const Vec3D shearVelocity = Vec3D(0, 0, 0);

    //define the simulation
    ShearCellWall dpm(speciesType, compressiveStress, shearVelocity, relativeCohesionStiffness, slidingFriction,
                      rollingFriction, shearRate);

    //reset the output name with the new set parameter values
    std::stringstream restartName;
    restartName << std::setprecision(6) << std::scientific << "ShearCellWallFE"
                << speciesTypeString
                << "_rC" << relativeCohesionStiffness
                << "_mus" << slidingFriction
                << "_mur" << rollingFriction;
    dpm.setName(restartName.str());

    //dpm.setName("ShearCellWallFE");
    dpm.removeOldFiles();
    dpm.setParticlesWriteVTK(false);
    dpm.setWallsWriteVTK(false);
    dpm.dataFile.setFileType(FileType::NO_FILE);
    dpm.eneFile.setFileType(FileType::ONE_FILE);
    dpm.fStatFile.setFileType(FileType::NO_FILE);
    dpm.restartFile.setFileType(FileType::ONE_FILE);
    dpm.setSaveCount(1000);
    //dpm.eneFile.setSaveCount(10);
    dpm.restartFile.setSaveCount(NEVER);
    //helpers::writeToFile(dpm.getName() + ".gnu",
                         //"p [][0:1.2] 'ShearCellWallFE.ene' u 2:($6/$4) t 'p/pG' w l, '' u 2:($8/$4) t 'tau/pg' w l, '' u 2:10 t 'shear'");
    dpm.solve();
    for (int i = 0; i < compressiveStress.size(); i++) {
        logger(INFO, "Sigma % Tau % Phi %", compressiveStress[i], dpm.getShearStress()[i],
               atan(dpm.getShearStress()[i] / compressiveStress[i]) * 180. / constants::pi);
    }

    logger(INFO, "Type gnuplot %.gnu to plot stresses", dpm.getName());

    //store the final data to a string
    std::stringstream ss1;
    for (int i = 0; i < compressiveStress.size(); i++) {
        ss1 << compressiveStress[i] << " " << dpm.getShearStress()[i] << " ";
    }

    //write the final data to a file
    std::cout << "Results:\n" << ss1.str();
    helpers::writeToFile(restartName.str() + ".txt", ss1.str());

    return 0;
}
