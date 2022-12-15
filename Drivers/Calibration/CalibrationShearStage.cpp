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
#include <Mercury3D.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Boundaries/StressStrainControlBoundary.h>
#include <Boundaries/CubeInsertionBoundary.h>
#include "Boundaries/LeesEdwardsBoundary.h"
#include "Calibration.h"
#include <iomanip>
#include <utility>
#include <cmath>
#include <MercuryTime.h>
using mathsFunc::square;
using mathsFunc::cubic;
using constants::pi;

enum class Stage {
    Decompress = 0, //Now we read in the presheared-relaxed sample and decompress the the normal stress for shear point/incipient flow
    Shear = 1 //Once the normal stress is reached, we do shear and find the peak failure value of shear stress
};

class ShearStage : public Calibration
{
    Mdouble maxShearStress_ = 0;
    Mdouble bulkDensity_;
    Mdouble timeOfMaxShearStress_ = 0;
    Stage stage_ = Stage::Decompress;
    Mdouble compressiveStress_;
    const Mdouble strainRate_;
    Mdouble gain_;
    std::ofstream fileOut_;
    StressStrainControlBoundary* boundary_;
    
public:
    
    /**
     *  Here we start the second part which reads in the presheared sample and do 3 shear points
     * @param stressGoal the stress tensor that the material should be subjected to
     * @param strainRate the applied strain rate
     * @param gainFactor the sensitivity by which volume is modified to achieve the stress, or strain, goal
     * @param isStrainRateControlled if this is switched off, then stress is applied and strain is controlled
     */
    ShearStage(int argc, char* argv[]) : Calibration (argc, argv),
        compressiveStress_(helpers::readFromCommandLine(argc,argv,"-compression",1000)),
        strainRate_(helpers::readFromCommandLine<double>(argc,argv,"-inertialNumber",1e-3)/psd.getVolumeDx(50)*sqrt(compressiveStress_/particleSpecies->getDensity())),
        gain_(helpers::readFromCommandLine<double>(argc,argv,"-gain",0.01))
    {
        //restart from Precompression
        setName("CalibrationShearCell" + param);
        readRestartFile();
        setRestarted(false);
        
        // set new name
        setName(getName()+"-"+std::to_string((int)compressiveStress_)+"Pa");

        // open .stress file
        fileOut_.open(getName()+".stress");
        fileOut_ << "NTime SolidFraction NormalStress ShearStress StrainRate EneRatio ShearStage\n";
        logger(INFO,"Writing statistics to: %.stress",getName());
    
        // write gnu script to view stress file
        helpers::writeToFile(getName() + ".stress.gnu",
                             "set key autotitle columnheader\n"
                             "p [1000:][] '" + getName()  + ".stress' u 1:3 w l, '' u 1:4 w l, '' u 1:($7*500) w l\n"
        );

        // xballs arguments
        setXBallsAdditionalArguments("-3dturn 1 -v0 -solidf");
        
        //Set up the new boundaries based on user input
        logger(INFO,"\nPhase 0: Decompression to % Pa, strainRate %",compressiveStress_, strainRate_);
        boundary_ = dynamic_cast<StressStrainControlBoundary*>(boundaryHandler.getLastObject());
        logger.assert_debug(boundary_,"StressStrainControlBoundary not found");
        Matrix3D stressGoal = Matrix3D(0,1e-8,0, 0,1,0, 0,0,0) * compressiveStress_;
        Matrix3D strainRate = boundary_->getStrainRate();
        Matrix3D gainFactor = Matrix3D(0,1,0, 0,1,0, 0,0,0) * gain_;
        bool isStrainRateControlled = true;
        boundary_->set(stressGoal, strainRate, gainFactor, isStrainRateControlled);
        //logger(INFO,"Integrated shift %",boundary_->getIntegratedShift());
    }
    
    void printTime() const override {
        if (getNumberOfTimeSteps()%20) return;
        const Matrix3D stress = getTotalStress();
        logger(INFO, " time % volFraction % EneRatio % pressure % shearStress % strainRate % maxShearStress %", getTime(),
               getTotalMass()/speciesHandler.getLastObject()->getDensity()/getTotalVolume(),
               getKineticEnergy() / getElasticEnergy(),
               stress.YY,
               stress.XY,
               boundary_->getStrainRate().XY,
               maxShearStress_
        );
    }
    
    void actionsAfterTimeStep() override
    {
        static unsigned timeOfSwitch = 0;
        const Matrix3D stress = getTotalStress();
        Mdouble shearStress = std::fabs(stress.XY);
        // write to file
        fileOut_ << getNumberOfTimeSteps() << ' '
                << particleHandler.getVolume()/getTotalVolume() << ' '
                << stress.YY << ' '
                << stress.XY << ' '
                << boundary_->getStrainRate().XY << ' '
                << getKineticEnergy()/getElasticEnergy() << ' '
                << (int) stage_ << std::endl;
        // find max shear stress
        if (stress.XY < 0 && shearStress > maxShearStress_) {
            timeOfMaxShearStress_ = getTime();
            maxShearStress_ = shearStress;
            bulkDensity_ = getTotalMass()/getTotalVolume();
        }
        
        //check every 50 timesteps
        static unsigned count = 0;
        if (++count>490) count=0; else return;
        //First stage: decompress to reach a lower normal stress at shear stage
        if (stage_==Stage::Decompress) {
            //if low kinetic energy and near goal compressive stress, go to shear stage
            if (getNumberOfTimeSteps() > 5000
                && fabs(stress.YY-compressiveStress_)<0.01*compressiveStress_) {
                //Shear: constant shear, controlled pressure
                Matrix3D stressGoal = Matrix3D(0,0,0, 0,1,0, 0,0,0) * compressiveStress_;
                Matrix3D strainRate = Matrix3D(0,1,0, 0,0,0, 0,0,0) * strainRate_;
                Matrix3D gainFactor = Matrix3D(0,0,0, 0,1,0, 0,0,0) * gain_;
                bool isStrainRateControlled = true;
                boundary_->set(stressGoal, strainRate, gainFactor, isStrainRateControlled);
                logger(INFO,"\nPhase 1: Switching from decompress to shear stage");
                stage_ = Stage::Shear;
                timeOfSwitch = getNumberOfTimeSteps();
                maxShearStress_ = 0.0;
            }
        } else if (stage_==Stage::Shear) { //now the shear stage is turned on and we shear until the peak value is detected
            //if enough time has elapsed since both decompression (20k time steps) and reaching the max shear stress (+25% of time), stop
            if (getTime() > 1.5*timeOfMaxShearStress_ && getNumberOfTimeSteps() > timeOfSwitch + 20000) {
                logger(INFO, "Time % %", getTime(), timeOfMaxShearStress_);
                setTimeMax(getTime());
            }
        } else {
            logger(ERROR,"This should never be called");
        }
    }
    
    Mdouble getMaxShearStress() const {
        return maxShearStress_;
    }
    
    //write output values
    void actionsAfterSolve() override {
        logger(INFO,"Max shear stress %, bulk density %", maxShearStress_, bulkDensity_);
        helpers::writeToFile(getName()+".out",
                             "maxShearStress " + helpers::to_string(maxShearStress_) + "\n"
                                      "bulkDensity " + helpers::to_string(bulkDensity_)
        );
    }
};

int main(int argc, char* argv[]) {
    ShearStage dpm(argc,argv);
    dpm.solve();
    return 0;
}