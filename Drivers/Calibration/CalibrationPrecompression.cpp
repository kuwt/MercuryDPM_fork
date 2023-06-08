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
    PreShear = 0, //We compress the samples from low density to high density with fixed normal stress and fixed shear rate, until steady state is reached
    ShearBack = 1, //Once the steady state is reached, we shear back until shear stress is zero
    ShearBackRelaxation = 2, //Then we relax the sample with normal stress constant to reach a full relaxed state
};

void writeHelpers(const std::string& name) {
    helpers::writeToFile(
    name+".stress.py",
    "import matplotlib.pyplot as plt\n"
    "import numpy as np\n"
    "import sys\n"
    "\n"
    "#read from file\n"
    "fileName = str(sys.argv[0])[:-3]\n"
    "print(fileName)\n"
    "file = open(fileName)\n"
    "lines = file.readlines()\n"
    "print(lines[0])\n"
    "data = [line.split() for line in lines[2:]]\n"
    "data = [[float(d) for d in dat] for dat in data]\n"
    "nTime = [dat[0] for dat in data]\n"
    "solidFraction = [dat[1] for dat in data]\n"
    "normalStress = [dat[2] for dat in data]\n"
    "shearStress = [dat[3] for dat in data]\n"
    "strainRate = [dat[4] for dat in data]\n"
    "eneRatio = [dat[5] for dat in data]\n"
    "shearStage = [dat[6] for dat in data]\n"
    "\n"
    "# Plot\n"
    "fig, axs = plt.subplots(2,3)\n"
    "axs[0][0].plot(nTime, solidFraction,'-')\n"
    "axs[0][0].set_title('solidFraction')\n"
    "axs[1][0].plot(nTime, normalStress,'-')\n"
    "axs[1][0].set_title('normalStress')\n"
    "if len(nTime) >= 100: axs[1][0].set_ylim(min(normalStress[100:]),max(normalStress[100:]))\n"
    "axs[0][1].plot(nTime, shearStress,'-')\n"
    "axs[0][1].set_title('shearStress')\n"
    "if len(nTime) >= 100: axs[0][1].set_ylim(min(shearStress[100:]),max(shearStress[100:]))\n"
    "axs[1][1].plot(nTime, shearStage,'-')\n"
    "axs[1][1].set_title('shearStage')\n"
    "axs[1][2].plot(nTime, strainRate,'-')\n"
    "axs[1][2].set_title('strainRate')\n"
    "axs[0][2].semilogy(nTime, eneRatio,'-')\n"
    "axs[0][2].set_title('eneRatio')\n"
    "plt.tight_layout()\n"
    "plt.show()");
}

/**
 * This simulation consists of two parts:
 * A)
 * Particles, with given volume fracion, are
 * placed in a cubic box of dimensions Xmax, Ymax, Zmax
 * and period boundaries are applied in all directions.
 * With a given strain-rate, all 6 period boundaries are moved inwards
 * isotropically to compress the particle gas to target normal stresses.
 * B)
 * After the system is relaxed, the samples are sheared follow the procedure of
 * Schulze Ring Shear Tester, for details please see the explaination of each
 * stage.
 * C)
 * After all the stages are finished, we will have one preshear point and
 * three shear points, so on a Normal_stress/Shear_stress plane, we could
 * first find a linearized yield locus, then reconstruct the mohr's circle
 * to extract the maximum preconsolidation stress sigma_1 and
 * unconfined yield stress sigma_c, thus the FFc =  sigma_1/sigma_c
 * Note:
 * at the end of this code, user can set compressiveStress for preshear, and normal stress points for shear.
 * For normal stress shear points, you have to set ShearStageData.compressiveStress = 0.8 * compressiveStress,
 * as an example. User can also give input material name like APAPP with possiblity of changing
 * restitutionCoefficient, relativeCohesionStiffness, slidingFriction and rollingFriction.
 **/
class PreCompression : public Calibration
{
    Stage stage_ = Stage::PreShear;
    Mdouble preCompressionShearStress_ = 0;
    Mdouble preCompressionCompressiveStress_ = 0;
    Mdouble preCompressionbulkDensity_ = 0;
    const Mdouble initialVolumeFraction_ = 0.6;
    const Mdouble compressiveStress_;
    const Mdouble gain_;
    const Mdouble strainRate_;
    const bool smallSample_;
    std::ofstream fileOut_;
    StressStrainControlBoundary* boundary_;
public:

    /**
     * Constructor.
     * @param stressGoal the stress tensor that the material should be subjected to
     * @param strainRate the applied strain rate
     * @param gainFactor the sensitivity by which volume is modified to achieve the stress, or strain, goal
     * @param isStrainRateControlled if this is switched off, then stress is applied and strain is controlled
     */
    PreCompression(int argc, char* argv[]) :
        Calibration(argc,argv), //the constructor of the calibration class sets the species and psd
        compressiveStress_(helpers::readFromCommandLine(argc,argv,"-normalStress",1000)),
        gain_(helpers::readFromCommandLine<double>(argc,argv,"-gain",0.01)),
        strainRate_(helpers::readFromCommandLine<double>(argc,argv,"-inertialNumber",1e-3)
            /psd.getVolumeDx(50)*sqrt(compressiveStress_/particleSpecies->getDensity())),
        smallSample_(helpers::readFromCommandLine(argc,argv,"-smallSample"))
    {
        // write initial values to the screen
        logger(INFO,"\nStarting precompression: Goal pressure %, d50 %, strain rate %, inertial number %",
               compressiveStress_,
               psd.getVolumeDx(50),
               strainRate_,
               strainRate_*psd.getVolumeDx(50)/sqrt(compressiveStress_/particleSpecies->getDensity())
               );

        //set name
        setName("CalibrationShearCell" + param);
        
        // open file for stat output
        fileOut_.open(getName() + ".stress");
        fileOut_ << "NTime SolidFraction NormalStress ShearStress StrainRate EneRatio ShearStage\n";
        writeHelpers(getName());
        logger(INFO,"Writing statistics to: %.stress",getName());

        // write gnu script to view stress file
        helpers::writeToFile(getName() + ".stress.gnu",
            "set key autotitle columnheader\n"
            "p [1000:][] '" + getName()  + ".stress' u 1:3 w l, '' u 1:4 w l, '' u 1:($7*500) w l\n"
        );
    
        // define output file settings
        setOutput(helpers::readFromCommandLine(argc,argv,"-fullOutput"));

        // set xballs parameters
        setXBallsAdditionalArguments("-3dturn 1 -v0 -solidf");
        
        // no gravity
        setGravity({0,0,0}); //set gravity in the system

        // create cubic volume that is 15 mean particle diameters side length at a volume fraction of 0.6
        double vf = cbrt(0.6 / initialVolumeFraction_); //so the rve is applicable to 0.6 vloume fraction
        Mdouble lengthRVE = std::max(vf*15.0*psd.getVolumeDx(50),2.0*2.0*psd.getMaxRadius());
        // (use size of 10 mean particle diameters if small sample)
        if (smallSample_) lengthRVE = vf*5.0*psd.getVolumeDx(50);
        logger(INFO,"Goal RVE size = % D50", helpers::round(lengthRVE/psd.getVolumeDx(50)/vf,1));
        setMax(lengthRVE*Vec3D(0.5,0.5,0.5));
        setMin(-getMax());
        
        // simulate forever (will be terminated in actionsAfterTimeStep)
        setTimeMax(1000);

        // set stressStrainControlBoundary based on user input
        logger(INFO,"\nStage 0 (preshear): shearing with controlled compressive stress % at strain rate %",compressiveStress_, strainRate_);
        StressStrainControlBoundary ssc;
        Matrix3D stressGoal = Matrix3D(1,0,0, 0,1,0, 0,0,1) * compressiveStress_;
        Matrix3D strainRate = Matrix3D(0,1,0,0,0,0,0,0,0) * strainRate_;
        Matrix3D gainFactor = Matrix3D(1,0,0, 0,1,0, 0,0,1) * gain_; //a bit more gain than usual
        bool isStrainRateControlled = true;
        ssc.setHandler(&boundaryHandler);
        ssc.set(stressGoal, strainRate, gainFactor, isStrainRateControlled);
        boundary_ = boundaryHandler.copyAndAddObject(ssc);
    }

    void setupInitialConditions() override {
        // add particles
        CubeInsertionBoundary insertion;
        SphericalParticle particle(speciesHandler.getObject(0));
        insertion.set(&particle, 100000, getMin(), getMax(), Vec3D(0, 0, 0), Vec3D(0, 0, 0), 0, 0);
        insertion.setPSD(psd);
        insertion.setInitialVolume(initialVolumeFraction_ * getTotalVolume());
        insertion.setCheckParticleForInteraction(false);
        insertion.checkBoundaryBeforeTimeStep(this);
        logger(INFO, "Inserted % particles (volfrac % of %)", particleHandler.getSize(),particleHandler.getVolume()/getTotalVolume(), initialVolumeFraction_);
        writePSDToFile();
    }

    void printTime() const override {
        if (getNumberOfTimeSteps()%20>0) return;
        const Matrix3D stress = getTotalStress();
        static Time time;
        std::stringstream out;
        out << std::fixed << std::setprecision(2)
            << "nTime " << getNumberOfTimeSteps()/1000 << "k"
            << " solidFr " << particleHandler.getVolume()/getTotalVolume()
            << " ene " << getKineticEnergy() / getElasticEnergy()
            << " pressure " << stress.YY
            << " shearStress " << stress.XY
            << " strainRate " << boundary_->getStrainRate().XY
        << " cpu " << time.toc()
        << " stage " << (int)stage_;
        logger(INFO, "%",out.str());
        time.tic();
    }
    
    // compute mean
    template <size_t S>
    Mdouble computeMean(std::array<Mdouble,S> vec) {
        Mdouble mean = 0;
        for (auto val : vec) mean += val;
        mean /= vec.size();
        return mean;
    }
    
    // compute standard deviation
    template <size_t S>
    Mdouble computeStd(std::array<Mdouble,S> vec, Mdouble mean) {
        Mdouble variance = 0;
        for (auto val : vec) variance += mathsFunc::square(val - mean);
        variance /= vec.size();
        return sqrt(variance);
    }
    
    // compute coefficient of variation and mean
    template <size_t S>
    std::pair<Mdouble,Mdouble> computeMeanAndCov(std::array<Mdouble,S> vec) {
        const Mdouble mean = computeMean(vec);
        const Mdouble std = computeStd(vec,mean);
        const Mdouble cov = std / fabs(mean);
        return {mean, cov};
    }
    
    void actionsAfterTimeStep() override
    {
        // check only every 100 timpsteps
        const unsigned step = 200;
        if (getNumberOfTimeSteps()%step) return;

        // write to .stress file
        const Matrix3D stress = getTotalStress();
        std::stringstream out;
        out << getNumberOfTimeSteps() << ' '
                 << getTotalMass()/speciesHandler.getLastObject()->getDensity()/getTotalVolume() << ' '
                 << stress.YY << ' '
                 << stress.XY << ' '
                 << boundary_->getStrainRate().XY << ' '
            << getKineticEnergy()/getElasticEnergy() << ' '
            << (int) stage_ << std::endl;
        fileOut_ << out.str();
        //std::cout << out.str();
    
        // when checking for steadiness, we average over a number of time steps (nSteps) uniformly distributed over a time interval (nTimeStepsAveraged)
        const size_t nTimeStepsAveraged = 80000; //making the shear stage twice as long (should be even longer)
        const size_t nSteps = nTimeStepsAveraged / step;

        // Start with preshear stage
        if (stage_==Stage::PreShear) {
            static std::array <Mdouble,nSteps> shearStress; //store the last shear stress values
            static std::array <Mdouble,nSteps> compressiveStress; //store the last compressive stress values
            
            //track how often this function has been called
            static size_t countPreShear = 0;
            //store the current shear and compressive stress
            shearStress[countPreShear % nSteps] = stress.XY;
            compressiveStress[countPreShear % nSteps] = stress.YY;
            countPreShear++;

            // compute mean shear stress
            auto [meanShearStress, covShearStress] = computeMeanAndCov(shearStress);

            // compute variance pressure
            auto [meanCompressiveStress, covCompressiveStress] = computeMeanAndCov(compressiveStress);
            
            // parameter: tolerance for CoV (starts as 3%, increases with time)
            const Mdouble covTolerance = (smallSample_?0.2:0.05) * ((double)getNumberOfTimeSteps() / (double)nTimeStepsAveraged);

            // end preShear stage if cov falls below tolerance
            if (countPreShear >= nSteps && //run for at least nSteps
                covCompressiveStress < covTolerance && //check that compressive stress is stable
                covShearStress < covTolerance) //check that shear stress is stable
            {
                // shear back; constant small negative shear rate, controlled pressure
                logger(INFO, "Shear stress steady (cov=%<tol=%). Pressure steady (cov=%<tol=%)."
                             "Mean compressiveStress %, mean shearStress %\n",
                       covShearStress, covTolerance, covCompressiveStress, covTolerance,
                       preCompressionCompressiveStress_, preCompressionShearStress_);
                
                // next stage
                stage_ = Stage::ShearBack;
                
                // store these values (for output and used in stage 2 end criterium)
                preCompressionShearStress_ = fabs(meanShearStress);
                preCompressionCompressiveStress_ = meanCompressiveStress;
                preCompressionbulkDensity_ = getTotalMass()/getTotalVolume();
                
                // apply a small negative strain rate
                Matrix3D stressGoal = Matrix3D(0, 0,0, 0,1,0, 0,0,0) * compressiveStress_;
                Matrix3D strainRate = Matrix3D(0,-0.2,0, 0,0,0, 0,0,0) * strainRate_;
                Matrix3D gainFactor = Matrix3D(0, 0,0, 0,1,0, 0,0,0) * gain_;
                bool isStrainRateControlled = true;
                boundary_->set(stressGoal, strainRate, gainFactor, isStrainRateControlled);

                // shear back; constant small negative shear rate, controlled pressure
                
                logger(INFO, "\nStage 1 (shearBack): shearing at strain rate % with controlled YY-stress %",strainRate.XY, stressGoal.YY);
            } else if (countPreShear % 20 == 0) {
                if (countPreShear >= nSteps) {
                    logger(INFO, "PreShear continues: cov=(shear %, compressive %) > tol=%", covShearStress, covCompressiveStress, covTolerance);
                } else {
                    logger(INFO, "PreShear continues: #time steps % < averaging interval %", countPreShear, nSteps);
                }
            }
        } else if (stage_==Stage::ShearBack) {
            // start next stage if shear stress reaches zero
            if (stress.XY>0) {
                // next stage
                stage_ = Stage::ShearBackRelaxation;
                // apply no strain, vanishing shear stress goal
                Matrix3D stressGoal = Matrix3D(0,1e-8,0, 0,1,0, 0,0,0) * compressiveStress_;
                Matrix3D strainRate = Matrix3D(0,0,0, 0,0,0, 0,0,0) * strainRate_;
                Matrix3D gainFactor = Matrix3D(0,1,0, 0,1,0, 0,0,0) * gain_;
                bool isStrainRateControlled = true;
                boundary_->set(stressGoal, strainRate, gainFactor, isStrainRateControlled);
                logger(INFO, "\nStage 2 (shearBackRelaxation): controlled XY and YY-stress % %",stressGoal.XY, stressGoal.YY);
            }
        } else if (stage_==Stage::ShearBackRelaxation) {
            // track how many time steps since beginning of shear relaxation
            static size_t nTime0 = getNumberOfTimeSteps();
            size_t nTime = getNumberOfTimeSteps()-nTime0;
            // track how often this function has been called
            static size_t count = 0;
            //store the last shear stress values
            static std::array <Mdouble,nSteps> shearStress;
            shearStress[count % nSteps] = stress.XY;
            count++;
    
            // compute variance shear stress
            Mdouble varianceShearStress = 0;
            for (auto s : shearStress) varianceShearStress += mathsFunc::square(s);
            varianceShearStress /= shearStress.size();
            // coefficient of variation shear stress
            Mdouble covShearStress = sqrt(varianceShearStress) / preCompressionShearStress_;
            
            // tolerance for CoV (starts as 5%, increases with time)
            Mdouble covTolerance = (smallSample_?0.2:0.05) * ((double)nTime / (double)nTimeStepsAveraged);

            //extra relaxation to get the final sample before shear
            if (count >= nSteps && //run for at least nSteps
                covShearStress < covTolerance) { //check that compressive stress is stable
                setTimeMax(getTime());
                logger(INFO,"Relaxation completed. Ready for Stage 3: shear.\n"
                            "Compressive stress %, shear stress %", stress.YY, stress.XY);
            } else if(nTime%nTimeStepsAveraged==0) {
                if (count < nSteps) {
                    logger(INFO, "ShearBackRelaxation continues: #time steps % < averaging interval %", count, nSteps);
                } else {
                    logger(INFO, "ShearBackRelaxation continues: cov=% > tol=%", covShearStress, covTolerance);
                }
            }
        }
        else {
            logger(ERROR,"This should never be called");
        }
        
        printTime();
    }
    
    //write output values
    void actionsAfterSolve() override {
        logger(INFO,"Pre-compression shear stress %, Pre-compression bulk density %", preCompressionShearStress_,
               preCompressionbulkDensity_);
        helpers::writeToFile(getName()+".out",
                             "preCompressionShearStress " + helpers::toString(preCompressionShearStress_) + "\n"
            "preCompressionBulkDensity " + helpers::toString(preCompressionbulkDensity_)
            );
    }
};

int main(int argc, char* argv[])
{
    //Run pre-compression simulation (compression, preshear, shear and shearBack stages)
    PreCompression preCompression(argc,argv);
    preCompression.solve();
    return 0;
}