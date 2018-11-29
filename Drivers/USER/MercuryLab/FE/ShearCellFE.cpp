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

enum class Stage {
    Compression, //here we compress the samples from low density to high density with fixed normal stress sigma_yy = experimental pre-shear normal stress
    PreShear, //Then we keep the normal stress constant and shear until the steady state
    ShearBack, //Once the steady state is reached, we shear back until shear stress is zero
    ShearBackRelaxation, //Then we relax the sample with normal stress constant to reach a full relaxed state
    Decompress, //Now we read in the pre-sheared-relaxed sample and decompress the the normal stress for shear point/incipient flow
    Shear //Once the normal stress is reached, we do shear and find the peak failure value of shear stress
};

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
 * After all the stages are finished, we will have one pre-shear point and
 * three shear points, so on a Normal_stress/Shear_stress plane, we could
 * first find a linearized yield locus, then reconstruct the mohr's circle
 * to extract the maximum preconsolidation stress sigma_1 and
 * unconfined yield stress sigma_c, thus the FFc =  sigma_1/sigma_c
 * Note:
 * at the end of this code, user can set compressiveStress for pre-shear, and normal stress points for shear.
 * For normal stress shear points, you have to set ShearStageData.compressiveStress = 0.8 * compressiveStress,
 * as an example. User can also give input material name like APAPP with possiblity of changing
 * restitutionCoefficient, relativeCohesionStiffness, slidingFriction and rollingFriction.
 **/
class StressStrainControl : public Mercury3D
{
public:
    
    /**
     *
     * @param stressGoal the stress tensor that the material should be subjected to
     * @param strainRate the applied strain rate
     * @param gainFactor the sensitivity by which volume is modified to achieve the stress, or strain, goal
     * @param isStrainRateControlled if this is switched off, then stress is applied and strain is controlled
     */
    StressStrainControl(SpeciesType type, const Mdouble compressiveStress, const Mdouble& shearRate)
     :  shearRate_(shearRate), stage_(Stage::Compression), compressiveStress_(compressiveStress), preShearShearStress_(0), preShearCompressiveStress_(0),
        relTolerance_(0.001)
    {
        const Mdouble initialVolumeFraction = 0.25;
        setXBallsAdditionalArguments("-3dturn 1 -v0 -solidf");

        //add species
        auto psd = getPSD(type);
        logger.assert(hasPSD(type),"Material type % not allowed",type);
        auto species = addSingleSpecies(type, getMedianParticleRadius(psd), getMinParticleRadius(psd), *this, true);
        logger(INFO, "%", *species);
    
        setName("ShearCellFE");
        setGravity({0,0,0}); //set gravity in the system
        // create initial cubic volume of 20 mean particle diameters side length
        setMax(25.0*getMedianParticleRadius(psd)*Vec3D(1,1,1));
        setMin(-getMax());
        setTimeMax(100);
        //setTimeStep(1e6);
        setSaveCount(1000);
        
        //Set up the new boundaries based on user input
        StressStrainControlBoundary boundary;
        
        //Here we setup the stressstraincontrol boundary with no shear, only compression
        Matrix3D stressGoal = Matrix3D(1,0,0, 0,1,0, 0,0,1) * compressiveStress;
        Matrix3D strainRate = Matrix3D(0,0,0,0,0,0,0,0,0);
        Matrix3D gainFactor = Matrix3D(1,0,0, 0,1,0, 0,0,1) * 1e-2;
        bool isStrainRateControlled = true;
        boundary.setHandler(&boundaryHandler);
        boundary.set(stressGoal, strainRate, gainFactor, isStrainRateControlled);
        boundaryHandler.copyAndAddObject(boundary);

        //add particles
        CubeInsertionBoundary insertion;
        BaseParticle particle(speciesHandler.getObject(0));
        insertion.set(&particle, 100000, getMin(), getMax(), Vec3D(0, 0, 0), Vec3D(0, 0, 0), 0, 0);
        insertion.setPSD(getPSD(type));
        insertion.setInitialVolume(initialVolumeFraction*getTotalVolume());
        insertion.checkBoundaryBeforeTimeStep(this);
        logger(INFO,"Inserted % particles",particleHandler.getSize());
        //write(std::cout, false);
    }
    
    void printTime() const override {
        logger(INFO, "time %\tsolidFraction %\tEneRatio %\tpressure %\tshearStress %\tshearRate %\t", getTime(),
                getTotalMass()/speciesHandler.getLastObject()->getDensity()/getTotalVolume(),
                getKineticEnergy() / getElasticEnergy(),
                getTotalStress().YY,
                getTotalStress().XY,
                dynamic_cast<const StressStrainControlBoundary*>(boundaryHandler.getLastObject())->getStrainRate().XY
                );
    }
    
    void actionsAfterTimeStep() override
    {
        //check only every 10 timpsteps
        static unsigned count = 0;
        if (++count>9) count=0; else return;
        //Initialize the compression stage
        if (stage_==Stage::Compression) {
            //stop compression stage after the bulk is relaxed
            if (getKineticEnergy()<0.02*getElasticEnergy() && fabs(getTotalStress().YY-compressiveStress_)<0.01*compressiveStress_)
            {
                stage_ = Stage::PreShear;
                logger(INFO,"Switching from compression to preShear stage");
                //Set up the new boundaries based on user input
                auto boundary = dynamic_cast<StressStrainControlBoundary*>(boundaryHandler.getLastObject());
                logger.assert(boundary,"StressStrainControlBoundary not found");
    
                Matrix3D stressGoal = Matrix3D(0,0,0, 0,1,0, 0,0,0) * compressiveStress_;
                Matrix3D strainRate = Matrix3D(0,1,0, 0,0,0, 0,0,0) * shearRate_;
                Matrix3D gainFactor = Matrix3D(0,0,0, 0,1,0, 0,0,0) * 1;
                bool isStrainRateControlled = true;
                boundary->set(stressGoal, strainRate, gainFactor, isStrainRateControlled);
            }
        } else if (stage_==Stage::PreShear) { //switch to pre-shear stage after compression
            //here we store last 50 shear stress values from last 500 timesteps, and check the mean, STD
            static size_t count = 0;
            const size_t num = 50; //number of values stored
            static std::array <Mdouble,num> shearStress;
            static std::array <Mdouble,num> compressiveStress;
            
            Matrix3D stress = getTotalStress();
            shearStress[count] = stress.XY;
            compressiveStress[count] = stress.YY;
            if (count==num-1) count=0; else count++;

            Mdouble meanShearStress = 0;
            for (auto s : shearStress) {
                meanShearStress += s;
            }
            meanShearStress /= shearStress.size();
    
            Mdouble var = 0;
            for (auto s : shearStress) {
                var += mathsFunc::square(s-meanShearStress);
            }
            var /= shearStress.size();
            
            Mdouble relStd = sqrt(var)/fabs(meanShearStress);

            //Here the tolerance increases with time
            if (getTime() <= 2) {
                relTolerance_= 0.001;
            }else if (getTime() > 2 && getTime() <= 3) {
                relTolerance_= 0.005;
            }else if (getTime() > 2 && getTime() <= 3) {
                relTolerance_= 0.01;
            }else if (getTime() > 3 && getTime() <= 4) {
                relTolerance_= 0.04;
            }else if (getTime() > 4 && getTime() <= 5) {
                relTolerance_= 0.08;
            }else {
                relTolerance_= 0.1;
            }
            //logger(INFO,"relative std dev %\t relTolerance %\t", relStd, relTolerance_);
            
            // stop simulation if relative std deviation of shear stress falls below 5%
            if (relStd<relTolerance_) {
                stage_ = Stage::ShearBack;
                //store the mean shear stress at the end of the preshear stage
                preShearShearStress_ = meanShearStress;
                //store the mean compressive stress at the end of the preshear stage
                preShearCompressiveStress_ = 0;
                for (auto c : compressiveStress) {
                    preShearCompressiveStress_ += c;
                }
                preShearCompressiveStress_ /= compressiveStress.size();

                logger(INFO,"Switching from preShear to shearBack stage (compressiveStress % shearStress %)",preShearCompressiveStress_,preShearShearStress_);
                Matrix3D stressGoal = Matrix3D(0, 0,0, 0,1,0, 0,0,0) * compressiveStress_;
                Matrix3D strainRate = Matrix3D(0,-0.2,0, 0,0,0, 0,0,0) * shearRate_;
                Matrix3D gainFactor = Matrix3D(0, 0,0, 0,1,0, 0,0,0) * 1;
                bool isStrainRateControlled = true;
                auto boundary = dynamic_cast<StressStrainControlBoundary*>(boundaryHandler.getLastObject());
                logger.assert(boundary,"StressStrainControlBoundary not found");
                boundary->set(stressGoal, strainRate, gainFactor, isStrainRateControlled);
            }
        } else if (stage_==Stage::ShearBack) { //we shear back to shear stress reaches close to zero
            if (getTotalStress().XY>0) {
                stage_ = Stage::ShearBackRelaxation;
                logger(INFO,"Switching from ShearBack steage to ShearBackRelaxation stage");
                Matrix3D stressGoal = Matrix3D(0,0.001,0, 0,1,0, 0,0,0) * compressiveStress_;
                Matrix3D strainRate = Matrix3D(0,0,0, 0,0,0, 0,0,0) * shearRate_;
                Matrix3D gainFactor = Matrix3D(0,1,0, 0,1,0, 0,0,0) * 1;
                bool isStrainRateControlled = true;
                auto boundary = dynamic_cast<StressStrainControlBoundary*>(boundaryHandler.getLastObject());
                logger.assert(boundary,"StressStrainControlBoundary not found");
                boundary->set(stressGoal, strainRate, gainFactor, isStrainRateControlled);
                logger(INFO,"Entering the relaxation of shearBack stage (compressiveStress % shearStress %)",preShearCompressiveStress_,preShearShearStress_);
            }
        } else if (stage_==Stage::ShearBackRelaxation)
        {
            if (fabs(getTotalStress().YY-compressiveStress_)<0.001*compressiveStress_) { //extra relaxation to get the final sample before shear
                setTimeMax(getTime());
                logger(INFO,"Ready for shear stages");
                logger(INFO, "Relaxation of shearBack stage (compressiveStress % shearStress %)",
                       getTotalStress().YY, getTotalStress().XY);
            }
        }
        else {
            logger(ERROR,"This should never be called");
        }

    }
    //Two functions to pass the stored preshear point to the main function
    Mdouble getPreShearShearStress() {return preShearShearStress_;}
    Mdouble getPreShearCompressiveStress() {return preShearCompressiveStress_;}

private:

    const Mdouble shearRate_;
    Stage stage_;
    Mdouble compressiveStress_;
    Mdouble preShearShearStress_;
    Mdouble preShearCompressiveStress_;
    Mdouble relTolerance_;
};

class ShearStage : public Mercury3D
{
public:

    /**
     *  Here we start the second part which reads in the pre-sheared sample and do 3 shear points
     * @param stressGoal the stress tensor that the material should be subjected to
     * @param strainRate the applied strain rate
     * @param gainFactor the sensitivity by which volume is modified to achieve the stress, or strain, goal
     * @param isStrainRateControlled if this is switched off, then stress is applied and strain is controlled
     */
    ShearStage(SpeciesType type, const Mdouble compressiveStress, const Mdouble& shearRate, const std::string& newName)
            :  shearRate_(shearRate), compressiveStress_(compressiveStress), maxShearStress_(0), timeOfMaxShearStress_(0), stage_(Stage::Decompress)
    {
        readRestartFile(newName);
        setRestarted(false);
        std::stringstream ss1;
        ss1 << std::fixed << std::setprecision(0) << compressiveStress;
        const auto newStressName = ss1.str();
        setName(newName+"-"+newStressName+"Pa");
        logger(INFO,"Name %",getName());

        setXBallsAdditionalArguments("-3dturn 1 -v0 -solidf");
        //Set up the new boundaries based on user input
        auto boundary = dynamic_cast<StressStrainControlBoundary*>(boundaryHandler.getLastObject());
        logger.assert(boundary,"StressStrainControlBoundary not found");
        Matrix3D stressGoal = Matrix3D(0,0.001,0, 0,1,0, 0,0,0) * compressiveStress_;
        Matrix3D strainRate = Matrix3D(0,0,0, 0,0,0, 0,0,0) * shearRate_;
        Matrix3D gainFactor = Matrix3D(0,1,0, 0,1,0, 0,0,0) * 1;
        bool isStrainRateControlled = true;
        boundary->set(stressGoal, strainRate, gainFactor, isStrainRateControlled);
    }

    void printTime() const override {
        logger(INFO, "time % %\tsolidFraction %\tEneRatio %\tpressure %\tshearStress %\tshearRate %\t", getTime(), getTimeMax(),
               getTotalMass()/speciesHandler.getLastObject()->getDensity()/getTotalVolume(),
               getKineticEnergy() / getElasticEnergy(),
               getTotalStress().YY,
               getTotalStress().XY,
               dynamic_cast<const StressStrainControlBoundary*>(boundaryHandler.getLastObject())->getStrainRate().XY
        );
        logger(INFO,"maxShearStress %",maxShearStress_);
    }

    void actionsAfterTimeStep() override
    {
        static unsigned timeOfSwitch = 0;

        // write out max shear stress
        Mdouble shearStress = std::fabs(getTotalStress().XY);
        if (getTotalStress().XY < 0 && shearStress > maxShearStress_) {
            timeOfMaxShearStress_ = getTime();
            maxShearStress_ = shearStress;
        }
        //check every 50 timesteps
        static unsigned count = 0;
        if (++count>49) count=0; else return;
        //initialize the decompression to reach a lower normal stress at shear stage
        if (stage_==Stage::Decompress) {
            if (getKineticEnergy()<1e-3*getElasticEnergy() && fabs(getTotalStress().YY-compressiveStress_)<0.001*compressiveStress_) {
                //Set up the new boundaries based on user input
                auto boundary = dynamic_cast<StressStrainControlBoundary*>(boundaryHandler.getLastObject());
                logger.assert(boundary,"StressStrainControlBoundary not found");
                Matrix3D stressGoal = Matrix3D(0,0,0, 0,1,0, 0,0,0) * compressiveStress_;
                Matrix3D strainRate = Matrix3D(0,1,0, 0,0,0, 0,0,0) * shearRate_;
                Matrix3D gainFactor = Matrix3D(0,0,0, 0,1,0, 0,0,0) * 1;
                bool isStrainRateControlled = true;
                boundary->set(stressGoal, strainRate, gainFactor, isStrainRateControlled);
                logger(INFO,"Switching from decompress to shear stage");
                stage_ = Stage::Shear;
                timeOfSwitch = getNumberOfTimeSteps();
                maxShearStress_ = 0.0;
            }
        } else if (stage_==Stage::Shear) { //now the shear stage is turned on and we shear until the peak value is detected
            //stopping criterium: timeOfMaxShearStress<0.8*getTime()
            if (timeOfMaxShearStress_ < 0.8 * getTime() && getNumberOfTimeSteps() > timeOfSwitch + 20000) {
                logger(INFO, "Time % %", getTime(), timeOfMaxShearStress_);
                setTimeMax(getTime());
            }
        } else {
            logger(ERROR,"This should never be called");
        }
    }

    Mdouble getMaxShearStress() {return maxShearStress_;}

private:

    const Mdouble shearRate_;
    Mdouble compressiveStress_;
    Mdouble maxShearStress_;
    Mdouble timeOfMaxShearStress_;
    Stage stage_;
};

// Here is a templated function to do the linear regression fit of data points.
template <typename T>
std::vector<T> getLinearFit(const std::vector<T>& Xdata, const std::vector<T>& Ydata)
{
    T xSum = 0, ySum = 0, xxSum = 0, xySum = 0, slope, intercept;
    

    for (long i = 0; i < Ydata.size(); i++)
    {
        xSum += Xdata[i];
        ySum += Ydata[i];
        xxSum += Xdata[i] * Xdata[i];
        xySum += Xdata[i] * Ydata[i];
    }
    slope = (Ydata.size() * xySum - xSum * ySum) / (Ydata.size() * xxSum - xSum * xSum);
    intercept = (ySum - slope * xSum) / Ydata.size();
    std::vector<T> res;
    res.push_back(slope);
    res.push_back(intercept);
    return res;
}

//Here is a function to solve the quadratic equation
std::pair<double , double > solveQuad(double a, double b, double c)
{
    std::pair<double, double> roots; //The pair which will contain the results.
    //It is set to hold two floats.
    roots.first = ((-1*b)+sqrt((b*b)-4*a*c))/(2.0*a); //Set the first root
    roots.second = ((-1*b)-sqrt((b*b)-4*a*c))/(2.0*a); //Set the second root
    return roots;
}

int main(int argc, char* argv[])
{
    //first simulation runs the compression, preshear and shearBacl stages
    const Mdouble compressiveStress = 1e3; // [Pa] First compress.
    const Mdouble shearRate = 0.5; // [1/s] Then shear.
    
    //Reading in the five parameters/arguments, which can switch material and some of material properties
    SpeciesType speciesType;
    const std::string speciesTypeString = helpers::readFromCommandLine(argc,argv,"-speciesType",std::string("SD100"));
    if (speciesTypeString=="APAPP") speciesType = SpeciesType::APAPP;
    else if (speciesTypeString=="MPT") speciesType = SpeciesType::MPT;
    else if (speciesTypeString=="PH101") speciesType = SpeciesType::PH101;
    else if (speciesTypeString=="SD100") speciesType = SpeciesType::SD100;
    else logger(ERROR, "SpeciesType % not found",speciesType);
    const Mdouble restitutionCoefficient = helpers::readFromCommandLine(argc,argv,"-restitutionCoefficient",0.5);
    const Mdouble relativeCohesionStiffness = helpers::readFromCommandLine(argc,argv,"-relativeCohesionStiffness",0.1);
    const Mdouble slidingFriction = helpers::readFromCommandLine(argc,argv,"-slidingFriction",0.5);
    const Mdouble rollingFriction = helpers::readFromCommandLine(argc,argv,"-rollingFriction",0.1);
    
    //reset the output name with the new set parameter values
    std::stringstream restartName;
    //restartName << std::fixed << std::setprecision(5) << "ShearCellFE" << "_e" << restitutionCoefficient << "_rC" << relativeCohesionStiffness << "_mus" << slidingFriction << "_mur" << rollingFriction;
    restartName << std::setprecision(6) << std::scientific << "ShearCellFE" << "_e" << restitutionCoefficient << "_rC" << relativeCohesionStiffness << "_mus" << slidingFriction << "_mur" << rollingFriction;
    const auto newName = restartName.str();
    std::cout << newName << std::endl;
    //Now reset/modify the default species to the one user inputs
    StressStrainControl dpm(speciesType, compressiveStress, shearRate);
    modifySpecies(dpm.speciesHandler.getLastObject(), restitutionCoefficient, relativeCohesionStiffness, slidingFriction, rollingFriction, speciesType);
    dpm.setName(newName);
    logger(INFO,"Species: %", *dpm.speciesHandler.getLastObject());
    
    dpm.setParticlesWriteVTK(false);
    dpm.setWallsWriteVTK(false);
    dpm.dataFile.setFileType(FileType::NO_FILE);
    dpm.eneFile.setFileType(FileType::NO_FILE);
    dpm.fStatFile.setFileType(FileType::NO_FILE);
    dpm.restartFile.setFileType(FileType::ONE_FILE);
    dpm.restartFile.setSaveCount(100000000);
    dpm.solve();

    //then we run three shear simulations at different compressive stresses, restarting from the shearBack stage

    struct ShearStageData {
        Mdouble compressiveStress;
        Mdouble maxShearStress;
    };
    
    //Set the normal stress points for the shear stage
    std::array<ShearStageData,3> shearStageData;
    shearStageData[0].compressiveStress = 0.8 * compressiveStress;
    shearStageData[1].compressiveStress = 0.6 * compressiveStress;
    shearStageData[2].compressiveStress = 0.4 * compressiveStress;

    for (auto& data : shearStageData) {
        ShearStage dpm(speciesType, data.compressiveStress, shearRate,newName);
        dpm.restartFile.setFileType(FileType::NO_FILE);
        dpm.solve();
        data.maxShearStress = fabs(dpm.getMaxShearStress());
    }
    
    //store the final date to a string
    std::stringstream ss;
    ss << compressiveStress << " " << fabs(dpm.getPreShearShearStress()) << " ";
    for (const auto data : shearStageData) {
        ss << data.compressiveStress << " " << data.maxShearStress << " ";
    }

    std::cout << "Results:\n" << ss.str();
    helpers::writeToFile(newName+".txt", ss.str());
    
    //This part does the linear fit through the shear points and get the yield locus with slope_a and intercept_b
    //Then we calculate the sigma_1 and sigma_c to compute flow function FFc
    double slope_a, intercept_b;
    std::vector<double> Xdata{shearStageData[0].compressiveStress, shearStageData[1].compressiveStress, shearStageData[2].compressiveStress};
    std::vector<double> Ydata{shearStageData[0].maxShearStress, shearStageData[1].maxShearStress, shearStageData[2].maxShearStress};
    std::vector<double> linearReg = getLinearFit(Xdata,Ydata);
    slope_a = linearReg[0];
    intercept_b = linearReg[1];
    std::cout << "slope " << slope_a << " " << "Intercept " << intercept_b << std::endl;
    
    double sigma_ps = compressiveStress*1.02; //safety factor to avoid the pre-shear point on top of the yield locus
    double tau_ps =fabs(dpm.getPreShearShearStress());
    double sigma_center, radius_r1, radius_rc, FFc, sigma_c;
    
    double coeff_a, coeff_b, coeff_c;
    coeff_a = 1/(slope_a*slope_a+1);
    coeff_b = -2*(sigma_ps+slope_a*intercept_b/(slope_a*slope_a+1));
    coeff_c = (pow(sigma_ps,2)+pow(tau_ps,2))-pow(intercept_b,2)/(pow(slope_a,2)+1);
    double disc = (coeff_b*coeff_b)-(4*coeff_a*coeff_c); //The discriminant

    if(disc < 0) //If the discriminant is less than 0
    {
        std::cout << "Discriminant < 0. Therefore both roots are imaginary." << std::endl;
        sigma_center = NaN;
    }
    else if(disc == 0) //Use 'else if'. Not 'if' again.
    {
        std::cout << "Discriminant = 0. Therefore one repeated root." << std::endl;
        std::cout << "Result: " << -coeff_b/2.0*coeff_a;
        sigma_center = -coeff_b/2.0*coeff_a;
    }
    else //If the other conditions fail (i.e if discriminant > 0)
    {
        //Calculate the roots
        std::pair<double, double> result = solveQuad(coeff_a, coeff_b, coeff_c);

        //Output it
        std::cout << "Discriminant > 0. Therefore, two real roots." << std::endl;
        std::cout << "x1 = " << result.first << "\nx2 = " << result.second << std::endl;
        sigma_center = std::max(result.first,result.second);
    }
    
    radius_r1=fabs(slope_a*sigma_center+intercept_b)/sqrt(pow(slope_a,2)+1);
    radius_rc = (slope_a+sqrt(slope_a*slope_a+1))*intercept_b;
    sigma_c = 2*radius_rc;
    FFc = (sigma_center+radius_r1)/sigma_c;
    std::cout << "sigma_1 = " << sigma_center+radius_r1 << std::endl;
    std::cout << "sigma_c = " << sigma_c << std::endl;
    std::cout << "FFc = " << FFc << std::endl;
    
    std::stringstream ss1;
    ss1 << " " << FFc << " ";
    helpers::addToFile(newName+".txt", ss1.str());

    return 0;
}
