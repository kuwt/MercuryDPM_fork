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
#include <Species/LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies.h>
#include <Mercury3D.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Boundaries/StressStrainControlBoundary.h>
#include <Boundaries/CubeInsertionBoundary.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <Boundaries/DeletionBoundary.h>
#include "Boundaries/LeesEdwardsBoundary.h"
#include "addSpecies.h"
#include <CG/TimeSmoothedCG.h>
#include "CG/TimeAveragedCG.h"

/**
 * This simulation consists of two parts:
 * A)
 * Particles, with given volume fracion, are
 * placed in a Rotating Drum. Then relaxed until
 * all the particles are not moving
 * B)
 * Start rotating the drum with a given Froude number Fr
 * until steady state and then extract the dynamic angle of repose.
 * Note:
 * user can set drumFillFraction, froudeNumber, and give
 * input material name like APAPP with possiblity of changing
 * restitutionCoefficient, relativeCohesionStiffness,
 * slidingFriction and rollingFriction.
 **/
class GranuDrum : public Mercury3D
{
public:

    /**
     *
     * @param stressGoal the stress tensor that the material should be subjected to
     * @param strainRate the applied strain rate
     * @param gainFactor the sensitivity by which volume is modified to achieve the stress, or strain, goal
     * @param isStrainRateControlled if this is switched off, then stress is applied and strain is controlled
     */
    GranuDrum(SpeciesType type)
    {
        //add species
        auto psd = getPSD(type);
        logger.assert(hasPSD(type),"Material type % not allowed",type);
        auto species = addSingleSpecies(type, getMedianParticleRadius(psd), getMinParticleRadius(psd), *this, true);
        logger(INFO, "%", *species);
    
        drumRadius = 30.0*getMedianParticleRadius(psd);
        drumDepth = 0.3*drumRadius;
        rotRate = froudeNumber*sqrt(9.81/drumRadius);

        setName("GranuDrum");
        setGravity({0,0,-9.8}); //set gravity in the system

        // Adding Geometry
		wallHandler.clear();
		AxisymmetricIntersectionOfWalls cylinder;
		cylinder.setSpecies(species);
		cylinder.setAxis({0.0,1.0,0.0});
		cylinder.addObject({1.0,0.0,0.0},{drumRadius,0.0,0.0});
		wallHandler.copyAndAddObject(cylinder);

		InfiniteWall w0;
		w0.setSpecies(species);
		w0.set(Vec3D(0.,-1.,0.),Vec3D(0.0,-0.5*drumDepth,0.0));
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D(0.,1.,0.),Vec3D(0.0,0.5*drumDepth,0.0));
		wallHandler.copyAndAddObject(w0);
        logger(INFO,"Walls have been placed");

        // Add Species
        setMax({drumRadius,0.5*drumDepth,drumRadius});
        setMin({-drumRadius,-0.5*drumDepth,-drumRadius});
        setTimeMax(30);
        setSaveCount(2000);

        ///Commented is my regular way of inputting particles
        //Mdouble currFillVolume = 0.0;
        //int failCounter;
        //BaseParticle p0;
        //Mdouble R,C;
        //Mdouble drumFillVolume = drumFillFraction*drumDepth*constants::pi*pow(drumRadius,2);

		//while (currFillVolume < drumFillVolume)
		//{
        //    failCounter = 0;
        //    do
        //    {
        //        p0.setSpecies(speciesS1);
        //        p0.setRadius(radiusS1);
        //        R = random.getRandomNumber(2*radiusS1,drumRadius-2.0*p0.getRadius()); // Find a random radial position within the cylinder
        //        C = random.getRandomNumber(0.0*constants::pi,2.0*constants::pi); // Find a random angular position within the cylinder
        //        // Set the position in the cylinder with a random z-coordinate
        //        p0.setPosition(Vec3D(R*cos(C),random.getRandomNumber(getYMin()+p0.getRadius(),getYMax()-p0.getRadius()),R*sin(C)));
        //        p0.setVelocity(Vec3D(0.0,0.0,0.0)); // Set initial particle velocity
        //        failCounter++; // Adding 1 to failcounter
        //            if (failCounter==1000)
        //            {
        //                break; // If we failed 1000 times to find a non-contact position we break the placement loop (makes sure that no infinite loops are created in combination with line 130)
        //            }
        //    } while(!checkParticleForInteraction(p0));
        //    particleHandler.copyAndAddObject(p0);
        //    numS1Inserted++;
        //    hGridRebuild();
		//}

      /// Ive got problems here as I cant get the desired volume fill fraction set by drumFillFraction
        CubeInsertionBoundary insertion;
        BaseParticle particle(speciesHandler.getObject(0));
        insertion.set(&particle, 10, getMin(), getMax(), Vec3D(0, 0, 0), Vec3D(0, 0, 0), 0, 0);
        insertion.setPSD(getPSD(type));
        insertion.setInitialVolume(1*drumFillFraction*constants::pi*pow(drumRadius,2)*drumDepth);
        //insertion.setInitialVolume(100); // No clue what this line does, and couldnt find documentation
        //insertion.checkBoundaryBeforeTimeStep(this);

        currFillFrac = getTotalMass() / getDensity(type) / (constants::pi*pow(drumRadius,2)*drumDepth);
        boundaryHandler.copyAndAddObject(insertion);

        logger(INFO,"Inserted % particles",particleHandler.getSize());
        logger(INFO,"Volume fill fraction %",currFillFrac);
        particleDensity = getDensity(type);
        //BoundaryHandler.removeLastObject();
        //logger(INFO,"Removed insertion boundary");
    }
    
//    void computeExternalForces (BaseParticle * CI) override
//    {
//        DPMBase::computeExternalForces (CI);
//        if (fillStage && !CI->isFixed())
//        {
//            // Background damping
//            CI->addForce(-dampingCoefficient * CI->getVelocity());
//            //logger(INFO,"damping coefficient %",dampingCoefficient);
//        }
//    }

    void printTime() const override {
        logger(INFO, "t % N % COM %", getTime(), particleHandler.getNumberOfObjects(), getCentreOfMass().Z);
    }

    void actionsAfterTimeStep() override
    {
        //check only every 1000th time
        static unsigned count = 0;
        if (++count>999) count=0; else return;


        if (placeStage) {
            currFillFrac = getTotalMass() / particleDensity / (constants::pi*pow(drumRadius,2)*drumDepth); // Inserted volume divided by drum volume
            logger(INFO,"currFillFrac %, totalMass %",currFillFrac,getTotalMass());
            if (currFillFrac > drumFillFraction && placeStage == true) {
                logger(INFO,"Inserted % particles",particleHandler.getSize());
                logger(INFO,"Volume fill fraction %",currFillFrac);
                boundaryHandler.removeLastObject();
                placeStage = false;
                fillStage = true;
            }
        }

        if (fillStage) {
            if (getKineticEnergy()<1e-5*getElasticEnergy())
            {
                fillStage = false;
                logger(INFO,"Particles have settled, starting rotation");
                logger(INFO,"rpm %",rotRate);

                wallHandler.getObject(0)->setAngularVelocity(Vec3D(0.0,rotRate,0.0));
				wallHandler.getObject(1)->setAngularVelocity(Vec3D(0.0,rotRate,0.0));
				wallHandler.getObject(2)->setAngularVelocity(Vec3D(0.0,rotRate,0.0));
				
            }
        }

    }

    Mdouble currFillFrac;
    Mdouble particleDensity;
    const Mdouble drumFillFraction = 0.25;
    Mdouble drumRadius;
    Mdouble drumDepth;
    bool placeStage = true;
    bool fillStage = false;

	const double froudeNumber = 0.06;
	Mdouble rotRate;// this is in rad/s, so no need to convert again
    Mdouble dampingCoefficient;
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

int main(int argc, char* argv[])
{
    //Reading in the five parameters
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
    
    std::stringstream restartName;
    //restartName << std::fixed << std::setprecision(5) << "DrumFE" << "_e" << restitutionCoefficient << "_rC" << relativeCohesionStiffness << "_mus" << slidingFriction << "_mur" << rollingFriction;
    restartName << std::setprecision(6) << std::scientific << "DrumFE" << "_e" << restitutionCoefficient << "_rC" << relativeCohesionStiffness << "_mus" << slidingFriction << "_mur" << rollingFriction;
    const auto newName = restartName.str();
    std::cout << newName << std::endl;
    
    GranuDrum dpm(speciesType);
    modifySpecies(dpm.speciesHandler.getLastObject(), restitutionCoefficient, relativeCohesionStiffness, slidingFriction, rollingFriction, speciesType);
//    auto species = dynamic_cast<LinearPlasticViscoelasticFrictionSpecies*>(dpm.speciesHandler.getLastObject());
//    dpm.dampingCoefficient = 0.1*species->getDissipation();
    dpm.setName(newName);
    logger(INFO,"Species: %", *dpm.speciesHandler.getLastObject());
    dpm.removeOldFiles();
    dpm.setParticlesWriteVTK(false);
    dpm.setWallsWriteVTK(false);
    dpm.dataFile.setFileType(FileType::ONE_FILE);
    dpm.eneFile.setFileType(FileType::NO_FILE);
    dpm.fStatFile.setFileType(FileType::NO_FILE);
    dpm.restartFile.setFileType(FileType::ONE_FILE);
    dpm.restartFile.setSaveCount(100000000);

    {
        TimeAveragedCG<CGCoordinates::XZ, CGFunctions::Heaviside> c;
        c.setN({50, 1, 50});
        c.setWidth(1.5e-3);
        c.setTimeMin(10);
        c.statFile.setSaveCount(25);
        c.statFile.setName(dpm.getName() + "XZ.TA.stat");
        dpm.cgHandler.copyAndAddObject(c);
    }

    dpm.solve();
    
    static const int max_line = 65536;
    std::ifstream statFile(dpm.getName() + "XZ.TA.stat");
    statFile.ignore(max_line,'\n');
    statFile.ignore(max_line,'\n');
    std::vector<double> col_x, col_z;
    double aa, bb, cc, dd, ee;
    
    while (statFile >> aa >> bb >> cc >> dd >> ee){
        if(relativeCohesionStiffness <= 0.1 && dd>0.2 && dd<0.35 && (bb*bb+cc*cc)<0.0002){
            col_x.push_back(bb);
            col_z.push_back(cc);
        } else if (relativeCohesionStiffness > 0.1 && relativeCohesionStiffness <= 0.5 && dd>0.15 && dd<0.3 && (bb*bb+cc*cc)<0.00016) {
            col_x.push_back(bb);
            col_z.push_back(cc);
        } else if (relativeCohesionStiffness > 0.5 && dd>0.12 && dd<0.16 && (bb*bb+cc*cc)<0.0002) {
            col_x.push_back(bb);
            col_z.push_back(cc);
        }
        statFile.ignore(max_line,'\n');
    }
    
    double slope_a, intercept_b, angle_drum;

    std::vector<double> linearReg = getLinearFit(col_x,col_z);
    slope_a = linearReg[0];
    intercept_b = linearReg[1];
    angle_drum = atan(fabs(slope_a))*180/constants::pi;
    std::cout << "slope " << slope_a << " " << "Intercept " << intercept_b << std::endl;
    std::cout << "Angle " << angle_drum << std::endl;
    
    std::stringstream ss;
    ss << angle_drum;
    
    std::cout << "Results:\n" << ss.str();
    helpers::writeToFile(newName+".txt", ss.str());
    
    return 0;
}
