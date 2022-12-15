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
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <Boundaries/DeletionBoundary.h>
#include <Species/LinearPlasticViscoelasticFrictionReversibleAdhesiveSpecies.h>
#include "Boundaries/LeesEdwardsBoundary.h"
#include "Calibration.h"

/**
 * This simulation consists of two parts:
 * A)
 * Particles, with given volume fracion, are
 * placed in a cylindral column with a bottom, no lid.
 * Then all particles are relaxed.
 * B)
 * Remove the cylindrical wall and wait until the particles
 * are not moving anymore, then we extract the angle of repose.
 * Note:
 * user can set baseRadius, fillHeight, and give
 * input material name like APAPP with possiblity of changing
 * restitutionCoefficient, relativeCohesionStiffness,
 * slidingFriction and rollingFriction.
 **/
class GranuHeap : public Calibration
{
    bool fillStage = true;
    Mdouble angleOfRepose_ = 0;

public:
    Mdouble baseRadius;
    Mdouble fillHeight;

    GranuHeap(int argc, char* argv[]) : Calibration(argc,argv)
    {
        //set name
        setName("CalibrationHeap" + param);
        logger(INFO,"Name: %",getName());

        // define -test behavior
        if (helpers::readFromCommandLine(argc,argv,"-test")) {
            angleOfRepose_ = 45.0;
            writeResults();
            exit(EXIT_SUCCESS);
        }

        //set additional parameters
        const Mdouble initialVolumeFraction = 0.3;
        double d50 = psd.getVolumeDx(50);
        double dMean = psd.getNumberDx(50);
        double dMax = 2.0*psd.getMaxRadius();
        logger(INFO,"dMin % dMean % d50 %  dMax %", 2.0*psd.getMinRadius(), dMean, d50, dMax);
        baseRadius = std::max(std::max(15.0*dMean,3.0*d50),dMax);
        fillHeight = 3.0*baseRadius;
        logger(INFO,"baseRadius = % d50", round(baseRadius/psd.getNumberDx(50)*10)/10);
        setGravity({0,0,-9.8}); //set gravity in the system

        AxisymmetricIntersectionOfWalls base;
        base.setSpecies(frictionalWallSpecies);
        base.setAxis({0,0,1});
        base.addObject({-1,0,0},{baseRadius,0,0});
        base.addObject({0,0,-1},{0,0,0});
        wallHandler.copyAndAddObject(base);

        AxisymmetricIntersectionOfWalls casing;
        casing.setSpecies(frictionlessWallSpecies);
        casing.setAxis({0,0,1});
        casing.addObject({1,0,0},{baseRadius,0,0});
        wallHandler.copyAndAddObject(casing);
   
        setMax({baseRadius,baseRadius,fillHeight});
        setMin({-baseRadius,-baseRadius,0});
        setTimeMax(1000);
        //setSaveCount(100);

        //add particles
        CubeInsertionBoundary insertion;
        SphericalParticle particle(speciesHandler.getObject(0));
        insertion.setPSD(psd);
        double cylinderVolume = getTotalVolume()*constants::pi/4.0;
        insertion.setInitialVolume(initialVolumeFraction*cylinderVolume);
        Vec3D max = getMax();
        do {
            insertion.set(&particle, 10000, getMin(), max, Vec3D(0, 0, 0), Vec3D(0, 0, 0), 0, 0);
            insertion.checkBoundaryBeforeTimeStep(this);
            max.Z *= 1.05;
            logger(INFO,"Inserted % particles, fill fraction %",particleHandler.getSize(),particleHandler.getVolume()/cylinderVolume);
        } while (insertion.getVolumeOfParticlesInserted()<initialVolumeFraction*cylinderVolume);
        writePSDToFile();
    }

    void printTime() const override {
        logger(INFO, "nt % COM % Ene % N %", getNumberOfTimeSteps(), getCentreOfMass().Z, getKineticEnergy()/getElasticEnergy(), particleHandler.getNumberOfObjects());
    }
    
    void actionsAfterTimeStep() override
    {
        //check only every 1000th time
        static unsigned count = 0;
        if (++count>999) count=0; else return;
        static bool timeSwitch = true;
        double fillTime;

        if (fillStage) {
            if (getKineticEnergy()<1e-2*getElasticEnergy())
            {
                fillStage = false;
                fillTime = getTime();
                logger(INFO,"remove outer cylinder");
                wallHandler.removeLastObject();
                logger(INFO,"remove particles below z=0");
                DeletionBoundary boundary;
                boundary.set({0,0,-1},0);
                boundaryHandler.copyAndAddObject(boundary);
            }
        } else if (getTime()>2.0*fillTime) {
            //why is it measured twice?
            if (getKineticEnergy()<1e-2*getElasticEnergy() && timeSwitch) {
                Mdouble COM = getCentreOfMass().Z;
                Mdouble height = 4.0* COM;
                angleOfRepose_ = atan(height/baseRadius)*180.0/pi;
                logger(INFO,"AoR=%, based on conical heap shape",angleOfRepose_);
                setTimeMax(getTime()+0.05);
                timeSwitch = false;
            } else if (getKineticEnergy()<1e-2*getElasticEnergy() && !timeSwitch) {
                Mdouble COM = getCentreOfMass().Z;
                Mdouble height = 4.0* COM;
                // angle of repose based on the shape of a pyramid
                angleOfRepose_ = atan(height/baseRadius)*180.0/pi;
                logger(INFO,"AoR=%, based on conical heap shape",angleOfRepose_);
            }
        }
    }

    void writeResults() {
        std::stringstream ss;
        ss << angleOfRepose_;
        std::cout << "Results:\n" << ss.str();
        helpers::writeToFile(getName()+".txt", ss.str());
    }
};

int main(int argc, char* argv[])
{
    GranuHeap dpm(argc,argv);
    dpm.solve();
    dpm.writeResults();
    return 0;
}