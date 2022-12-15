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
#include <Boundaries/PeriodicBoundary.h>
#include <Boundaries/CubeInsertionBoundary.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <Boundaries/DeletionBoundary.h>
#include "Calibration.h"

/**
 * This simulation consists of two parts:
 * A)
 * Particles, with given volume fraction, are
 * placed in a Rotating Drum. Then relaxed until
 * all the particles are not moving
 * B)
 * Start rotating the drum with a given Froude number Fr
 * until steady state and then extract the dynamic angle of repose.
 * Note:
 * user can set fillFraction, froude, restitutionCoefficient, relativeCohesionStiffness,
 * slidingFriction and rollingFriction.
 **/
class GranuDrum : public Calibration
{
    Mdouble drumRadius;
    Mdouble drumDepth;
    bool placeStage = true;
    Mdouble rotRate;// this is in rad/s, so no need to convert again
    double angle_drum;

public:

    /*
     * Calibration(argc,argv) already sets the species, type and psd.
     */
    GranuDrum(int argc, char* argv[]) : Calibration(argc,argv)
    {
        const double froude = helpers::readFromCommandLine(argc, argv, "-froude", 0.06);
        const Mdouble fillFraction = helpers::readFromCommandLine(argc, argv, "-fillFraction", 0.25);
        const double revolutions = helpers::readFromCommandLine(argc, argv, "-revolutions", 2.0);
        const double relDrumRadius = helpers::readFromCommandLine(argc, argv, "-relDrumRadius", 15.0);
        const double relDrumWidth = helpers::readFromCommandLine(argc, argv, "-relDrumWidth", 4.5);

        //set name
        setName("CalibrationDrum" + param);
        logger(INFO,"Name: %",getName());
        removeOldFiles();

        // define -test behavior
        if (helpers::readFromCommandLine(argc,argv,"-test")) {
            angle_drum = 45.0;
            writeResults();
            exit(EXIT_SUCCESS);
        }


        //set drum radius, depth
        drumRadius = std::max(relDrumRadius*psd.getVolumeDx(50),4.0*psd.getMaxRadius());
        if (drumRadius>4.0*psd.getMaxRadius()) {
            logger(INFO, "drumRadius = %*d50", relDrumRadius);
        } else {
            logger(INFO, "drumRadius = 2*d100");
        }
        drumDepth = std::max(relDrumWidth*psd.getVolumeDx(50),2.4*psd.getMaxRadius());
        if (drumDepth>2.4*psd.getMaxRadius()) {
            logger(INFO, "drumDepth = %*d50", relDrumWidth);
        } else {
            logger(INFO, "drumDepth = 1.2*d100");
        }
        // domain
        setMax({drumRadius,0.5*drumDepth,drumRadius});
        setMin({-drumRadius,-0.5*drumDepth,-drumRadius});

        //set gravity and rotation rate
        setGravity({0,0,-9.8}); //set gravity in the system
        rotRate = froude*sqrt(getGravity().getLength()/drumRadius);
        //rotate twice
        setTimeMax(2.0*pi/rotRate*revolutions);
        logger(INFO,"Simulating % revolutions at Froude number %",revolutions, froude);

        // adding cylindrical wall
		AxisymmetricIntersectionOfWalls cylinder;
		cylinder.setSpecies(frictionalWallSpecies);
		cylinder.setAxis({0.0,1.0,0.0});
		cylinder.addObject({1.0,0.0,0.0},{drumRadius,0.0,0.0});
		cylinder.setAngularVelocity(Vec3D(0.0,rotRate,0.0));
		wallHandler.copyAndAddObject(cylinder);

		//adding side walls
		InfiniteWall w0;
		w0.setSpecies(frictionlessWallSpecies);
        w0.setAngularVelocity(Vec3D(0.0,rotRate,0.0));
		w0.set(Vec3D(0.,-1.,0.),Vec3D(0.0,-0.5*drumDepth,0.0));
		wallHandler.copyAndAddObject(w0);
		w0.set(Vec3D(0.,1.,0.),Vec3D(0.0,0.5*drumDepth,0.0));
		wallHandler.copyAndAddObject(w0);
        logger(INFO,"Walls have been placed");

        /// Ive got problems here as I cant get the desired volume fill fraction set by fillFraction
        CubeInsertionBoundary insertion;
        SphericalParticle particle(speciesHandler.getObject(0));
        insertion.set(&particle, 100, getMin(), getMax(), Vec3D(0, 0, 0), Vec3D(0, 0, 0), 0, 0);
        insertion.setPSD(psd);
        insertion.setInitialVolume(1*fillFraction*constants::pi*pow(drumRadius,2)*drumDepth);
        insertion.setCheckParticleForInteraction(false);
        insertion.setHandler(&boundaryHandler);
        insertion.checkBoundaryBeforeTimeStep(this);

        if (helpers::readFromCommandLine(argc,argv,"-output")) {
            setSaveCount(static_cast<unsigned>(getTimeMax()/getTimeStep()/100.0));
        }
    }

    void printTime() const override {
        logger(INFO, "t % N % AoR % del %", getTime()/getTimeMax(), particleHandler.getNumberOfObjects(), atan(getCentreOfMass().X/getCentreOfMass().Z)*180.0/constants::pi, interactionHandler.getMeanOverlap()/particleHandler.getMeanRadius());
    }

    void writeResults() {
        // center of mass; both x and z-values  should be negative
        Vec3D com = getCentreOfMass();
        angle_drum = atan(com.X/com.Z)*180/constants::pi;
        logger(INFO,"Angle % COM (%, %)",angle_drum,com.X,com.Z);
        double r = sqrt(com.X * com.X + com.Z * com.Z)/drumRadius;
        if (r<0.2) {
            logger(WARN,"Center of mass is too close to the center, |COM|/R=%, setting angle to 90",r);
            angle_drum=90;
        }
        if (com.X>0 || com.Z>0) {
            logger(WARN,"Center of mass (% %) is not in fourth quadrant, setting angle to 90",com.X,com.Z);
            angle_drum=90;
        }
        helpers::writeToFile(getName()+".txt", helpers::to_string(angle_drum));
    }
};

int main(int argc, char* argv[])
{
    GranuDrum dpm(argc,argv);
    //dpm.stats();
    dpm.solve();
    //dpm.restart();
    dpm.writeResults();
    return 0;
}
