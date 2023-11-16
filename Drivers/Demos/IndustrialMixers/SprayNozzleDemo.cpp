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
#include <Species/HertzianViscoelasticMindlinLiquidMigrationWilletSpecies.h>
#include <Boundaries/CubeInsertionBoundary.h>
#include <Boundaries/DropletBoundary.h>
#include <Particles/LiquidFilmParticle.h>

/**
 * A new demo for the insertion of liquid droplets using nozzle.
 * Note that you can switch on or off the nozzle by changing the value of NozzleDemo::addNozzle.
 *
 * A video of this demo can be found here <https://youtu.be/Ffu063_qLz0>.
 */

class NozzleDemo : public Mercury3D {

    // Group id of the rotating geometry
    unsigned rotatingWallID = 0;

    // Pointers to the insertion boundaries
    CubeInsertionBoundary* insertionBoundary;
    DropletBoundary* dropletBoundary = nullptr;

    // Whether the nozzle gets added to the mixer
    bool addNozzle = true;


    void setupInitialConditions () override {
        // Set name of output files
        setName("NozzleDemo");
        setGravity(Vec3D(0,-9.8,0));
        // Remove existing vtk files
        removeOldFiles();
        // Set domain size
        setMin({-0.08,-0.08,0.0});
        setMax({0.08,0.08,0.085});
        // Set time step, final time and how often to output
        setTimeMax(10.0);
        setSaveCount(300);
        // Output files: wall-vtu and particle-vtu files
        wallHandler.setWriteVTK(true);
        setParticlesWriteVTK(true);
        //Contact law and density
        HertzianViscoelasticMindlinLiquidMigrationWilletSpecies species;
        // \todo check effective elastic modulus and shear modulus
        Mdouble elasticModulus = 1.325e6;
        Mdouble shearModulus = 4e5;
        Mdouble poisson = elasticModulus/(2*shearModulus)-1;
        Mdouble effectiveElasticModulus=elasticModulus/(2*poisson);
        Mdouble effectiveShearModulus=effectiveElasticModulus / 2 * (1 + poisson);
        species.setEffectiveElasticModulusAndRestitutionCoefficient(effectiveElasticModulus, 0.5);
        species.setEffectiveShearModulus(effectiveShearModulus);
        species.setSlidingFrictionCoefficient(0.61);
        species.setLiquidBridgeVolumeMax(1e-8);
        species.setDensity(900);
        speciesHandler.copyAndAddObject(species);
        // Introduce a rotating wall
        Mdouble wallScaleFactor = 1e-3; // Scale used in the stl file (mm)
        Vec3D shift = {0,0,0};
        Vec3D velocity = {0,0,0};
        rotatingWallID = wallHandler.readTriangleWall("RotatingDrum.stl",speciesHandler.getObject(0), wallScaleFactor,shift,velocity,Vec3D(0,0,0));
        // Add particles
        CubeInsertionBoundary boundary;
        LiquidFilmParticle particle;
        Mdouble radius = 400e-5;
        double insertionVolume = 0.2/species.getDensity();
        logger(INFO,"Insertion volume: %", insertionVolume);

        particle.setSpecies(speciesHandler.getObject(0));
        particle.setRadius(radius);
        boundary.set(particle, 1000, Vec3D(0.018,0.003,0), Vec3D(0.053,0.05,0.085),Vec3D(0,0,0),Vec3D(0,0,0));
        boundary.setInitialVolume(insertionVolume);
        insertionBoundary = boundaryHandler.copyAndAddObject(boundary);

        if (addNozzle) {
            // Volumetric flow rate of nozzle (m^3/s)
            double flowRate = 1e-5; //10 ml/s  =1e-5
            // Radius of droplet (m)
            double dropletRadius = 1e-3;
            // Volume of droplet (m^3)
            double dropletVolume = std::pow(2.0 * dropletRadius, 3) * constants::pi / 6.0;
            // Angle in radians where the nozzle is located
            double nozzleAngle = 25. * constants::pi / 180.;
            double cosNozzleAngle = cos(nozzleAngle);
            double sinNozzleAngle = sin(nozzleAngle);
            // Angle in radians of the flat-fan, i.e. how wide the spray is distributed in z-direction
            double sprayAngle = 25. * constants::pi / 180.;
            double sinSprayAngle = sin(sprayAngle);
            // Speed of each droplet at nozzle exit (m/s)
            double dropletSpeed = 1;
            // Radial distance of nozzle from rotation axis
            double cylinderRadius = 0.07;
            // Add a DropletBoundary that acts as a flat-fan spray nozzle
            DropletBoundary d;
            d.setGenerateDroplets(
                    [this, flowRate, dropletRadius, dropletVolume, dropletSpeed, sinNozzleAngle, cosNozzleAngle, sinSprayAngle, cylinderRadius](
                            DropletBoundary &d) {
                        static double accumulatedDropletVolume = 0;
                        while (accumulatedDropletVolume < flowRate * getTime()) {
                            Vec3D position = Vec3D(-sinNozzleAngle * cylinderRadius, cosNozzleAngle * cylinderRadius,
                                                   0.5 * (getZMax() + getZMin()));
                            Vec3D velocity = Vec3D(cosNozzleAngle, -sinNozzleAngle, random.getRandomNumber(-1, 1)) * sinSprayAngle * dropletSpeed;
                            d.droplets_.emplace_back(DropletBoundary::Droplet{position, velocity, dropletRadius});
                            accumulatedDropletVolume += dropletVolume;
                        }
                    });
            logger(INFO, "Insert droplet boundary");
            dropletBoundary = boundaryHandler.copyAndAddObject(d);
            // Turn on vtk output for boundaries
            boundaryHandler.setWriteVTK(true);
        }

        // \todo check
        Mdouble dt = 0.2*helpers::getRayleighTime(radius, shearModulus, poisson, species.getDensity());
        setTimeStep(dt); //Find the rayleigh time step
        logger(INFO,"Time step: %", dt);
    }

    void printTime() const override {
        logger(INFO,"t=%, tmax=%, N=%, ND=%", getTime(), getTimeMax(), particleHandler.getSize(),dropletBoundary?dropletBoundary->droplets_.size():0);
    }

    // Start the rotation after the particles have been inserted
    void actionsAfterTimeStep() override {
        static bool rotationStarted = false;
        Vec3D angularVelocity = Vec3D(0,0,4.0/3.0*constants::pi);
        if (!rotationStarted and insertionBoundary->getInitialVolume()<=insertionBoundary->getVolumeOfParticlesInserted()) {
            rotationStarted = true;
            logger(INFO,"Starting rotation");
            for (const auto wall : wallHandler) {
                if (wall->getGroupId()==rotatingWallID) {
                    wall->setAngularVelocity(angularVelocity);
                }
            }
        }

    }
};

int main(int argc, char* argv[])
{

    // Problem setup
    NozzleDemo problem;
    // Start solving in time
    problem.solve();
    logger(INFO,"Load %Wall_*.vtu in paraview to see the wall geometry",problem.getName());

    return 0;
}
