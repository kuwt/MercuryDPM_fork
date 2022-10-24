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
#include "Mercury3D.h"
#include "Walls/IntersectionOfWalls.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Species/LinearViscoelasticFrictionLiquidMigrationLSSpecies.h"
#include "Particles/LiquidFilmParticle.h"

/**
 * This code simulates a 3D hour glass, i.e. a cylindrical domain with a neck in the middle.
 * Particles are inserted into the upper half of the domain and flow into the lower half due to gravity.
 * If the friction is high enough, arching is observed at the neck, where particles interlock to obstruct the flow.
 *
 * This code illustrates the effect of friction in DPM simulations.
 * It also illustrates how to use AxisymmetricIntersectionOfWalls to set up axisymmetrical shapes, i.e. the outer cylinder and the neck.
 */
class AoRTest : public Mercury3D
{
public:

    void setupInitialConditions() override {
        Vec3D mid = {
            (getXMin() + getXMax()) / 2.0,
            (getYMin() + getYMax()) / 2.0,
            (getZMin() + getZMax()) / 2.0};
        const Mdouble halfWidth = (getXMax() - getXMin()) / 2.0;

        //species information
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionLiquidMigrationLSSpecies());
        species->setDensity(2650);
        species->setStiffness(1500);
        species->setDissipation(0.002);
        const Mdouble mass = species->getMassFromRadius(radius);
        species->setRestitutionCoefficient(0.2, mass);//need information here
        species->setSlidingStiffness(0.4*species->getStiffness());
        species->setSlidingDissipation(0.4*species->getDissipation());
        species->setSlidingFrictionCoefficient(0.4);
        species->setSlidingFrictionCoefficientStatic(0.4);
        species->setRollingStiffness(1500);
        species->setRollingDissipation(0.002);
        species->setRollingFrictionCoefficient(0.2);

        species->setLiquidBridgeVolumeMax(mass * 0.1 * 2/1000.0);
        species->setLiquidBridgeVolumeMin(0.0);
        species->setDistributionCoefficient(0.5);
        species->setSurfaceTension(0.07);
        species->setContactAngle(constants::pi/18);
        species->setViscosity(0.00085);
        //tangential (rolling) torques
//    species.setRollingFrictionCoefficient(0.2);
//    species.setRollingStiffness(1.2e4);
//    species.setRollingDissipation(6.3e-2);
        //normal (torsion/spin) torques
//    species.setTorsionFrictionCoefficient(0.1);
//    species.setTorsionStiffness(1.2e4);
//    species.setSlidingDissipation(6.3e-2);
//        InfiniteWall w2;
//        w2.setSpecies(speciesHandler.getObject(0));
//        w2.set(Vec3D(1, 0, 0), Vec3D(0.0075, 0, 0));
//        wallHandler.copyAndAddObject(w2);
//        w2.set(Vec3D(-1, 0, 0), Vec3D(-0.0025, 0, 0));
//        wallHandler.copyAndAddObject(w2);
//        w2.set(Vec3D(0, -1, 0), Vec3D( 0, -0.0025, 0));
//        wallHandler.copyAndAddObject(w2);
//        w2.set(Vec3D(0, 1, 0), Vec3D( 0, 0.0075, 0));
//        wallHandler.copyAndAddObject(w2);



        //! [CST:outer]
        AxisymmetricIntersectionOfWalls w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.setPosition(Vec3D(mid.X, mid.Y, 0));
        //w1.setAxis(Vec3D(0, 0, 1));
        w0.setOrientation(Vec3D(0,0,1));
        w0.addObject(Vec3D(1, 0, 0), Vec3D((getXMax() - getXMin()) / 2.0, 0, mid.Z/2.0));
        w0.addObject(Vec3D(0,0,1), Vec3D((getXMax() - getXMin()) / 2.0,0, mid.Z/2.0));
        wallHandler.copyAndAddObject(w0);
        //! [CST:outer]

        //top wall
//        InfiniteWall w0;
//        w0.setSpecies(speciesHandler.getObject(0));
//        w0.set(Vec3D(0, 0, 1), Vec3D(mid.X, mid.Y, getZMax()));
//        wallHandler.copyAndAddObject(w0);

        //bottom wall
        InfiniteWall w1;
        w1.setSpecies(speciesHandler.getObject(0));
        w1.set(Vec3D(0, 0, -1), Vec3D(mid.X, mid.Y, mid.Z/2.0));
        wallHandler.copyAndAddObject(w1);


        LiquidFilmParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        for (Mdouble z = mid.Z/2.0;
            particleHandler.getNumberOfObjects() <= N;
            z += 2.0 * radius)
        {
            for (Mdouble r = halfWidth - radius; r > 0; r -= 1.999 * radius)
            {
                for (Mdouble c = 2.0 * radius; c <= 2 * constants::pi * r; c += 2.0 * radius)
                {
                    p0.setRadius(radius);
                    p0.setPosition(Vec3D(mid.X + r * mathsFunc::sin(c / r),
                                         mid.Y + r * mathsFunc::cos(c / r),
                                         z + p0.getRadius()));
                    const Mdouble vz = random.getRandomNumber(-0.03, 0);
                    const Mdouble vx = random.getRandomNumber(-0.03, 0.03);
                    p0.setVelocity(Vec3D(vx,0.0,vz));

                    Mdouble mass_lfp = 4.0/3.0 * constants::pi * mathsFunc::cubic(p0.getRadius()) * species->getDensity();
                    p0.setLiquidVolume(0.1 * mass_lfp/1000.0);
                    particleHandler.copyAndAddObject(p0);

                }
            }
        }
    }

    //Initially, a wall is inserted in the neck of the hourglass to prevent particles flowing through.
    //This wall is moved to form the base of the hourglass at time 0.9
    void actionsAfterTimeStep() override {
        if (getTime() < 0.05 && getTime() + getTimeStep() > 0.05)
        {
            logger(INFO, "Move the bottom wall");
            dynamic_cast<InfiniteWall*>(wallHandler.getLastObject())->set(Vec3D(0, 0, -1), Vec3D(0, 0, getZMin()));
//            wallHandler.removeLastObject();
//            InfiniteWall w3;
//            w3.setSpecies(speciesHandler.getObject(0));
//            w3.set(Vec3D(0, 0, 1), Vec3D(0,0,0));
//            wallHandler.copyAndAddObject(w3);
        }
    }

    Mdouble radius;
    unsigned int N;
};

int main(int argc, char *argv[])
{
    //logger(INFO, "Hourglass Simulation");
    // note: the scaling of the hourglass is based on Stefan Luding's demo code MDCLR/DEMO/W7

    //all parameters should be set in the main file
    //here, we use SI units (generally, other, scaled units are possible)

    //create an instance of the class and name it
    AoRTest AoR;
    AoR.setName("AngleofReposeTest_rollingfriction_02");


    //these parameters are needed in setupInitialConditions()

    //specify geometry
    //specify dimensions of the hourglass
    const Mdouble width = 5e-3; // 0.5cm
    const Mdouble height = 1.5e-2; // 1cm
    //set domain accordingly (domain boundaries are not walls!)
    AoR.setXMin(0.0);
    AoR.setXMax(width);
    AoR.setYMin(0.0);
    AoR.setYMax(width);
    AoR.setZMin(0.0);
    AoR.setZMax(height);

    //specify how big the wedge of the contraction should be
//    const Mdouble contractionWidth = 1.25e-2; //1.25cm
//    const Mdouble contractionHeight = 2.5e-2; //2.5cm
//    AoR.contractionWidth = contractionWidth;
//    AoR.contractionHeight = contractionHeight;

    //specify particle properties

    //these parameters are needed in setupInitialConditions()
    AoR.radius = 0.00025; // 0.25mm

    //specify body forces
    AoR.setGravity(Vec3D(0.0, 0.0, -9.8));

    //test normal forces
//    const Mdouble minParticleMass = species.getDensity() * 4.0 / 3.0 * constants::pi * mathsFunc::cubic(AoR.minParticleRadius);
//    //Calculates collision time for two copies of a particle of given dissipation_, k, effective mass
//    logger(INFO, "minParticleMass = %", minParticleMass);
//    //Calculates collision time for two copies of a particle of given dissipation_, k, effective mass
//    const Mdouble tc = species.getCollisionTime(minParticleMass);
//    logger(INFO, "tc  = %", tc);
//    //Calculates restitution coefficient for two copies of given dissipation_, k, effective mass
//    const Mdouble r = species.getRestitutionCoefficient(minParticleMass);
//    logger(INFO, "restitution coefficient = %", r);

    //set other simulation parameters
    AoR.setTimeStep(1e-6);
    AoR.setTimeMax(0.2);
    AoR.setSaveCount(1000);


    AoR.setXBallsAdditionalArguments("-v0 -solidf");
    AoR.N = 1500; //number of particles
    //uncomment next two line 2 to create paraview files
    AoR.setParticlesWriteVTK(true);
    AoR.setWallsWriteVTK(FileType::ONE_FILE);
    AoR.setInteractionsWriteVTK(true);
    AoR.setFileType(FileType::ONE_FILE);

    AoR.solve(argc, argv);
    return 0;
}

