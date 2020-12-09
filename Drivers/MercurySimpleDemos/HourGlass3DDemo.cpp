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
#include "Species/LinearViscoelasticFrictionSpecies.h"

/**
 * This code simulates a 3D hour glass, i.e. a cylindrical domain with a neck in the middle.
 * Particles are inserted into the upper half of the domain and flow into the lower half due to gravity.
 * If the friction is high enough, arching is observed at the neck, where particles interlock to obstruct the flow.
 *
 * This code illustrates the effect of friction in DPM simulations.
 * It also illustrates how to use AxisymmetricIntersectionOfWalls to set up axisymmetrical shapes, i.e. the outer cylinder and the neck.
 */
class HourGlass : public Mercury3D
{
public:

    void setupInitialConditions() override {
        Vec3D mid = {
            (getXMin() + getXMax()) / 2.0,
            (getYMin() + getYMax()) / 2.0,
            (getZMin() + getZMax()) / 2.0};
        const Mdouble halfWidth = (getXMax() - getXMin()) / 2.0;

        //! [CST:outer]
        AxisymmetricIntersectionOfWalls w1;
        w1.setSpecies(speciesHandler.getObject(0));
        w1.setPosition(Vec3D(mid.X, mid.Y, 0));
        w1.setAxis(Vec3D(0, 0, 1));
        w1.addObject(Vec3D(1, 0, 0), Vec3D((getXMax() - getXMin()) / 2.0, 0, 0));
        wallHandler.copyAndAddObject(w1);
        //! [CST:outer]

        //! [CST:neck]
        AxisymmetricIntersectionOfWalls w2;
        w2.setSpecies(speciesHandler.getObject(0));
        w2.setPosition(Vec3D(mid.X, mid.Y, 0));
        w2.setAxis(Vec3D(0, 0, 1));
        std::vector<Vec3D> points(3);
        //define the neck as a prism through  corners of your prismatic wall in clockwise direction
        points[0] = Vec3D(halfWidth, 0.0, mid.Z + contractionHeight);
        points[1] = Vec3D(halfWidth - contractionWidth, 0.0, mid.Z);
        points[2] = Vec3D(halfWidth, 0.0, mid.Z - contractionHeight);
        w2.createOpenPrism(points);
        wallHandler.copyAndAddObject(w2);
        //! [CST:neck]

        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0, 0, -1), Vec3D(0, 0, mid.Z));
        wallHandler.copyAndAddObject(w0);

        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        for (Mdouble z = mid.Z + contractionHeight;
            particleHandler.getNumberOfObjects() <= N;
            z += 2.0 * maxParticleRadius)
        {
            for (Mdouble r = halfWidth - maxParticleRadius; r > 0; r -= 1.999 * maxParticleRadius)
            {
                for (Mdouble c = 2.0 * maxParticleRadius; c <= 2 * constants::pi * r; c += 2.0 * maxParticleRadius)
                {
                    p0.setRadius(random.getRandomNumber(minParticleRadius, maxParticleRadius));
                    p0.setPosition(Vec3D(mid.X + r * mathsFunc::sin(c / r),
                                         mid.Y + r * mathsFunc::cos(c / r),
                                         z + p0.getRadius()));
                    const Mdouble vz = random.getRandomNumber(-1, 0);
                    const Mdouble vx = random.getRandomNumber(-1, 1);
                    p0.setVelocity(Vec3D(vx,0.0,vz));
                    particleHandler.copyAndAddObject(p0);
                }
            }
        }
    }

    //Initially, a wall is inserted in the neck of the hourglass to prevent particles flowing through.
    //This wall is moved to form the base of the hourglass at time 0.9
    void actionsAfterTimeStep() override {
        if (getTime() < 0.9 && getTime() + getTimeStep() > 0.9)
        {
            logger(INFO, "Shifting bottom wall downward");
            dynamic_cast<InfiniteWall*>(wallHandler.getLastObject())->set(Vec3D(0, 0, -1), Vec3D(0, 0, getZMin()));
        }
    }

    Mdouble contractionWidth;
    Mdouble contractionHeight;
    Mdouble minParticleRadius;
    Mdouble maxParticleRadius;
    unsigned int N;
};

int main(int argc, char *argv[])
{
    logger(INFO, "Hourglass Simulation");
    // note: the scaling of the hourglass is based on Stefan Luding's demo code MDCLR/DEMO/W7

    //all parameters should be set in the main file
    //here, we use SI units (generally, other, scaled units are possible)

    //create an instance of the class and name it
    HourGlass HG;
    HG.setName("HourGlass3DDemo");


    //these parameters are needed in setupInitialConditions()

    //specify geometry
    //specify dimensions of the hourglass
    const Mdouble width = 10e-2; // 10cm
    const Mdouble height = 60e-2; // 60cm
    //set domain accordingly (domain boundaries are not walls!)
    HG.setXMin(0.0);
    HG.setXMax(width);
    HG.setYMin(0.0);
    HG.setYMax(width);
    HG.setZMin(0.0);
    HG.setZMax(height);

    //specify how big the wedge of the contraction should be
    const Mdouble contractionWidth = 2.5e-2; //2.5cm
    const Mdouble contractionHeight = 5e-2; //5cm
    HG.contractionWidth = contractionWidth;
    HG.contractionHeight = contractionHeight;

    //specify particle properties

    //these parameters are needed in setupInitialConditions()
    HG.minParticleRadius = 6e-3; // 6mm
    HG.maxParticleRadius = 10e-3; //10mm

    //specify body forces
    HG.setGravity(Vec3D(0.0, 0.0, -9.8));

    //make the species of the particle and wall
    LinearViscoelasticFrictionSpecies species;
    species.setDensity(2000);
    //specify contact properties
    //normal forces
    species.setStiffness(1e5);
    species.setDissipation(0.63);
    //tangential (sliding) forces
    species.setSlidingFrictionCoefficient(0.5);
    species.setSlidingStiffness(1.2e4);
    species.setSlidingDissipation(0.16);
    //tangential (rolling) torques
    species.setRollingFrictionCoefficient(0.2);
    species.setRollingStiffness(1.2e4);
    species.setRollingDissipation(6.3e-2);
    //normal (torsion/spin) torques
    species.setTorsionFrictionCoefficient(0.1);
    species.setTorsionStiffness(1.2e4);
    species.setSlidingDissipation(6.3e-2);
    HG.speciesHandler.copyAndAddObject(species);

    //test normal forces
    const Mdouble minParticleMass = species.getDensity() * 4.0 / 3.0 * constants::pi * mathsFunc::cubic(HG.minParticleRadius);
    //Calculates collision time for two copies of a particle of given dissipation_, k, effective mass
    logger(INFO, "minParticleMass = %", minParticleMass);
    //Calculates collision time for two copies of a particle of given dissipation_, k, effective mass
    const Mdouble tc = species.getCollisionTime(minParticleMass);
    logger(INFO, "tc  = %", tc);
    //Calculates restitution coefficient for two copies of given dissipation_, k, effective mass
    const Mdouble r = species.getRestitutionCoefficient(minParticleMass);
    logger(INFO, "restitution coefficient = %", r);

    //set other simulation parameters
    HG.setTimeStep(tc / 50.0);
    HG.setTimeMax(10.0);
    HG.setSaveCount(500);
    HG.setXBallsAdditionalArguments("-v0 -solidf");
    HG.N = 600; //number of particles
    //uncomment next two line 2 to create paraview files
    HG.setParticlesWriteVTK(true);
    HG.setWallsWriteVTK(FileType::ONE_FILE);

    HG.solve(argc, argv);
    return 0;
}

