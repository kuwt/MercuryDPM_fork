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

//Tutorial 11 introduces the axisymmetrical walls. Two asymmetrical walls are created into a
//cylindrical boundary.

//! [T11:headers]
#include "Mercury3D.h"
#include "Walls/IntersectionOfWalls.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Species/LinearViscoelasticSpecies.h"
//! [T11:headers]

//! [T11:class]
class Tutorial11 : public Mercury3D
{
public:
    //! [T11:constructor]
    Tutorial11(const Mdouble width, const Mdouble height){
        logger(INFO, "Tutorial 11");
        setName("Tutorial11");
        setFileType(FileType::ONE_FILE);
        restartFile.setFileType(FileType::ONE_FILE);
        fStatFile.setFileType(FileType::NO_FILE);
        eneFile.setFileType(FileType::NO_FILE);
        setParticlesWriteVTK(true);
        setWallsWriteVTK(FileType::ONE_FILE);
        setXBallsAdditionalArguments("-v0 -solidf");

        //specify body forces
        setGravity(Vec3D(0.0, 0.0, -9.8));

        //set domain accordingly (domain boundaries are not walls!)
        setXMin(0.0);
        setXMax(width);
        setYMin(0.0);
        setYMax(width);
        setZMin(0.0);
        setZMax(height);
    }
    //! [T11:constructor]

    //! [T11:initialConditions]
    void setupInitialConditions() override {
        Vec3D mid = {
                (getXMin() + getXMax()) / 2.0,
                (getYMin() + getYMax()) / 2.0,
                (getZMin() + getZMax()) / 2.0};
        const Mdouble halfWidth = (getXMax() - getXMin()) / 2.0;

        //! [T11:Cylinder]
        //Cylinder
        AxisymmetricIntersectionOfWalls w1;
        w1.setSpecies(speciesHandler.getObject(0));
        w1.setPosition(Vec3D(mid.X, mid.Y, 0));
        w1.setAxis(Vec3D(0, 0, 1));
        w1.addObject(Vec3D(1, 0, 0), Vec3D((getXMax() - getXMin()) / 2.0, 0, 0));
        wallHandler.copyAndAddObject(w1);
        //! [T11:Cylinder]

        //! [T11: Second Asymmetric wall]
        //Cone
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
        //! [T11: Second Asymmetric wall]

        //! [T11: simple wall]
        //Flat surface
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0, 0, -1), Vec3D(0, 0, mid.Z));
        wallHandler.copyAndAddObject(w0);
        //! [T11: simple wall]

        //! [T11: particle species]
        const int N = 600; //maximum number of particles
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
        //! [T11: particle species]
    }
    //! [T11:initialConditions]

    //! [T11:functiontime]
    //Initially, a wall is inserted in the neck of the hourglass to prevent particles flowing through.
    //This wall is moved to form the base of the hourglass at time 0.9
    void actionsAfterTimeStep() override {
        if (getTime() < 0.9 && getTime() + getTimeStep() > 0.9)
        {
            logger(INFO, "Shifting bottom wall downward");
            dynamic_cast<InfiniteWall*>(wallHandler.getLastObject())->set(Vec3D(0, 0, -1), Vec3D(0, 0, getZMin()));
        }
    }
    //! [T11:functiontime]

    Mdouble contractionWidth{};
    Mdouble contractionHeight{};
    Mdouble minParticleRadius{};
    Mdouble maxParticleRadius{};

};
//! [T11:class]

//! [T11:main]
int main(int argc, char *argv[])
{
    //![T11: domain]
    Mdouble width = 10e-2; // 10cm
    Mdouble height = 60e-2; // 60cm

    //Point the object HG to class
    Tutorial11 HG(width,height);

    //Specify particle radius:
    //these parameters are needed in setupInitialConditions()
    HG.minParticleRadius = 6e-3; // 6mm
    HG.maxParticleRadius = 10e-3; //10mm
    //![T11: domain]

    //![T11: other specifications]
    //Specify the number of particles

    //specify how big the wedge of the contraction should be
    const Mdouble contractionWidth = 2.5e-2; //2.5cm
    const Mdouble contractionHeight = 5e-2; //5cm
    HG.contractionWidth = contractionWidth;
    HG.contractionHeight = contractionHeight;
    //![T11: other specifications]

    //![T11: species configuration]
    //make the species of the particle and wall
    LinearViscoelasticSpecies species;
    species.setDensity(2000);

    species.setStiffness(1e5);
    species.setDissipation(0.63);
    HG.speciesHandler.copyAndAddObject(species);

    //![T11: species configuration]

    //![T11: tests]
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
    //![T11: tests]

    //![T11: time]
    //time integration parameters
    HG.setTimeStep(tc / 10.0);
    HG.setTimeMax(5.0);
    HG.setSaveCount(500);
    //![T11: time]

    //![T11: solve]
    HG.solve(argc, argv);
    //![T11: solve]
    return 0;
}
//! [T11:main]
