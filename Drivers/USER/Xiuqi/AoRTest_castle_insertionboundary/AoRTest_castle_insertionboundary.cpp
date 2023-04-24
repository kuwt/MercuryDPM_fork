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
#include "Boundaries/CubeInsertionBoundary.h"
#include <CMakeDefinitions.h>
#include <random>

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
        const Mdouble mass_max = species->getMassFromRadius(radius_m+3*radius_dev);//approximate maximum radius
        //species->setRestitutionCoefficient(0.7, mass_max);//needs information here
        species->setSlidingStiffness(0.29*species->getStiffness());
        species->setSlidingDissipation(0.29*species->getDissipation());
        species->setSlidingFrictionCoefficient(0.4);
        species->setSlidingFrictionCoefficientStatic(0.4);

        species->setLiquidBridgeVolumeMax(mass_max * 0.2 * 2/1000.0);
        species->setLiquidBridgeVolumeMin(0.0);
        species->setDistributionCoefficient(1);
        species->setSurfaceTension(0.07);
        species->setContactAngle(constants::pi/18);
        species->setViscosity(0.00085);
        //tangential (rolling) torques
        species->setRollingStiffness(0.4*species->getStiffness());//1.2e4);
        species->setRollingDissipation(0.4*species->getDissipation());//6.3e-2);
        species->setRollingFrictionCoefficient(0.5);
        //normal (torsion/spin) torques
        species->setTorsionStiffness(0.4*species->getStiffness());//1.2e4);
        species->setSlidingDissipation(0.4*species->getDissipation());//6.3e-2);
        species->setTorsionFrictionCoefficient(0.5);

        //bottom wall
        InfiniteWall w1;
        w1.setSpecies(speciesHandler.getObject(0));
        w1.set(Vec3D(0, 0, -1), Vec3D(mid.X, mid.Y, getZMin()));
        wallHandler.copyAndAddObject(w1);

        //! [CST:outer]
        AxisymmetricIntersectionOfWalls w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.setPosition(Vec3D(mid.X, mid.Y, 0));
        //w1.setAxis(Vec3D(0, 0, 1));
        w0.setOrientation(Vec3D(0,0,1));
        w0.addObject(Vec3D(1, 0, 0), Vec3D((getXMax() - getXMin()) / 2.0, 0, getZMin()));
        w0.addObject(Vec3D(0,0,1), Vec3D((getXMax() - getXMin()) / 2.0,0, getZMin()));
        w0.addObject(Vec3D(0,0,-1),Vec3D((getXMax() - getXMin()) / 2.0,0, getZMax()));
        wallHandler.copyAndAddObject(w0);
        //! [CST:outer]


        LiquidFilmParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
//        std::default_random_engine gen;
//        std::normal_distribution<double> dis(radius_m, radius_dev);
//        p0.setRadius(dis(gen));
        particleHandler.copyAndAddObject(p0);

        CubeInsertionBoundary insertionBoundary;
        unsigned maxFail = 1; //insert as quick as possible: try every time step, until you maxFail=1 particle fails to be insertable (overlaps with another particle or wall)
        Vec3D posMin = {getXMin(),getYMin(),getZMin()*2.0/3.0};
        Vec3D posMax = {getXMax(),getYMin(),getZMax()};
        Vec3D velMin = {-0.005, 0, -0.005};
        Vec3D velMax = {0.005, 0, 0};
        insertionBoundary.set(&p0, maxFail, posMin, posMax, velMin, velMax);
        PSD psd;
        psd.setPSDFromCSV(getMercuryDPMSourceDir() + "/Drivers/USER/PSD_AoR.csv",
                          PSD::TYPE::CUMULATIVE_NUMBER_DISTRIBUTION);//mnt/C/dev/mercurydpm_git_xiuqi/   Drivers/USER/Xiuqi/PSD_INPUTDATA/PSD_AoR.csv
        insertionBoundary.setPSD(psd);
        insertionBoundary.setManualInsertion(true);
        insertionBoundary.setInitialVolume(2e-4);
        boundaryHandler.copyAndAddObject(insertionBoundary);

        //distribute the liquid volume of each particle
        for (BaseParticle* bp : particleHandler) {
            auto lfp = dynamic_cast<LiquidFilmParticle*>(bp);//cast the baseparticle to LF type
            logger.assert_debug(lfp != nullptr, "Your particles need to be of type LiquidFilmParticle for this code to work.");
            Mdouble mass_lfp = 4.0/3.0 * constants::pi * pow(lfp->getRadius(), 3.0) * species->getDensity();
            lfp->setLiquidVolume(0.095 * mass_lfp/1000.0);
        }
//        for (Mdouble z = mid.Z/2.0;
//            particleHandler.getNumberOfObjects() <= N;
//             z += 2.0 * (radius_m+3*radius_dev))
//        {
//            for (Mdouble r = halfWidth - (radius_m+3*radius_dev); r > 0; r -= 1.999 * (radius_m+3*radius_dev))
//            {
//                for (Mdouble c = 2.0 * (radius_m+3*radius_dev); c <= 2 * constants::pi * r; c += 2.0 * (radius_m+3*radius_dev))
//                {
//                    p0.setRadius(dis(gen));
//                    p0.setPosition(Vec3D(mid.X + r * mathsFunc::sin(c / r),
//                                         mid.Y + r * mathsFunc::cos(c / r),
//                                         z + p0.getRadius()));
//                    const Mdouble vz = random.getRandomNumber(-0.005, 0);
//                    const Mdouble vx = random.getRandomNumber(-0.005, 0.005);
//                    p0.setVelocity(Vec3D(vx,0.0,vz));
//
//                    Mdouble mass_lfp = 4.0/3.0 * constants::pi * mathsFunc::cubic(p0.getRadius()) * species->getDensity();
//                    p0.setLiquidVolume(0.2 * mass_lfp/1000.0);
//                    particleHandler.copyAndAddObject(p0);
//
//                }
//            }
//        }
    }

    //Initially, the cylinder is resting where it is set
    //The cylinder wall is moved upwards at time 0.1
    void actionsAfterTimeStep() override {
        if (getTime() < 0.1 && getTime() + getTimeStep() > 0.1)
        {
            logger(INFO, "Move the cylinder");
            dynamic_cast<AxisymmetricIntersectionOfWalls*>(wallHandler.getLastObject())->setVelocity(Vec3D(0,0,0.01));
        }
    }

    Mdouble radius_m, radius_dev;
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
    AoR.setName("AngleofReposeTest_castle_insertionboundary");


    //these parameters are needed in setupInitialConditions()

    //specify geometry
    //specify dimensions of the hourglass
    const Mdouble width = 0.1; // 100mm
    const Mdouble height = 0.15; // 150mm
    //set domain accordingly (domain boundaries are not walls!)
    AoR.setXMin(0.0);
    AoR.setXMax(width);
    AoR.setYMin(0.0);
    AoR.setYMax(width);
    AoR.setZMin(0.0);
    AoR.setZMax(height);

    //specify particle properties

    //these parameters are needed in setupInitialConditions()
//    AoR.radius_m = 0.00025; // 0.25mm
//    AoR.radius_dev = 0.00003;

    //specify body forces
    AoR.setGravity(Vec3D(0.0, 0.0, -9.8));

    //set other simulation parameters
    AoR.setTimeStep(4e-5);
    AoR.setTimeMax(10);
    AoR.setSaveCount(2000);


    AoR.setXBallsAdditionalArguments("-v0 -solidf");
    //AoR.N = 1800; //number of particles
    //uncomment next two line 2 to create paraview files
    AoR.setParticlesWriteVTK(true);
    AoR.setWallsWriteVTK(true);
    AoR.setInteractionsWriteVTK(true);
    AoR.setFileType(FileType::ONE_FILE);

    AoR.solve(argc, argv);
    return 0;
}

