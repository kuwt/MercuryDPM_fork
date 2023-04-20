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

#include "Mercury3D.h"
#include "Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h"
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>

/** This code creates a cylindrical container, and particles are inserted randomly. The simulation stop
 * when the equilibrium of the system is reached.
*/
class T1DensificationPolyamide : public Mercury3D{
public:

    T1DensificationPolyamide()
    {
        //Particle species and Wall species
        particleSpecies = speciesHandler.copyAndAddObject(ThermalSinterLinFrictionReversibleAdhesiveSpecies());
        wallSpecies = speciesHandler.copyAndAddObject(ThermalSinterLinFrictionReversibleAdhesiveSpecies());
        wallParticleSpecies = speciesHandler.getMixedObject(wallSpecies, particleSpecies);

        particleRadius = std::numeric_limits<double>::quiet_NaN();
        restitutionCoefficient = std::numeric_limits<double>::quiet_NaN();
    }

    void setupInitialConditions() override
    {
        //--------------------------------------------------
        //Walls species:
        InfiniteWall* baseWall = wallHandler.copyAndAddObject(InfiniteWall());
        baseWall->set(Vec3D(0.0,0.0,-1.0), Vec3D(0.0,0.0,getZMin()));
        baseWall->setSpecies(wallSpecies);

        AxisymmetricIntersectionOfWalls* sideWall = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
        sideWall->setPosition(Vec3D((getXMin()+getXMax())/2.0,(getYMin()+getYMax())/2.0,0));
        sideWall->setOrientation(Vec3D(0, 0, 1));
        sideWall->addObject(Vec3D(1,0,0),Vec3D((getXMax()-getXMin())/2.0,0,0));
        sideWall->setSpecies(wallSpecies);

        //--------------------------------------------------
        //Particle species:
        const Mdouble k1 = 40 * particleRadius; //[Stiffness depends on particle radius]
        const Mdouble pDepth = 0.1;

        Mdouble effectiveMass = particleSpecies->getMassFromRadius(particleRadius);

        particleSpecies->setStiffnessAndRestitutionCoefficient(k1, restitutionCoefficient, effectiveMass);
        particleSpecies->setPlasticParameters(k1, 30.0 * k1, k1, pDepth);

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //Friction contribution
        const Mdouble friction =  0.7;
        Mdouble k_s = 0.2*particleSpecies->getLoadingStiffness();
//        Mdouble k_r = 0.1*particleSpecies->getLoadingStiffness();
//        Mdouble k_o = 0.1*particleSpecies->getLoadingStiffness();

        Mdouble miu_s = 1.0*friction;
//        Mdouble miu_r = 0.1*friction;
//        Mdouble miu_o = 0.1*friction;

        Mdouble gamma_s = 0.2*particleSpecies->getDissipation();
//        Mdouble gamma_r = 0.05*particleSpecies->getDissipation();
//        Mdouble gamma_o = 0.05*particleSpecies->getDissipation();

        particleSpecies->setSlidingStiffness(k_s);
//        particleSpecies->setRollingStiffness(k_r);
//        particleSpecies->setTorsionStiffness(k_o);

        particleSpecies->setSlidingFrictionCoefficient(miu_s);
//        particleSpecies->setRollingFrictionCoefficient(miu_r);
//        particleSpecies->setTorsionFrictionCoefficient(miu_o);

        particleSpecies->setSlidingDissipation(gamma_s);
//        particleSpecies->setRollingDissipation(gamma_r);
//        particleSpecies->setTorsionDissipation(gamma_o);

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        //particleSpecies->setCollisionTimeAndRestitutionCoefficient(collisionTime, restitutionCoefficient, effectiveMass);
        const Mdouble collisionTime = 0.0002;
        logger(INFO, "Collision time %", collisionTime);

        setTimeStep(collisionTime/4000);

        //--------------------------------------------------
        //set wall properties based on particle species:
        wallParticleSpecies->setLoadingStiffness(particleSpecies->getLoadingStiffness());
        wallParticleSpecies->setPlasticParameters(10*k1, 100*k1, 0.0, pDepth);

        //--------------------------------------------------
        ThermalParticle P;
        P.setRadius(particleRadius);
        P.setSpecies(particleSpecies);

        Vec3D pos;
        Mdouble numberOfParticlesInserted = 0;
        Mdouble numberOfParticles = (getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin())/mathsFunc::cubic(2.0*particleRadius)/(4.0/3.14);
        std::cout << "Inserting " << std::floor(numberOfParticles) << " particles" << std::endl;
        hGridRebuild();
        while (numberOfParticlesInserted < numberOfParticles)
        {
            pos.X = random.getRandomNumber(getXMin()+P.getRadius(), getXMax()-P.getRadius());
            pos.Y = random.getRandomNumber(getYMin()+P.getRadius(), getYMax()-P.getRadius());
            pos.Z = random.getRandomNumber(getZMin()+P.getRadius(), 3.0*getZMax()-P.getRadius());
            Mdouble r=sqrt(mathsFunc::square(pos.X-sideWall->getPosition().X)+mathsFunc::square(pos.Y-sideWall->getPosition().Y));
            P.setPosition(pos);

            if (r<(getXMax()-getXMin())/2.0-P.getRadius() && checkParticleForInteraction(P))
            {
                particleHandler.copyAndAddObject(P);
                if (random.getRandomNumber(0,1)>0.5)
                {
                    P.setRadius(particleRadius);
                }
                else
                {
                    P.setRadius(particleRadius);
                }
                ++numberOfParticlesInserted;
                std::cout << ".";
            }
        }
        std::cout << std::endl;
    }

    //--------------------------------------------------
    //The simulation stops when the packing is relaxed (Ekin<1e-5*Eela)
    bool continueSolve() const override
    {
        static unsigned int counter = 0;
        if (++counter>100)
        {
            counter=0;
            if (getKineticEnergy()<1e-9*getElasticEnergy())
                return false;
        }
        return true;
    }

    //--------------------------------------------------
    //Print time and the energy of the system
    void printTime() const override
    {
        std::cout << "t=" << getTime() << " Ene " << getKineticEnergy()/getElasticEnergy() << std::endl;
    }

    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies; //pointer to the particle properties
    ThermalSinterLinFrictionReversibleAdhesiveMixedSpecies* wallParticleSpecies; //pointer to the wall-particle
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* wallSpecies; //pointer to the wall properties

    double particleRadius; //mean particle diameter
    Mdouble restitutionCoefficient; // dissipativeness of particles
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    T1DensificationPolyamide ic;
    ic.setName("DensificationPolyamideS0");
    ic.setGravity(Vec3D(0.0,0.0,-9.8));

    ic.particleRadius = 1.0e-6;
    ic.particleSpecies->setDensity(1000.0);
    ic.restitutionCoefficient = 0.00001; // dissipative particles

    ic.setXMin(0.0);
    ic.setYMin(0.0);
    ic.setZMin(0.0);
    ic.setXMax(10e-6);
    ic.setYMax(10e-6);
    ic.setZMax(15e-6);

    ic.setFileType(FileType::ONE_FILE);
    ic.setParticlesWriteVTK(true);
    ic.wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);

    ic.setSaveCount(200);
    ic.setTimeMax(0.4);
    ic.solve();

    std::cout << "Execute 'DensificationPolyamideS0.gnu' to view output" << std::endl;
    helpers::writeToFile("DensificationPolyamideS0.gnu",
                         "set xlabel 'displacement'\n"
                         "set ylabel 'force'\n"
                         "plot 'DensificationPolyamideS0.fstat' u 7:9 w lp\n"
    );

    helpers::writeToFile("EnergyDensificationPolyamide.gnu",
                         "set xlabel 't'\n"
                         "set ylabel 'Energy'\n"
                         "plot 'DensificationPolyamideS0.ene' u 1:2 w lp\n"
    );
//    gnuplot> p ’FreeFallSelfTest.ene’ u 1:2 w l title ’Ene_grav’, \
//’’ u 1:3 w l title ’Ene_kin’, ’’ u 1:5 w l title ’Ene_ela’, \
//’’ u 1:($2+$3+$5) w l title ’Ene_tot’
}
