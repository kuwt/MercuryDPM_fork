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
#include <Species/LinearPlasticViscoelasticSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>

/** This code creates a cylindrical container, inserts particles and lets them settle.
*/
class InitialConditions : public Mercury3D{
public:

    //setting default values for member variables (the actual values should be set in the main function)
    InitialConditions()
    {
        //creating pointers to species for particle-particle interactions and particle-wall interactions
        particleSpecies = speciesHandler.copyAndAddObject(LinearPlasticViscoelasticSpecies());
        wallSpecies = speciesHandler.copyAndAddObject(LinearPlasticViscoelasticSpecies());
        wallParticleSpecies = speciesHandler.getMixedObject(wallSpecies, particleSpecies);

        particleDiameter = std::numeric_limits<double>::quiet_NaN();
        collisionTime = std::numeric_limits<double>::quiet_NaN();
        restitutionCoefficient = std::numeric_limits<double>::quiet_NaN();
    }

    //set slave variables (i.e. compute stiffness, dissipation, time step, wall and particle positions)
    void setupInitialConditions() override
	{
        //set contact properties:
        //- calculate stiffness and dissipation
        Mdouble effectiveMass = 0.5*particleSpecies->getMassFromRadius(0.5*particleDiameter);
        std::cout << "Mass" << 2.0*effectiveMass << std::endl;
        particleSpecies->setCollisionTimeAndRestitutionCoefficient(collisionTime, restitutionCoefficient, 2.0*effectiveMass);
        setTimeStep(collisionTime/40.0); //note the time step is slightly too large, but fine for relaxation

        //- set plasticity properties
        particleSpecies->setUnloadingStiffnessMax(5.0 * particleSpecies->getLoadingStiffness());
        particleSpecies->setCohesionStiffness(particleSpecies->getLoadingStiffness());
        particleSpecies->setPenetrationDepthMax(0.5); //this is the fluid overlap divided by the effective particle diameter (i.e. the particle radius)

        //set wall properties
        wallParticleSpecies->setLoadingStiffness(particleSpecies->getLoadingStiffness());
        wallParticleSpecies->setDissipation(particleSpecies->getDissipation());
        wallParticleSpecies->setUnloadingStiffnessMax(particleSpecies->getLoadingStiffness());
        wallParticleSpecies->setCohesionStiffness(0.0);
        wallParticleSpecies->setPenetrationDepthMax(particleSpecies->getPenetrationDepthMax());

        // the walls are set based on xMin, xMax, ..., zMax
        InfiniteWall* baseWall = wallHandler.copyAndAddObject(InfiniteWall());
        baseWall->set(Vec3D(0.0,0.0,-1.0), Vec3D(0.0,0.0,getZMin()));
        baseWall->setSpecies(wallSpecies);

        AxisymmetricIntersectionOfWalls* sideWall = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
        sideWall->setPosition(Vec3D((getXMin()+getXMax())/2.0,(getYMin()+getYMax())/2.0,0));
        sideWall->setOrientation(Vec3D(0, 0, 1));
        sideWall->addObject(Vec3D(1,0,0),Vec3D((getXMax()-getXMin())/2.0,0,0));
        sideWall->setSpecies(wallSpecies);

        // check if the mean particle diameter is set
        if (particleDiameter==0.0)
        {
            std::cerr << "Error in InitialConditions::setupInitialConditions: The particle diameter is not set." << std::endl;
            exit(-1);
        }

        // Inserting particles; 50% of the particles are large (r_s=0.6*d), 50% small (r_s=0.4*d)
        BaseParticle P;
        P.setRadius(0.6*particleDiameter);
        P.setSpecies(particleSpecies);

        Vec3D pos;
        Mdouble numberOfParticlesInserted = 0;
        Mdouble numberOfParticles = (getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin())/mathsFunc::cubic(particleDiameter)/(4.0/3.14);
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
                    P.setRadius(0.4*particleDiameter);
                }
                else
                {
                    P.setRadius(0.6*particleDiameter);
                }
                ++numberOfParticlesInserted;
                std::cout << ".";
            }
        }
        std::cout << std::endl;
    }

    //override continueSolve function such that the code stops when the packing is relaxed (Ekin<1e-5*Eela)
    bool continueSolve() const override
    {
        static unsigned int counter = 0;
        if (++counter>100)
        {
            counter=0;
            if (getKineticEnergy()<1e-5*getElasticEnergy())
                return false;
        }
        return true;
    }

    //override printTime function such that console output shows the state of relaxation (Ekin/Eela)
    void printTime() const override
    {
        std::cout << "t=" << getTime() << " Ene " << getKineticEnergy()/getElasticEnergy() << std::endl;
    }

    LinearPlasticViscoelasticSpecies* particleSpecies; //pointer to the particle properties
    LinearPlasticViscoelasticMixedSpecies* wallParticleSpecies; //pointer to the wall-particle collision properties
    LinearPlasticViscoelasticSpecies* wallSpecies; //pointer to the wall properties

    double particleDiameter; //mean particle diameter
    Mdouble collisionTime; //softness of particles
    Mdouble restitutionCoefficient; // dissipativeness of particles
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
	InitialConditions ic;
    ic.setName("InitialConditions");
    ic.setGravity(Vec3D(0.0,0.0,-9.8));
    ic.setFileType(FileType::ONE_FILE);
    ic.setXBallsAdditionalArguments(" -v0 -solidf ");
    ic.setSaveCount(200);
    ic.setTimeMax(1e20); //run forever
    ic.particleDiameter = 2.0e-3;
    ic.particleSpecies->setDensity(2000.0);
    ic.collisionTime = 0.002; //soft particles
    ic.restitutionCoefficient = 0.001; // dissipative particles
    
    ic.setXMin(0.0);
    ic.setYMin(0.0);
    ic.setZMin(0.0);
    ic.setXMax(27e-3);
    ic.setYMax(27e-3);
    ic.setZMax(27e-3);

    ic.solve();
}
