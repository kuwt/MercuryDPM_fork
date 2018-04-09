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
#include <Species/LinearPlasticViscoelasticFrictionSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>

/** This code creates a cylindrical container, inserts particles and lets them settle.
*/
class InitialConditions : public Mercury3D{
public:

    //set slave variables (i.e. compute stiffness, dissipation, timestep, wall and particle positions)
    void setupInitialConditions() override
	{
        //creating pointers to species for particle-particle interactions and particle-wall interactions
        LinearPlasticViscoelasticFrictionSpecies s;
        s.setDensity(capRockDensity);
        s.setCollisionTimeAndRestitutionCoefficient(capRockCollisionTime, restitutionCoefficient, s.getDensity()*3.14/6.0*mathsFunc::cubic(capRockDiameter));
        //- set plasticity properties
        s.setUnloadingStiffnessMax(2.0 * s.getLoadingStiffness());
        s.setCohesionStiffness(s.getLoadingStiffness());
        s.setPenetrationDepthMax(0.5); //this is the fluid overlap divided by the effective particle diameter (i.e. the particle radius)
        s.setSlidingFrictionCoefficient(0.5);
        s.setSlidingStiffness(2.0/7.0*s.getLoadingStiffness());
        s.setSlidingDissipation(2.0/7.0*s.getDissipation());
        s.setRollingFrictionCoefficient(0.5);
        s.setRollingStiffness(2.0/5.0*s.getLoadingStiffness());
        s.setRollingDissipation(2.0/5.0*s.getDissipation());
        s.setTorsionFrictionCoefficient(0.5);
        s.setTorsionStiffness(2.0/5.0*s.getLoadingStiffness());
        s.setTorsionDissipation(2.0/5.0*s.getDissipation());
        capRockSpecies = speciesHandler.copyAndAddObject(s); //pointer to the particle properties
        
        s.setDensity(pureSilicaDensity);
        s.setCollisionTimeAndRestitutionCoefficient(pureSilicaCollisionTime, restitutionCoefficient, s.getDensity()*3.14/6.0*mathsFunc::cubic(pureSilicaDiameter));
        // s.setLoadingStiffness(pureSilicaStiffness);
        // s.setDissipation(pureSilicaDissipation);
        //- set plasticity properties
        s.setUnloadingStiffnessMax(2.0 * s.getLoadingStiffness());
        s.setCohesionStiffness(s.getLoadingStiffness());
        s.setPenetrationDepthMax(0.5); //this is the fluid overlap divided by the effective particle diameter (i.e. the particle radius)
        s.setSlidingStiffness(2.0/7.0*s.getLoadingStiffness());
        s.setSlidingDissipation(2.0/7.0*s.getDissipation());
        s.setRollingStiffness(2.0/5.0*s.getLoadingStiffness());
        s.setRollingDissipation(2.0/5.0*s.getDissipation());
        s.setTorsionStiffness(2.0/5.0*s.getLoadingStiffness());
        s.setTorsionDissipation(2.0/5.0*s.getDissipation());
        pureSilicaSpecies = speciesHandler.copyAndAddObject(s); //pointer to the particle properties
        
        s.setDensity(pureCalciumCarbonateDensity);
        s.setCollisionTimeAndRestitutionCoefficient(pureCalciumCarbonateCollisionTime, restitutionCoefficient, s.getDensity()*3.14/6.0*mathsFunc::cubic(pureCalciumCarbonateCollisionTime));
        //- set plasticity properties
        s.setUnloadingStiffnessMax(2.0 * s.getLoadingStiffness());
        s.setCohesionStiffness(s.getLoadingStiffness());
        s.setPenetrationDepthMax(0.5); //this is the fluid overlap divided by the effective particle diameter (i.e. the particle radius)
        s.setSlidingStiffness(2.0/7.0*s.getLoadingStiffness());
        s.setSlidingDissipation(2.0/7.0*s.getDissipation());
        s.setRollingStiffness(2.0/5.0*s.getLoadingStiffness());
        s.setRollingDissipation(2.0/5.0*s.getDissipation());
        s.setTorsionStiffness(2.0/5.0*s.getLoadingStiffness());
        s.setTorsionDissipation(2.0/5.0*s.getDissipation());
        pureCalciumCarbonateSpecies = speciesHandler.copyAndAddObject(s); //pointer to the particle properties

        setTimeStep(capRockCollisionTime/50.0); //note the time step is slightly too large, but fine for relaxation

        // the walls are set based on xMin, xMax, ..., zMax
        InfiniteWall* baseWall = wallHandler.copyAndAddObject(InfiniteWall());
        baseWall->set(Vec3D(0.0,0.0,-1.0), Vec3D(0.0,0.0,getZMin()));
        baseWall->setSpecies(capRockSpecies);

        PeriodicBoundary* xPeriodicBoundary = boundaryHandler.copyAndAddObject(PeriodicBoundary());
        xPeriodicBoundary->set(Vec3D(1.0,0.0,0.0), getXMin(), getXMax());
        PeriodicBoundary* yPeriodicBoundary = boundaryHandler.copyAndAddObject(PeriodicBoundary());
        yPeriodicBoundary->set(Vec3D(0.0,1.0,0.0), getYMin(), getYMax());
        
        // check if the mean particle diameter is set
        if (capRockDiameter==0.0 || pureSilicaDiameter==0.0 || pureCalciumCarbonateDiameter==0.0)
        {
            std::cerr << "Error in InitialConditions::setupInitialConditions: The particle diameter is not set." << std::endl;
            exit(-1);
        }

        BaseParticle P;
        Vec3D pos;
        
        // Inserting particles; 50% of the particles are large (r_s=0.6*d), 50% small (r_s=0.4*d)
        P.setRadius(0.5*capRockDiameter);
        P.setSpecies(capRockSpecies);
        Mdouble numberOfParticlesInserted = 0;
        Mdouble numberOfParticles = 0.5*0.64*(getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin())/(3.14/6.0*mathsFunc::cubic(capRockDiameter));
        std::cout << "Inserting " << std::floor(numberOfParticles) << " particles" << std::endl;
        hGridRebuild();
        while (numberOfParticlesInserted < numberOfParticles)
        {
            pos.X = random.getRandomNumber(getXMin()+P.getRadius(), getXMax()-P.getRadius());
            pos.Y = random.getRandomNumber(getYMin()+P.getRadius(), getYMax()-P.getRadius());
            pos.Z = random.getRandomNumber(1.5*getZMax()+P.getRadius(), 3.0*getZMax()-P.getRadius());
            P.setPosition(pos);

            if (checkParticleForInteraction(P))
            {
                particleHandler.copyAndAddObject(P);
                ++numberOfParticlesInserted;
                std::cout << "C";
            }
        }
        std::cout << std::endl;

        // Inserting porous media particles;
        P.setRadius(0.5*pureSilicaDiameter);
        P.setSpecies(pureSilicaSpecies);
        numberOfParticlesInserted = 0;
        numberOfParticles = 0.5*0.64*(getXMax()-getXMin())*(getYMax()-getYMin())*(getZMax()-getZMin())
            /(3.14/6.0*0.5*(mathsFunc::cubic(pureCalciumCarbonateDiameter)+mathsFunc::cubic(pureSilicaDiameter)));
        std::cout << "Inserting " << std::floor(numberOfParticles) << " particles" << std::endl;
        hGridRebuild();
        while (numberOfParticlesInserted < numberOfParticles)
        {
            pos.X = random.getRandomNumber(getXMin()+P.getRadius(), getXMax()-P.getRadius());
            pos.Y = random.getRandomNumber(getYMin()+P.getRadius(), getYMax()-P.getRadius());
            pos.Z = random.getRandomNumber(getZMin()+P.getRadius(), 1.5*getZMax()-P.getRadius());
            P.setPosition(pos);

            if (checkParticleForInteraction(P))
            {
                particleHandler.copyAndAddObject(P);
                if (random.getRandomNumber(0,1)>0.5)
                {
                    P.setRadius(0.5*pureSilicaDiameter);
                    P.setSpecies(pureSilicaSpecies);
                    std::cout << "s";
                }
                else
                {
                    P.setRadius(0.5*pureCalciumCarbonateDiameter);
                    P.setSpecies(pureCalciumCarbonateSpecies);
                    std::cout << "c";
                }
                ++numberOfParticlesInserted;
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

    LinearPlasticViscoelasticFrictionSpecies* capRockSpecies; //pointer to the particle properties
    LinearPlasticViscoelasticFrictionSpecies* pureSilicaSpecies; //pointer to the particle properties
    LinearPlasticViscoelasticFrictionSpecies* pureCalciumCarbonateSpecies; //pointer to the particle properties
    
    Mdouble capRockDiameter; //mean particle diameter
    Mdouble pureSilicaDiameter; //mean particle diameter
    Mdouble pureCalciumCarbonateDiameter; //mean particle diameter

    Mdouble capRockDensity; //
    Mdouble pureSilicaDensity; //
    Mdouble pureCalciumCarbonateDensity; //

    Mdouble capRockCollisionTime; //softness of particles
    Mdouble pureSilicaCollisionTime; //softness of particles
    Mdouble pureCalciumCarbonateCollisionTime; //softness of particles

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
    ic.capRockDiameter = 5.0e-3;
    ic.pureSilicaDiameter = 2.5e-3;
    ic.pureCalciumCarbonateDiameter = 2.0e-3;
    ic.capRockCollisionTime = 0.001; //
    ic.pureSilicaCollisionTime = 0.004; //
    ic.pureCalciumCarbonateCollisionTime = 0.005; //
    ic.restitutionCoefficient = 0.001; // dissipative particles
    ic.capRockDensity = 2500; //
    ic.pureSilicaDensity = 2600; //
    ic.pureCalciumCarbonateDensity = 2700; //
 
    ic.setXMin(0.0);
    ic.setYMin(0.0);
    ic.setZMin(0.0);
    ic.setXMax(27e-3);
    ic.setYMax(27e-3);
    ic.setZMax(27e-3);

    ic.solve();
}
