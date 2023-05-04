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

#include "Mercury3D.h"
#include "Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"

/** This program defines a cylindrical container, then particles are inserted.
 * The simulation  stops at an equilibrium state.
*/

// Main class:
class S1_Insertion : public Mercury3D{
public:
    explicit S1_Insertion(Mdouble particleRadius)
    {
        //Global parameters
        pRadius = particleRadius;

        setName("S1_InsertionScaledMass");
        setGravity(Vec3D(0.0,0.0,-9.8));

        setFileType(FileType::ONE_FILE);
//        setParticlesWriteVTK(true);
//        wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);

        setXBallsAdditionalArguments("-solidf -v0");

        //Boundary Parameters:
        setXMin(0.0);
        setYMin(0.0);
        setZMin(0.0);
        setXMax(350e-6);
        setYMax(350e-6);
        setZMax(500e-6);

        //--------------------------------------------------

        particleSpecies = speciesHandler.copyAndAddObject(ThermalSinterLinFrictionReversibleAdhesiveSpecies());
        //Particle species
        particleSpecies-> setDensity(930.0); //https://www.eos.info/en/additive-manufacturing/3d-printing-plastic/sls-polymer-materials/polyamide-pa-12-alumide
        const Mdouble E_Modulus = 1650e6; //[Pa] http://www.matweb.com/search/DataSheet.aspx?MatGUID=0e37a459c4eb452faa9d92659f9a0ccc&ckck=1
        const Mdouble eta = 0.35; //Poisson ratio

        const Mdouble E_reduced = 1.0/(2.0*((1.0 - eta*eta)/E_Modulus + (1.0-eta*eta)/E_Modulus));
        const Mdouble Stiffness = E_reduced*pRadius;

        const Mdouble k1 = 1.0e5; //[Stiffness depends on particle radius]
        const Mdouble pDepth = 0.1;
        const Mdouble restitutionCoefficient = 0.65; //
        const Mdouble effectiveMass = particleSpecies->getMassFromRadius(pRadius);

        particleSpecies->setStiffnessAndRestitutionCoefficient(k1, restitutionCoefficient, effectiveMass);
        particleSpecies->setPlasticParameters(k1, 10.0*k1, 0.0, pDepth);
        particleSpecies->setSinterAdhesion(0.0);
        particleSpecies->setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);
//        particleSpecies->setSinterAdhesion(0.001*k1);
        //--------------------------------------------------
        //Wall species
        wallSpecies = speciesHandler.copyAndAddObject(ThermalSinterLinFrictionReversibleAdhesiveSpecies());
        wallSpecies-> setDensity(2500.0e6);
        const Mdouble YoungM = 1.65e9; //[Pa] Young's Modulus
        Mdouble k_wall = YoungM;
        wallSpecies->setPlasticParameters(k_wall,10*k_wall, 0.0, 0.001);

        //wallParticleSpecies = speciesHandler.getMixedObject(particleSpecies,wallSpecies)->mixAll(particleSpecies,wallSpecies);
        speciesHandler.getMixedObject(particleSpecies,wallSpecies)->mixAll(particleSpecies,wallSpecies);
    }

    void scaleMass(Mdouble scale)
    {
        scaleFactor = scale;
        logger(INFO,"Scaling density, inverse gravity, square of timestep, and square of diss. coeff. up by a factor %",scale);
        particleSpecies->setDensity(scale * particleSpecies->getDensity());
        setGravity(getGravity()/scale);
        //setTimeStep(sqrt(scale) * getTimeStep());
        particleSpecies->setDissipation(sqrt(scale) * particleSpecies->getDissipation());
    }

    void setupInitialConditions() override
    {
        //Walls:
        InfiniteWall baseWall;
        baseWall.setSpecies(wallSpecies);
        baseWall.set(Vec3D(0.0,0.0,-1.0), Vec3D(0.0,0.0,getZMin()));
        wallHandler.copyAndAddObject(baseWall);

        AxisymmetricIntersectionOfWalls sideWall;
        sideWall.setSpecies(wallSpecies);
        sideWall.setPosition(Vec3D((getXMin()+getXMax())/2.0,(getYMin()+getYMax())/2.0,0));
        sideWall.setOrientation(Vec3D(0, 0, 1));
        sideWall.addObject(Vec3D(1,0,0),Vec3D((getXMax()-getXMin())/2.0,0,0));
        wallHandler.copyAndAddObject(sideWall);

        //--------------------------------------------------
        hGridRebuild();
        ThermalParticle P;
        P.setRadius(pRadius);
        P.setSpecies(particleSpecies);
        Mdouble volumeOfParticlesMax;
        Mdouble domainRadius = getXMax()/2.0;
        Mdouble domainHeight = getZMax();

        volumeOfParticlesMax = 0.64*constants::pi*domainRadius*domainRadius*domainHeight;

        Mdouble maxHeight=getZMin()+P.getRadius();
        while (volumeOfParticles < volumeOfParticlesMax)
        {
            Vec3D pos;
            pos.Z = random.getRandomNumber(getZMin()+P.getRadius(), maxHeight);
            pos.X = random.getRandomNumber(-domainRadius+P.getRadius(), domainRadius-P.getRadius());
            pos.Y = random.getRandomNumber(-domainRadius+P.getRadius(), domainRadius-P.getRadius());
            P.setPosition(pos);

            if (checkParticleForInteraction(P))
            {
                particleHandler.copyAndAddObject(P);
                volumeOfParticles += P.getVolume();
                P.setRadius(random.getRandomNumber(1.0,1.0)*pRadius);
                if (particleHandler.getNumberOfObjects()%100==0) std::cout << "." << std::flush;
            }
            else {
                maxHeight+=0.05*pRadius;
            }
        }
        std::cout << std::endl;
        logger(INFO,"Inserted % particles",particleHandler.getNumberOfObjects());
        logger(INFO,"Total Volume particles",volumeOfParticles);
        std::cout << std::endl;

        // Integration variables
        Mdouble minRadius = particleHandler.getSmallestParticle()->getRadius();
        Mdouble minMass = particleSpecies->getMassFromRadius(minRadius);
        Mdouble oldGravity = getGravity().Z;
        setGravity(getGravity()*1e5);//set to 5e5/(H/d)
        const Mdouble collisionTime = particleSpecies->getCollisionTime(minMass);
        setTimeStep(4.0*particleSpecies->computeTimeStep(minMass));
        logger(INFO,"Time step: %", getTimeStep());
        logger(INFO,"Fudge factors: Gravity scaled by factor %, timestep scaled by %",getGravity().Z/oldGravity,getTimeStep()/particleSpecies->computeTimeStep(minMass));

        std::cout << "Particles inserted = " << particleHandler.getNumberOfObjects()<< std::endl;
        std::cout << "Mean radius = " << particleHandler.getMeanRadius()<< std::endl;
        std::cout << "Volume inserted = " << particleHandler.getVolume()<< std::endl;
        std::cout << "Mass inserted = " << particleHandler.getMass()/scaleFactor<< std::endl;
        std::cout << "Density inserted = " << particleHandler.getMass()/(scaleFactor*particleHandler.getVolume())<< std::endl;
    }
    //--------------------------------------------------
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

    //--------------------------------------------------
    // Returns the height of the highest particle
    Mdouble getSurfaceHeight() const
    {
        Mdouble height = 0.0;


        for (const BaseParticle* p : particleHandler)
        {
            double newHeight = p->getPosition().Z;
            if (height<newHeight) height = newHeight;
        }
        return height+particleHandler.getMeanRadius();
    }

    //--------------------------------------------------
    // To compute the mean coordination number.
    void actionsAfterTimeStep() override
    {
        newHeight = getSurfaceHeight();

        for (int i = particleHandler.getNumberOfObjects()-1; i>=0; i--)
        {
            meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
        }
        meanCoordinationNumber /= particleHandler.getNumberOfObjects();

    }
    //--------------------------------------------------
    // Print parameters of the system.
    void printTime() const override
    {

        Mdouble volSystem = constants::pi*(std::pow(getXMax()/2,2))*newHeight;
        Mdouble volTotalParticles = particleHandler.getVolume();
        Mdouble volParticle = (4.0/3.0)*constants::pi*(std::pow(particleHandler.getMeanRadius(),3.0));
        Mdouble densityParticle = particleSpecies->getDensity()/scaleFactor;

        Mdouble massParticle = densityParticle*volParticle;
        Mdouble massTotalParticles = massParticle*particleHandler.getNumberOfObjects();
        Mdouble massTotalParticles2 =  particleHandler.getMass()/scaleFactor;

        Mdouble volParticlesPlusVoidInSystem = constants::pi*(std::pow(getXMax()/2/2,2.0))*newHeight;

        Mdouble bulkDensity = massTotalParticles/(volParticlesPlusVoidInSystem);
        Mdouble relativeDensity = bulkDensity/450.0;

        Mdouble volumeFraction = volParticlesPlusVoidInSystem/volSystem;

        std::cout
        << "t=" << std::setprecision(3) << std::left<< std::setw(6)<< getTime()
        << " tmax= " << std::setprecision(3) << std::left<< std::setw(6)<< getTimeMax()
        << " t_step= " << std::setprecision(3) << std::left<< std::setw(6)<< getTimeStep()
        << " t_step+= " << std::setprecision(3) << std::left<< std::setw(6)<< getNextTime()/4
        << " volumeFractionSys=" << std::setprecision(3) << std::left<< std::setw(6)<< volTotalParticles/volSystem
        << " Porosity=" << std::setprecision(3) << std::left<< std::setw(6)<< volTotalParticles/volSystem
        << " MassTotalParticles=" << std::setprecision(3) << std::left<< std::setw(6)<< massTotalParticles2
        << " Boundary Height= " << getZMax()
        << " MaxHeightParticle pos= " << getSurfaceHeight()
        << " Loading Stiffness= " << particleSpecies->getLoadingStiffness()
        << std::endl;
    }

private:
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies; //pointer to the particle properties
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* wallSpecies; //pointer to the wall properties

    Mdouble pRadius;
    Mdouble meanCoordinationNumber = 0.0;
    Mdouble volumeOfParticles = 0.0;
    Mdouble newHeight = 0.0;
    Mdouble scaleFactor = 0.0;

};

// Main function:
int main(int argc UNUSED, char *argv[] UNUSED)
{
    Mdouble particleRadius = 33.0e-6;
    S1_Insertion ic(particleRadius);
    ic.scaleMass(1e12);

    //This helper is to see the Fn vs Overlap, rate of overlap.
    std::cout << "Execute gnuplot, 'load S1_InsertionScaledMass.gnu' to view output" << std::endl;
    helpers::writeToFile("S1_Insertion.gnu",
                         "set xlabel 'displacement [{/Symbol d}]'\n"
                         "set ylabel 'force [f^n]'\n"
                         "set grid\n"
                         "plot 'S1_Insertion.fstat' u 7:9 w lp\n"
    );

    //Simulation set-up:
    ic.setSaveCount(1000);
    ic.setTimeMax(1e20);
    ic.removeOldFiles();
    ic.solve();

    return 0;
}
