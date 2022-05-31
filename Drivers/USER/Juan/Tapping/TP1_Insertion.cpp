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
#include "Walls/InfiniteWall.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include <Boundaries/CubeDeletionBoundary.h>


/** This program defines a cylindrical container, then particles are inserted.
 * The simulation  stops at an equilibrium state.
*/

using constants::pi;
using mathsFunc::cubic;

// Main class:
class S1_Insertion : public Mercury3D{

private:

    int setTypeOfInsertion = 1 ; //0: Random, 1:CubicInsertion
    int setTypeOfDomain = 0 ; //0: cylindrical, 1: cubical
    CubeInsertionBoundary* insb;

    Mdouble radius_ = 0.0;

    Mdouble volumeOfParticles = 0.0;
    Mdouble newHeight = 0.0;

    Mdouble meanCoordinationNumber = 0.0;
    Mdouble scalarNormalForce = 0.0;

    Mdouble volFraction_ = 0.0;

    bool removed_insb = false;

    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies_; //pointer to the particle properties
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* wallSpecies; //pointer to the wall properties

public:
    explicit S1_Insertion(Mdouble radius, ParticleSpecies& particleSpecies, Mdouble volumeFraction)
    {
        volFraction_ = volumeFraction;
        radius_ = radius;
        //-------------------------------------------------
        speciesHandler.copyAndAddObject(particleSpecies); //particle-particle interactions
        particleSpecies_ = dynamic_cast<ThermalSinterLinFrictionReversibleAdhesiveSpecies*>(speciesHandler.getObject(0));

        //WallSpecies. This parameters are always constant.
        wallSpecies = speciesHandler.copyAndAddObject(ThermalSinterLinFrictionReversibleAdhesiveSpecies());
        wallSpecies-> setDensity(4000);
//        wallSpecies->setDissipation(particleSpecies->getDissipation());
        const Mdouble YoungM = 1.65e9; //[Pa] Young's Modulus
        Mdouble k_wall = YoungM;
        wallSpecies->setPlasticParameters(100*k_wall,200*k_wall, 0.0, 0.001);
        speciesHandler.getMixedObject(0,1)->mixAll(&particleSpecies,wallSpecies); //particle-wall interactions (particle-roller)
    }

    void setupInitialConditions() override
    {
        //Walls:
        InfiniteWall w;
        w.setSpecies(wallSpecies);
        w.set({0,0,-1},{0,0,getZMin()});
        wallHandler.copyAndAddObject(w);

        if(setTypeOfDomain == 0)
        {
            logger(INFO,"Cilindrical domain ...");
            AxisymmetricIntersectionOfWalls sideWall;
            sideWall.setSpecies(wallSpecies);
            sideWall.setPosition(Vec3D((getXMin()+getXMax())/2.0,(getYMin()+getYMax())/2.0,0));
            sideWall.setOrientation(Vec3D(0, 0, 1));
            sideWall.addObject(Vec3D(1,0,0),Vec3D((getXMax()-getXMin())/2.0,0,0));
            wallHandler.copyAndAddObject(sideWall);
        }else if(setTypeOfDomain == 1){
            logger(INFO,"Cubical domain ...");
            w.set({0,-1,0},{0,getYMin(),0});
            wallHandler.copyAndAddObject(w);
            w.set({0,1,0},{0,getYMax(),0});
            wallHandler.copyAndAddObject(w);
            w.set({-1,0,0},{getXMin(),0,0});
            wallHandler.copyAndAddObject(w);
            w.set({1,0,0},{getXMax(),0,0});
            wallHandler.copyAndAddObject(w);
        }else{
            logger(INFO,"Error, please select type of boundary between 0 or 1 ...");
            exit(-1);
        }

        //--------------------------------------------------
        //++++++++++To store particle information++++++++++++++
        logger(INFO,"Adding particles ...");
        hGridRebuild();

        if(setTypeOfInsertion ==0)
        {
            logger(INFO,"Adding particles randomly ...");
            ThermalParticle P;
            P.setSpecies(speciesHandler.getObject(0));
            P.setRadius(radius_);

            Mdouble volumeOfParticlesMax;

            Mdouble domainRadius = (getXMax()-getXMin())/2;

            Mdouble domainHeight = (getZMax()-getZMin());

            Mdouble domainLength = (getXMax()-getXMin());
            Mdouble domainwidth = (getYMax()-getYMin());

            if(setTypeOfDomain == 0)
            {   //Cylinder
                volumeOfParticlesMax =  volFraction_*constants::pi*(std::pow(domainRadius,2.0))*domainHeight;
            }else if(setTypeOfDomain == 1)
            {   //Cube
                volumeOfParticlesMax =  volFraction_*domainLength*domainwidth*domainHeight;
            }

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
                    P.setRadius(random.getRandomNumber(1.0,1.0)*radius_);
                    if (particleHandler.getNumberOfObjects()%100==0) std::cout << "." << std::flush;
                }
                else {
                    maxHeight+=0.0001*radius_;
                }
            }
            std::cout << std::endl;
            logger(INFO,"Inserted % particles",particleHandler.getNumberOfObjects());
            logger(INFO,"Total Volume particles",volumeOfParticles);
            std::cout << std::endl;

        }else if(setTypeOfInsertion==1)
        {
            logger(INFO,"Adding particles using cubic insertion boundary ...");
            ThermalParticle particle;
            insb = new CubeInsertionBoundary;

            Vec3D max = getMax();

            insb->set(&particle,100,getMin(),Vec3D(getXMax(), getYMax(), 1.05*getZMax()),Vec3D(0, 0, 0),Vec3D(0, 0, 0));
            PSD psd;
            psd.setDistributionUniform(0.8*radius_,1.2*radius_, 50);
            insb->setPSD(psd);
            insb = boundaryHandler.copyAndAddObject(insb);

        }else{
            logger(INFO,"Error, please select type of insertion between 0 or 1 ...");
            exit(-1);
        }
    }
    //--------------------------------------------------
    //override continueSolve function such that the code stops when the packing is relaxed (Ekin<1e-5*Eela)
    bool continueSolve() const override
    {
        static unsigned int counter = 0;
        if (++counter>100)
        {
            counter=0;
            if (getKineticEnergy()<1e-6*getElasticEnergy())
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
        return height;
    }

    void actionsBeforeTimeStep()override {
        if(getKineticEnergy()<1e-6*getElasticEnergy() && !removed_insb && setTypeOfInsertion==1)
        {
            logger(INFO,"Stop inserting particles ...");
            boundaryHandler.removeObject(insb->getIndex());
            removed_insb = true;
        }
    }

    //--------------------------------------------------
    // To compute the mean coordination number.
    void actionsAfterTimeStep() override
    {
        newHeight = getSurfaceHeight();

        //To compute the normal force
        for (auto i : interactionHandler)
        {
            scalarNormalForce += Vec3D::dot(i->getForce(),i->getNormal());
        }

        //To compute the coordination number
        for (int i = particleHandler.getNumberOfObjects()-1; i>=0; i--)
        {
            meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
        }
        meanCoordinationNumber /= particleHandler.getNumberOfObjects();

        Mdouble volSystem = constants::pi*(std::pow(getXMax()/2,2))*getZMax();
        Mdouble volTotalParticles = particleHandler.getVolume();
        Mdouble volumeFraction = volTotalParticles/volSystem;

        //To stop to problem when the maximum volume fraction is reached
        if(volumeFraction >=volFraction_ && !removed_insb && setTypeOfInsertion==1)
        {
            logger(INFO,"Stop inserting particles ...");
            boundaryHandler.removeObject(insb->getIndex());
            removed_insb = true;
        }else if(getKineticEnergy()<1e-6*getElasticEnergy() &&!removed_insb && setTypeOfInsertion==1){
            logger(INFO,"Stop inserting particles ...");
            boundaryHandler.removeObject(insb->getIndex());
            removed_insb = true;
        }

    }

    //--------------------------------------------------
    //Function to override the output file with specific parameters
    void writeFstatHeader(std::ostream &os) const override
    {
        Mdouble volSystem = constants::pi*(std::pow(getXMax()/2,2))*getZMax();

        Mdouble volParticle = (4.0/3.0)*constants::pi*(std::pow(particleHandler.getMeanRadius(),3.0));
        Mdouble densityParticle = speciesHandler.getObject(0)->getDensity();

        Mdouble volTotalParticles = particleHandler.getVolume();
        Mdouble volumeFraction = volTotalParticles/volSystem;

        Mdouble massParticle = densityParticle*volParticle;
        Mdouble massTotalParticles = particleHandler.getMass();

        //Mdouble volParticlesPlusVoidInSystem = constants::pi*(std::pow(getXMax()/2,2.0))*getSurfaceHeight();
        //Mdouble bulkDensity = massTotalParticles/volParticlesPlusVoidInSystem;
        //Mdouble TheoDensity = 550.0;

        os << getTime() //[1] Current time
           << " " << getTimeMax() //[2] Max time
           << " " << volSystem //[3] Volume of the system
           << " " << getZMax() //[4] Height of the system
           << " " << particleHandler.getNumberOfObjects() //[5] Number of particles inserted
           << " " << particleHandler.getMeanRadius() // [6] Mean radius of a particle
           << " " << volParticle //[7] Mean volume of a particle
           << " " << densityParticle //[8] Mean density of a particle
           << " " << volTotalParticles //[9] Volume of particles inserted
           << " " << volumeFraction //[10] Volume fraction
           << " " << massTotalParticles //[11] total mass of particle inserted
           << " " << getSurfaceHeight() // [12] Max height reached by particles.
           << " " << getKineticEnergy() // [13] kinetic energy
           << " " << getElasticEnergy() //[14] Elastic energy
           << " " << scalarNormalForce //[15] Normal force
           << " " << meanCoordinationNumber // [16] meanCoordinationNumber
           << std::endl;
    }

    //--------------------------------------------------
    // Print parameters of the system.
    void printTime() const override
    {

        Mdouble volSystem = constants::pi*(std::pow(getXMax()/2,2))*getZMax();
        Mdouble volTotalParticles = particleHandler.getVolume();
        Mdouble volParticle = (4.0/3.0)*constants::pi*(std::pow(particleHandler.getMeanRadius(),3.0));
        Mdouble densityParticle = particleSpecies_->getDensity();

        Mdouble massParticle = densityParticle*volParticle;
        Mdouble massTotalParticles = massParticle*particleHandler.getNumberOfObjects();
        Mdouble massTotalParticles2 =  particleHandler.getMass();

        Mdouble volParticlesPlusVoidInSystem = constants::pi*(std::pow(getXMax()/2/2,2.0))*newHeight;

        Mdouble bulkDensity = massTotalParticles/(volParticlesPlusVoidInSystem);
        Mdouble relativeDensity = bulkDensity/450.0;

        Mdouble volumeFraction = volParticlesPlusVoidInSystem/volSystem;

        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << ", tmax=" << std::setprecision(3) << std::left << std::setw(3) << getTimeMax()
                  << ", No. Particles=" << std::setprecision(3) << std::left << std::setw(6) << particleHandler.getNumberOfObjects()
                  << ", MassTotalParticles=" << std::setprecision(3) << std::left<< std::setw(6)<< massTotalParticles2
                  << ", KineticEnergy=" << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()
                  << ", Volume inserted=" << particleHandler.getVolume()
                  << ", Volume Fraction=" << particleHandler.getVolume()/volSystem
                  << ", MaxHeightParticle pos=" << getSurfaceHeight()
                  << std::endl;
    }
};

// Main function:
int main(int argc UNUSED, char *argv[] UNUSED)
{
    //Set problem parameters:
    std::string setFilename = "TP1_PA12";

    Mdouble meanRadius = 30.0e-6; //

    //define domain size
    Mdouble domainLength = 6.0*meanRadius; //[m]
    Mdouble domainWidth = 6.0*meanRadius; //[m]
    Mdouble domainDepth = 4.0*domainLength; //[m]
    Mdouble volumeFraction = 0.58;

    //----------------------------------------------
    //define properties of particle-particle contacts
    Mdouble density = 1000.0; //[kg/m^3]

    Mdouble mass = density*(4.0/3.0)*pi*cubic(meanRadius); //mass of the particle
    //--------------------------------------------------
    // Setup parameters:
    const Mdouble E_Modulus = 2e4; //Young's modulus
    const Mdouble eta = 0.34; //Poisson ratio
    const Mdouble pDepth = 0.01; //Penetration depth
    const Mdouble E_reduced = 1.0/((1.0-eta*eta)/E_Modulus + (1.0-eta*eta)/E_Modulus);
    const Mdouble r_effective = (meanRadius*meanRadius)/(meanRadius+meanRadius);

    const Mdouble k1 = 4.0/3.0 * E_reduced*sqrt(r_effective);
    const Mdouble k2 = 10.0*k1; //Unloading stiffness
    const Mdouble kc = 0.0; //Cohesive stiffness
    const Mdouble restitutionCoefficient = 0.01; //

    //**********************************************************************
    //----------------------------------------------
    //create a contact law for particle-particle interactions
    ThermalSinterLinFrictionReversibleAdhesiveSpecies particleSpecies;

    //-------->[Set Plastic Contribution]
//    particleSpecies.setStiffnessAndRestitutionCoefficient(k1, restitutionCoefficient, mass);
    particleSpecies.setDensity(density);
    particleSpecies.setPlasticParameters(k1, k2, kc, pDepth);
    particleSpecies.setDissipation(0.00001);


    //Create a solver and run the commands in the constructor PowderBed::PowderBed
    S1_Insertion pb (meanRadius,particleSpecies,volumeFraction);

    Mdouble TimeStep = particleSpecies.getCollisionTime(mass)/50;
    //----------------------------------------------
    pb.setGravity({0,0,-9.81}); //mm/mu s^2 [Take care]
    pb.setMin({0,0,0});

    pb.setMax({2.0*domainLength,2.0*domainLength,domainDepth});
    pb.setMin({0.0,0.0,0.0});
    //----------------------------------------------

    pb.setSaveCount(100);
    pb.setTimeStep(TimeStep);
    pb.setTimeMax(2.0);
    //----------------------------------------------

    logger(INFO,"Time step: %", pb.getTimeStep());
    //----------------------------------------------
    //    pb.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    pb.setParticlesWriteVTK(false);
    pb.setWallsWriteVTK(false);

    pb.setName(setFilename);
    pb.setXBallsAdditionalArguments("-solidf -v0");
    pb.setFileType(FileType::ONE_FILE);

    pb.removeOldFiles();
    pb.solve();


    //This helper is to plot results in gnuplot.
    std::cout << "Execute gnuplot, 'plot S1_Insertion.fstat' u1:2 to view output" << std::endl;

    //plot "S1_Insertion.fstat" u 5:10 with lines lw 1.5 lt 1 dt 5 lc "black" title "Cylinder-500L"

    //plot "S1_Insertion.fstat" u 1:10 with linespoints ls 6 lc "black" pi 6 title "Cylinder"
    //set yrange [0:0.8]
    //set xlabel "time [s]"
    //set ylabel "Volume Fraction"

    //plot "S1_Insertion.fstat" u 5:10 with lines lw 1.5 lt 1 dt 5 lc "black" title "Cylinder-500L"
    //replot "S1_Insertion_428L.fstat" u 5:10 with lines lw 1.5 lt 1 dt 3 lc "black" title "Cylinder-428L"
    //replot "S1_Insertion_545L.fstat" u 5:10 with lines lw 1.5 lc "black" title "Cylinder-545L"

    return 0;
}