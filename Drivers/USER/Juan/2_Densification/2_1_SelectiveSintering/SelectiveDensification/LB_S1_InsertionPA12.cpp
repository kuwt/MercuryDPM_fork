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
#include <Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h>

#include <Walls/InfiniteWall.h>
#include "Boundaries/CubeInsertionBoundary.h"
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include "Walls/TriangleWall.h"
#include <Boundaries/CubeDeletionBoundary.h>

using constants::pi;
using mathsFunc::cubic;

class Insertion_S1 : public Mercury3D
{
private:

    ThermalSinterLinFrictionReversibleAdhesiveSpecies* wallSpecies; //pointer to the wall properties

    CubeInsertionBoundary* insb;
//    DeletionBoundary* delb;

    double vol_inserted;
    bool removed_insb = false;

    Mdouble radius_ = 0.0;
    Mdouble volumeOfParticles = 0.0;
    Mdouble volFraction_ = 0.0;

    Mdouble meanCoordinationNumber = 0.0;
    Mdouble scalarNormalForce = 0.0;

public:

    Insertion_S1 (Mdouble radius, ParticleSpecies& particleSpecies, Mdouble volFraction)
    {
        radius_ = radius;
        volFraction_ = volFraction;
        //-------------------------------------------------
        speciesHandler.copyAndAddObject(particleSpecies); //particle-particle interactions

        //WallSpecies. This parameters are always constant.
        wallSpecies = speciesHandler.copyAndAddObject(ThermalSinterLinFrictionReversibleAdhesiveSpecies());
        wallSpecies-> setDensity(4000);
//        wallSpecies->setDissipation(particleSpecies->getDissipation());
        const Mdouble YoungM = 1.65e9; //[Pa] Young's Modulus
        Mdouble k_wall = YoungM;
        wallSpecies->setPlasticParameters(100*k_wall,200*k_wall, 0.0, 0.001);

        //wallParticleSpecies = speciesHandler.getMixedObject(particleSpecies,wallSpecies)->mixAll(particleSpecies,wallSpecies);
        speciesHandler.getMixedObject(0,1)->mixAll(&particleSpecies,wallSpecies); //particle-wall interactions (particle-roller)
        //-------------------------------------------------
    }

    void setupInitialConditions() override
    {
        //++++++++++Wall information++++++++++++++
        logger(INFO,"Creating a base wall at z=0");
        InfiniteWall w;
        w.setSpecies(speciesHandler.getObject(1));
        w.set({0,0,-1},{0,0,getZMin()});
        wallHandler.copyAndAddObject(w);

        w.set({0,-1,0},{0,getYMin(),0});
        wallHandler.copyAndAddObject(w);
        w.set({0,1,0},{0,getYMax(),0});
        wallHandler.copyAndAddObject(w);
        w.set({-1,0,0},{getXMin(),0,0});
        wallHandler.copyAndAddObject(w);

        w.set({1,0,0},{getXMax(),0,0});
        wallHandler.copyAndAddObject(w);


//        logger(INFO,"Creating a wall at x=Max");
//        TriangleWall wall, wall2;
//        //auto species = speciesHandler.getObject(0);
//        wall.setSpecies(speciesHandler.getObject(0));
//        wall.setVertices(Vec3D(getXMax(),getYMax(),0.0),Vec3D(getXMax(),getYMax(),getZMax()),Vec3D(getXMax(),getYMin(),getZMax()));
//
//        wall2.setSpecies(speciesHandler.getObject(0));
//        wall2.setVertices(Vec3D(getXMax(),getYMin(),getZMax()),Vec3D(getXMax(),getYMin(),0),Vec3D(getXMax(),getYMax(),0));
//
//        wallHandler.copyAndAddObject(wall);
//        wallHandler.copyAndAddObject(wall2);

        //++++++++++To store particle information++++++++++++++
        //first way to add particles to the system
        logger(INFO,"Adding particles ...");
        hGridRebuild();
        ThermalParticle insertionBoundaryParticle;
        insertionBoundaryParticle.setTemperature(0.0);
        insertionBoundaryParticle.setRadius(radius_);

        insb = new CubeInsertionBoundary;
        insb->set(insertionBoundaryParticle,
                  100,
                  Vec3D(getXMin(), getYMin(), getZMin()),
                  Vec3D(getXMax(), getYMax(), getZMax()),
                  Vec3D(0, 0, 0),
                  Vec3D(0, 0, 0)
        );
        insb = boundaryHandler.copyAndAddObject(insb);

        //--------------------------------------------------
        //Second way to add particles to the system
        //logger(INFO,"Adding particles ...");

//        hGridRebuild();
//        ThermalParticle P;
//        P.setRadius(radius_);
//        P.setSpecies(speciesHandler.getObject(0));
//        Mdouble volumeOfParticlesMax;
//        Mdouble domainRadius = getXMax();
//        Mdouble domainHeight = getZMax();
//
//        volumeOfParticlesMax =  volFraction_*getXMax()*getYMax()*domainHeight;
//
//        Mdouble maxHeight=getZMin()+P.getRadius();
//        while (volumeOfParticles < volumeOfParticlesMax)
//        {
//            Vec3D pos;
//            pos.Z = random.getRandomNumber(getZMin()+P.getRadius(), maxHeight);
//            pos.X = random.getRandomNumber(-domainRadius+P.getRadius(), domainRadius-P.getRadius());
//            pos.Y = random.getRandomNumber(-domainRadius+P.getRadius(), domainRadius-P.getRadius());
//            P.setPosition(pos);
//
//            if (checkParticleForInteraction(P))
//            {
//                particleHandler.copyAndAddObject(P);
//                volumeOfParticles += P.getVolume();
//                P.setRadius(random.getRandomNumber(0.8,1.0)*radius_);
//                if (particleHandler.getNumberOfObjects()%100==0) std::cout << "." << std::flush;
//            }
//            else {
//                maxHeight+=0.0005*radius_;
//            }
//        }
//        std::cout << std::endl;
//        logger(INFO,"Inserted % particles",particleHandler.getNumberOfObjects());
//        logger(INFO,"Total Volume particles",volumeOfParticles);
//        std::cout << std::endl;
    }

    //override continueSolve function such that the code stops when the packing is relaxed (Ekin<1e-5*Eela)
    bool continueSolve() const override
    {
        static unsigned int counter = 0;
        if (++counter>100)
        {
            counter=0;
            if (getKineticEnergy()<1e-8*getElasticEnergy())
                return false;
        }

        return true;
    }

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

    //--------------------------------------------------
    void actionsAfterTimeStep() override
    {
        //To compute the normal force
        for (auto i : interactionHandler)
        {
            scalarNormalForce += Vec3D::dot(i->getForce(),i->getNormal());
        }

        //To get the coordination number
        for (int i = particleHandler.getNumberOfObjects()-1; i>=0; i--)
        {
            meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
        }
        meanCoordinationNumber /= particleHandler.getNumberOfObjects();


        //Mdouble volSystem = constants::pi*(std::pow(getXMax()/2,2))*getZMax();
        Mdouble volSystem = (getXMax()-getXMin())* (getYMax()-getYMin())*(getZMax()-getZMin());
        Mdouble volTotalParticles = particleHandler.getVolume();
        Mdouble volumeFraction = volTotalParticles/volSystem;

        if(volumeFraction >=volFraction_ && !removed_insb)
        {
            logger(INFO,"Stop inserting particles ...");
            boundaryHandler.removeObject(insb->getIndex());
            removed_insb = true;
        }

    }

    //--------------------------------------------------
    //Function to override the output file with specific parameters
    void writeFstatHeader(std::ostream &os) const override
    {
        Mdouble volSystem = (getXMax()-getXMin())* (getYMax()-getYMin())*(getZMax()-getZMin());

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

    void printTime() const override
    {
        Mdouble massTotalParticles =  particleHandler.getMass();
        Mdouble volSystem = (getXMax()-getXMin())* (getYMax()-getYMin())*(getZMax()-getZMin());

        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
        << ", tmax=" << std::setprecision(3) << std::left << std::setw(3) << getTimeMax()
        << ", No. Particles=" << std::setprecision(3) << std::left << std::setw(6) << particleHandler.getNumberOfObjects()
        << ", MassTotalParticles=" << std::setprecision(3) << std::left<< std::setw(6)<< massTotalParticles
        << ", KineticEnergy=" << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()
//        << ", ElasticEnergy=s" << std::setprecision(3) << std::left << std::setw(6) << getElasticEnergy()
        << ", Volume inserted=" << particleHandler.getVolume()
        << ", Volume Fraction=" << particleHandler.getVolume()/volSystem
        << ", MaxHeightParticle pos=" << getSurfaceHeight()
        << std::endl;
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    //Set problem parameters:
    std::string setFilename = "LB_InsertionPA12";

    //--------------------------s--------------------
    //define properties of particle-particle contacts
    Mdouble density = 1250.0; //[kg/m^3]
    Mdouble meanRadius = 125.0e-6; //[m]//[
    Mdouble mass = density*(4.0/3.0)*pi*cubic(meanRadius); //mass of the particle

    //define domain size
    Mdouble domainLength = 500e-5; //[m]
    Mdouble domainWidth = 500e-5; //[m]
    Mdouble domainDepth = 4.0*150e-6; //[m]
    Mdouble volumeFraction = 0.59; //Stop inserting particles when reaching this volume fraction

    //--------------------------------------------------
    // Setup parameters:
    // Set-up parameters:
    const Mdouble E_Modulus = 1700e6; //[Pa] Young's modulus for PS particles
    const Mdouble eta = 0.34; //Poisson ratio
    const Mdouble pDepth = 1.35; //Penetration depth
    const Mdouble E_reduced = 1.0/((1.0-eta*eta)/E_Modulus + (1.0-eta*eta)/E_Modulus);
    const Mdouble r_effective = (meanRadius*meanRadius)/(meanRadius+meanRadius);
//    const Mdouble E_reduced = 2.3e6; //[Pa]

    Mdouble k1 = 4.0/3.0 * E_reduced*sqrt(r_effective)*1e-4;
    const Mdouble k2 = 5.0*k1; //Unloading stiffness
    const Mdouble kc = 0.0; //Cohesive stiffness
    const Mdouble restitutionCoefficient = 0.1; //

    //**********************************************************************
    //----------------------------------------------
    //create a contact law for particle-particle interactions
    ThermalSinterLinFrictionReversibleAdhesiveSpecies particleSpecies;
    //----------------------------------------------

    //-------->[Set Plastic Contribution]
    particleSpecies.setStiffnessAndRestitutionCoefficient(k1, restitutionCoefficient, mass);
    particleSpecies.setPlasticParameters(k1, k2, kc, pDepth);
    particleSpecies.setDensity(density);
    particleSpecies.setHeatCapacity(1.0);

    //----------------------------------------------
    //Create a solver and run the commands in the constructor PowderBed::PowderBed
    Insertion_S1 pb (meanRadius,particleSpecies,volumeFraction);

    //----------------------------------------------
    pb.setGravity({0,0,-9.81}); //mm/mu s^2 [Take care]
    pb.setMin({0,0,0});
    pb.setMax({domainLength,domainWidth,domainDepth});
    //----------------------------------------------

    Mdouble TimeStep = particleSpecies.getCollisionTime(mass)/50;
    //----------------------------------------------
    pb.setSaveCount(100);
    pb.setTimeStep(TimeStep);
    pb.setTimeMax(0.2);

    logger(INFO,"Time step: %", pb.getTimeStep());
    //----------------------------------------------
//    pb.wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);
    pb.setParticlesWriteVTK(false);

    pb.setName(setFilename);
    pb.setXBallsAdditionalArguments("-solidf -v0");
    pb.setFileType(FileType::ONE_FILE);

    pb.removeOldFiles();
    pb.solve();

    return 0;
}

//plot "S1_InsertionPA12.fstat" u 5:10 with lines lw 1.5 lt 1 dt 5 lc "black" title "Cylinder-500L"