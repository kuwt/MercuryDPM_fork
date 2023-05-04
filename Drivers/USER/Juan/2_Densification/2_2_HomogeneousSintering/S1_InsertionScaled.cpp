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
#include <Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include "Boundaries/CubeInsertionBoundary.h"
#include <Boundaries/CubeDeletionBoundary.h>

using constants::pi;
using mathsFunc::cubic;

namespace units{
    //mass: 1 ug = 1 kg
    Mdouble mUnit = 1e-9;
    //length: 1 mm = 1e-3 m
    Mdouble lUnit = 1e-3;
    //time 1 ms = 1e-3 s
    Mdouble tUnit = 1e-3;
    //force: 1 uN = 1e-6N
    inline Mdouble fUnit(){return mUnit*lUnit/pow(tUnit,2);}
    //stiffness: 1 N/m = 1e-3 uN/mm
    inline Mdouble kUnit(){return fUnit()/lUnit;}
    //Young Modulus 1 N/m^2 = 1.0 uN/mm^2
    inline Mdouble youngUnit(){return fUnit()/pow(lUnit,2);}
    //stress: ug/(mm*ms^2) = 1 Pa
    inline Mdouble sigUnit(){return mUnit/(lUnit*pow(tUnit,2));}
    //density: 1 ug/mm^3 = 1 kg/m^3
    inline Mdouble rhoUnit(){ return  mUnit/pow(lUnit,3);}
    //Compliance units 1 (Pa s)^-1 = (ug/(mm*ms^2) * 1e-3 ms)^-1
    inline Mdouble complianceUnit(){return sigUnit()*tUnit;}
    //Velocity
    inline Mdouble velUnit(){return lUnit/tUnit;}
    //acceleration
    inline Mdouble accUnit(){return lUnit/pow(tUnit,2);}
    //simulation name
    inline Mdouble energyUnit(){return fUnit()* lUnit;}
    std::string name;
}

class NeckGrowth : public Mercury3D
{
private:

    Mdouble radius_ = 0.0;
    Mdouble newHeight = 0.0;
    Mdouble deltaR = 1.0e-4; //[check] Thermal expansion coefficient

    CubeInsertionBoundary* insb;
    DeletionBoundary* delb;

    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies_; //pointer to the particle properties
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* wallSpecies; //pointer to the wall properties
public:

    NeckGrowth (Mdouble radius, ParticleSpecies& particleSpecies)
    {
        radius_ = radius;
        //-------------------------------------------------
        speciesHandler.copyAndAddObject(particleSpecies); //particle-particle interactions
        //-------------------------------------------------
        particleSpecies_ = dynamic_cast<ThermalSinterLinFrictionReversibleAdhesiveSpecies*>(speciesHandler.getObject(0));

        //WallSpecies. This parameters are always constant.
        wallSpecies = speciesHandler.copyAndAddObject(ThermalSinterLinFrictionReversibleAdhesiveSpecies());
        wallSpecies-> setDensity(4000.0/units::rhoUnit());
//        wallSpecies->setDissipation(particleSpecies->getDissipation());
        const Mdouble YoungM = 1.65e9/units::youngUnit(); //[Pa] Young's Modulus
        Mdouble k_wall = YoungM;
        wallSpecies->setPlasticParameters(100*k_wall,200*k_wall, 0.0, 0.001);
        speciesHandler.getMixedObject(0,1)->mixAll(&particleSpecies,wallSpecies); //particle-wall interactions (particle-roller)
    }

    void setupInitialConditions() override {

        //++++++++++To store particle information++++++++++++++
        //Walls:
        InfiniteWall w;
        w.setSpecies(wallSpecies);
        w.set({0,0,-1},{0,0,getZMin()});
        wallHandler.copyAndAddObject(w);

        //Walls:
        wallHandler.copyAndAddObject(w);
        w.set({0,-1,0},{0,getYMin(),0});
        wallHandler.copyAndAddObject(w);
        w.set({0,1,0},{0,getYMax(),0});
        wallHandler.copyAndAddObject(w);
        w.set({-1,0,0},{getXMin(),0,0});
        wallHandler.copyAndAddObject(w);
        w.set({1,0,0},{getXMax(),0,0});
        wallHandler.copyAndAddObject(w);

//        //Cylinder
//        AxisymmetricIntersectionOfWalls sideWall;
//        sideWall.setSpecies(wallSpecies);
//        sideWall.setPosition(Vec3D((getXMin()+getXMax())/2.0,(getYMin()+getYMax())/2.0,0));
//        sideWall.setOrientation(Vec3D(0, 0, 1));
//        sideWall.addObject(Vec3D(1,0,0),Vec3D((getXMax()-getXMin())/2.0,0,0));
//        wallHandler.copyAndAddObject(sideWall);

        //--------------------------------------------------
        //++++++++++To store particle information++++++++++++++
        logger(INFO,"Adding particles ...");
        hGridRebuild();
        ThermalParticle insertionBoundaryParticle;
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
    }

    //--------------------------------------------------
    Mdouble getMeanRelativeContactRadius() const
    {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanRadius = particleHandler.getMeanRadius();
        return sqrt(meanOverlap/meanRadius);
    }

    //    Function to override the output file with specific parameters
    void writeFstatHeader(std::ostream &os) const override
    {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
//        Mdouble Force = dynamic_cast<const ThermalParticle*>(particleHandler.getObject(0))->getForce().X;
        Mdouble meanRadius = particleHandler.getMeanRadius();
        Mdouble meanContactRadius = sqrt(meanOverlap/meanRadius)*meanRadius;

        //Stress at neck
        Mdouble K_s = 2.0*meanRadius/(std::pow(meanContactRadius,2.0));
        Mdouble Shrinkage = (1.0/5.0)*std::pow(getMeanRelativeContactRadius(),2.0);

        Mdouble stress = K_s*particleSpecies_->getSurfTension();

//        Mdouble rho_p = 1.0/(1 + std::pow(1 - X,3.0));

        os << getTime()*units::tUnit/1.0e-9 //1
           << " " << 1.0 //2 getTemperature()
           << " " << 1.0 //dynamic_cast<const ThermalParticle*>(particleHandler.getObject(0))->getRadius()*units::lUnit//3
           << " " << particleSpecies_->getDensity()*units::rhoUnit()//4
           << " " << 1.0//Force*units::fUnit()//5
           << " " << stress*units::sigUnit() //6
           << " " << particleHandler.getMass()*units::mUnit//7
           << " " << particleSpecies_->getLoadingStiffness()*units::kUnit()//8
           << " " << meanOverlap*units::lUnit//9
           << " " << meanContactRadius*units::lUnit //10
           << " " << particleHandler.getMeanRadius()*units::lUnit //11
           << " " << getMeanRelativeContactRadius()//12 // Neckgrowth
           << " " << 1.0//13 [J]//getThermalEnergy()
           << " " << 1.0 //14 [W/s]getThermalEnergy()/(getTime()) //14 [W/s]
           << " " << getKineticEnergy()*units::energyUnit() //15
           << " " << getElasticEnergy()*units::energyUnit() //16
           << " " << (K_s * particleSpecies_->getSurfTension())*units::sigUnit() //17 Stress at the neckTip
           << " " << Shrinkage*units::lUnit //18 shrinkage
           << " " << K_s/units::lUnit//19
           << std::endl;
    }

    //--------------------------------------------------
    // To compute the mean coordination number.
    void actionsAfterTimeStep() override {
        newHeight = getSurfaceHeight();
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

int main(int argc UNUSED, char *argv[] UNUSED)
{
    //Set problem parameters:
    std::string setFilename = "S1_PA12_500R_500L_Scaled";
    Mdouble tMax = 1e-4/units::tUnit;
    Mdouble gravity = -9.81/units::accUnit();

    Mdouble RVESize = 60.0e-6/units::lUnit;
    Mdouble meanRadius = 0.5*RVESize; //[mm]//

    //define domain size
    Mdouble domainLength = 10*RVESize; //[m]
    Mdouble domainWidth = 10*RVESize; //[m]
    Mdouble domainDepth = 10*RVESize; //[m]

    //----------------------------------------------
    //define properties of particle-particle contacts
    Mdouble density = 1000.0/units::rhoUnit(); //[ug/mm^3]
    Mdouble mass = density*(4.0/3.0)*pi*cubic(meanRadius); //mass of the particle

    //--------------------------------------------------
    // Setup parameters:
    const Mdouble E_Modulus = 1650.0e9/units::youngUnit(); //Young's modulus
    const Mdouble eta = 0.34; //Poisson ratio
    const Mdouble E_reduced = 1.0/((1.0-eta*eta)/E_Modulus + (1.0-eta*eta)/E_Modulus);
    const Mdouble r_effective = (meanRadius*meanRadius)/(meanRadius+meanRadius);
    const Mdouble k1 = 4.0/3.0 * E_reduced*sqrt(r_effective);

    const Mdouble pDepth = 1.45; //Penetration depth
    const Mdouble k2 = 10.0*k1; //Unloading stiffness
    const Mdouble kc = 0.0; //Cohesive stiffness
    const Mdouble restitutionCoefficient = 0.1; //

    //----------------------------------------------
    //Virgin. deltaC = 2.8e-7, C1 = 1.71
    //1stAged: 2.9e-7, C1= 0.82
    //2ndAgeed: deltaC = 3.00e-7, C1 = 0.497
//    Mdouble setDeltaC = 2.8e-7/units::lUnit; //[mm]
//    Mdouble setC1 = 1.710/units::complianceUnit(); // [Pa s]^-1;
//    Mdouble surfaceTension = 0.041/units::kUnit();

    //create a contact law for particle-particle interactions
    ThermalSinterLinFrictionReversibleAdhesiveSpecies particleSpecies;
    //----------------------------------------------

    //-------->[Set Plastic Contribution]
    particleSpecies.setRestitutionCoefficient(restitutionCoefficient,mass);
    particleSpecies.setPlasticParameters(k1, k2, kc, pDepth);
    particleSpecies.setDensity(density);
//    particleSpecies.setSinterAdhesion(0.0016*k1);

    //-------------------
    //Sinter parameters
//    particleSpecies.setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);  //FRENKEL OR VISCOELASTIC_CONTACT
//    particleSpecies.setComplianceZero((1.0/(2.0*E_Modulus))); //Book: Adhesive Particles
//    particleSpecies.setSurfTension(surfaceTension); // [N/m] or [J/m^2] surface energy density or surface tension.

    //Parameters to calibrate:
//    particleSpecies.setFluidity(setC1);
//    particleSpecies.setSeparationDis(setDeltaC);

    //----------------------------------------------
    //Create a solver and run the commands in the constructor PowderBed::PowderBed
    NeckGrowth pb (meanRadius,particleSpecies);

    //----------------------------------------------
    pb.setGravity({0,0,gravity}); //mm/mu s^2 [Take care]
    pb.setMin({0.0,0.0,0.0});
    pb.setMax({domainLength,domainWidth,domainDepth});
    //----------------------------------------------

    Mdouble TimeStep = (particleSpecies.getCollisionTime(mass))*10.0;

    //----------------------------------------------
    pb.setSaveCount(1);
    pb.setTimeStep(TimeStep);
    pb.setTimeMax(tMax);

    logger(INFO,"Time step: %", pb.getTimeStep());
    //----------------------------------------------
//    pb.wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);
    pb.setParticlesWriteVTK(false);

    pb.setName(setFilename);
//    pb.setXBallsAdditionalArguments("-solidf -v0");
    pb.setFileType(FileType::ONE_FILE);

    pb.removeOldFiles();
//    pb.setXBallsAdditionalArguments("-solidf -v0 -cmode 8 -cmaxset 100 ");

    pb.solve();

    std::cout << "Execute gnuplot 'load 'Scaled_NeckGrowthPA12.gnu' ' to view output" << std::endl;
    helpers::writeToFile(setFilename + ".gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'a(t)/R'\n"
                         "plot 'Scaled_NeckGrowthPA12.fstat' u 1:12 title 'DEM simulation'  with lines linestyle 2\n"
                         //                         "replot 'PA12_Aged1_3.370e-5.txt' u ($1):($2) w p ls 1 title 'Experimental data 3.37 um' \n"
                         "replot 1.0 title 'a_{o} Limit' with lines linestyle 3 \n"
                         "replot 'PA12_2.584e-5.txt' u ($1):($2) w p ls 2 title 'Experimental data 2.584' \n"
                         "replot 'PA12_2.928e-5.txt' u ($1):($2) w p ls 3 title 'Experimental data 2.982 ' \n"
                         "replot 'PA12_3.044e-5.txt' u ($1):($2) w p ls 4 title 'Experimental data 3.044 '\n "
                         "replot 'PA12_3.210e-5.txt' u ($1):($2) w p ls 4 title 'Experimental data 3.210 ' \n"
                         "replot 'PA12_3.216e-5.txt' u ($1):($2) w p ls 4 title 'Experimental data 3.210 '\n "
                         "replot 'PA12_2.93e-5_Zhao.txt' u ($1):($2) w p ls 4 title 'Experimental data 2.93 Zao '\n "
                         "replot '2AgedPA12_2.3e-05.txt' u ($1):($2) w p ls 4 title 'Experimental data First aged 2.330 '\n"
                         "replot '2AgedPA12_3.8e-05.txt' u ($1):($2) w p ls 4 title 'Experimental data First aged 3.050 '\n"
//                         "replot '2AgedPA12_3.63e-05.txt' u ($1):($2) w p ls 4 title 'Experimental data First aged 3.37'\n "

    );


    std::cout << "Execute gnuplot 'load 'Scaled_NeckGrowthPA12_2.gnu' ' to view output" << std::endl;
    helpers::writeToFile("Scaled_NeckGrowthPA12_2.gnu",
                         "set xlabel 'a(t)/R'\n"
                         "set ylabel 'Q [J]'\n"
                         "plot 'Scaled_NeckGrowthPA12.fstat' u 12:13 title 'Heat energy'  with lines linestyle 2");


    return 0;
}
////Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
////For the list of developers, see <http://www.MercuryDPM.org/Team>.
////
////Redistribution and use in source and binary forms, with or without
////modification, are permitted provided that the following conditions are met:
////  * Redistributions of source code must retain the above copyright
////    notice, this list of conditions and the following disclaimer.
////  * Redistributions in binary form must reproduce the above copyright
////    notice, this list of conditions and the following disclaimer in the
////    documentation and/or other materials provided with the distribution.
////  * Neither the name MercuryDPM nor the
////    names of its contributors may be used to endorse or promote products
////    derived from this software without specific prior written permission.
////
////THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
////ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
////WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
////DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
////DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
////(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
////LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
////ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
////(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
////SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//#include "Mercury3D.h"
//#include "Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h"
//#include "Walls/InfiniteWall.h"
//#include "Walls/AxisymmetricIntersectionOfWalls.h"
//#include "Boundaries/CubeInsertionBoundary.h"
//#include <Boundaries/CubeDeletionBoundary.h>
//
//
///** This program defines a cylindrical container, then particles are inserted.
// * The simulation  stops at an equilibrium state.
//*/
//
//using constants::pi;
//using mathsFunc::cubic;
//
//namespace units{
//    //mass: 1 ug = 1 kg
//    Mdouble mUnit = 1e-9;
//    //length: 1 mm = 1e-3 m
//    Mdouble lUnit = 1e-3;
//    //time 1 ms = 1e-3 s
//    Mdouble tUnit = 1e-3;
//    //force: 1 uN = 1e-6N
//    inline Mdouble fUnit(){return mUnit*lUnit/pow(tUnit,2);}
//    //stiffness: 1 N/m = 1e-3 uN/mm
//    inline Mdouble kUnit(){return fUnit()/lUnit;}
//    //Young Modulus 1 N/m^2 = 1.0 uN/mm^2
//    inline Mdouble youngUnit(){return fUnit()/pow(lUnit,2);}
//    //stress: ug/(mm*ms^2) = 1 Pa
//    inline Mdouble sigUnit(){return mUnit/(lUnit*pow(tUnit,2));}
//    //density: 1 ug/mm^3 = 1 kg/m^3
//    inline Mdouble rhoUnit(){ return  mUnit/pow(lUnit,3);}
//    //Compliance units 1 (Pa s)^-1 = (ug/(mm*ms^2) * 1e-3 ms)^-1
//    inline Mdouble complianceUnit(){return sigUnit()*tUnit;}
//    //Velocity
//    inline Mdouble velUnit(){return lUnit/tUnit;}
//    //acceleration
//    inline Mdouble accUnit(){return lUnit/pow(tUnit,2);}
//    //simulation name
//    inline Mdouble energyUnit(){return fUnit()* lUnit;}
//    std::string name;
//}
//
//// Main class:
//class S1_Insertion : public Mercury3D{
//
//private:
//
//    CubeInsertionBoundary* insb;
//    DeletionBoundary* delb;
//
//    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies_; //pointer to the particle properties
//    ThermalSinterLinFrictionReversibleAdhesiveSpecies* wallSpecies; //pointer to the wall properties
//
//    Mdouble radius_ = 0.0;
//
//    Mdouble volumeOfParticles = 0.0;
//    Mdouble newHeight = 0.0;
//
//    Mdouble meanCoordinationNumber = 0.0;
//    Mdouble scalarNormalForce = 0.0;
//
//    Mdouble volFraction_ = 0.0;
//
//    bool removed_insb = false;
//
//public:
//    explicit S1_Insertion(Mdouble radius, ParticleSpecies& particleSpecies, Mdouble volumeFraction)
//    {
//	volFraction_ = volumeFraction;
//        radius_ = radius;
//        //-------------------------------------------------
//        speciesHandler.copyAndAddObject(particleSpecies); //particle-particle interactions
//        particleSpecies_ = dynamic_cast<ThermalSinterLinFrictionReversibleAdhesiveSpecies*>(speciesHandler.getObject(0));
//
//        //WallSpecies. This parameters are always constant.
//        wallSpecies = speciesHandler.copyAndAddObject(ThermalSinterLinFrictionReversibleAdhesiveSpecies());
//        wallSpecies-> setDensity(4000/units::rhoUnit());
////        wallSpecies->setDissipation(particleSpecies->getDissipation());
//        const Mdouble YoungM = 1.65e9/units::youngUnit(); //[Pa] Young's Modulus
//        Mdouble k_wall = YoungM;
//        wallSpecies->setPlasticParameters(100*k_wall,200*k_wall, 0.0, 0.001);
//        speciesHandler.getMixedObject(0,1)->mixAll(&particleSpecies,wallSpecies); //particle-wall interactions (particle-roller)
//    }
//
//    void setupInitialConditions() override
//    {
//        //Walls:
//        InfiniteWall w;
//        w.setSpecies(wallSpecies);
//        w.set({0,0,-1},{0,0,getZMin()});
//        wallHandler.copyAndAddObject(w);
//
//	    //Cylinder
//        AxisymmetricIntersectionOfWalls sideWall;
//        sideWall.setSpecies(wallSpecies);
//        sideWall.setPosition(Vec3D((getXMin()+getXMax())/2.0,(getYMin()+getYMax())/2.0,0));
//        sideWall.setOrientation(Vec3D(0, 0, 1));
//        sideWall.addObject(Vec3D(1,0,0),Vec3D((getXMax()-getXMin())/2.0,0,0));
//        wallHandler.copyAndAddObject(sideWall);
//
//        //--------------------------------------------------
//        //++++++++++To store particle information++++++++++++++
//        logger(INFO,"Adding particles ...");
//        hGridRebuild();
//        ThermalParticle insertionBoundaryParticle;
//
//        insb = new CubeInsertionBoundary;
//        insb->set(insertionBoundaryParticle,
//                100,
//                Vec3D(getXMin(), getYMin(), getZMin()),
//                Vec3D(getXMax(), getYMax(), getZMax()),
//                Vec3D(0, 0, 0),
//                Vec3D(0, 0, 0),
//                radius_*1.0,
//                radius_*1.0
//        );
//        insb = boundaryHandler.copyAndAddObject(insb);
//    }
//    //--------------------------------------------------
//    //override continueSolve function such that the code stops when the packing is relaxed (Ekin<1e-5*Eela)
//    bool continueSolve() const override
//    {
//        static unsigned int counter = 0;
//        if (++counter>100)
//        {
//            counter=0;
//            if (getKineticEnergy()<1e-5*getElasticEnergy()/units::energyUnit())
//                return false;
//        }
//        return true;
//    }
//
//    //--------------------------------------------------
//    // Returns the height of the highest particle
//    Mdouble getSurfaceHeight() const
//    {
//        Mdouble height = 0.0;
//
//        for (const BaseParticle* p : particleHandler)
//        {
//            double newHeight = p->getPosition().Z;
//            if (height<newHeight) height = newHeight;
//        }
//        return height;
//    }
//
//    void actionsBeforeTimeStep()override {
//        if(getKineticEnergy()<1e-5*getElasticEnergy()/units::energyUnit() && !removed_insb)
//        {
//            logger(INFO,"Stop inserting particles ...");
//            boundaryHandler.removeObject(insb->getIndex());
//            removed_insb = true;
//        }
//    }
//
//    //--------------------------------------------------
//    // To compute the mean coordination number.
//    void actionsAfterTimeStep() override
//    {
//        newHeight = getSurfaceHeight();
//
//        //To compute the normal force
//        for (auto i : interactionHandler)
//        {
//                scalarNormalForce += Vec3D::dot(i->getForce(),i->getNormal());
//	    }
//
//        //To compute the coordination number
//        for (int i = particleHandler.getNumberOfObjects()-1; i>=0; i--)
//        {
//            meanCoordinationNumber += (particleHandler.getObject(i) -> getInteractions()).size();
//        }
//        meanCoordinationNumber /= particleHandler.getNumberOfObjects();
//
//        Mdouble volSystem = constants::pi*(std::pow(getXMax()/2,2))*getZMax();
//        Mdouble volTotalParticles = particleHandler.getVolume();
//        Mdouble volumeFraction = volTotalParticles/volSystem;
//
//        //To stop to problem when the maximum volume fraction is reached
//        if(volumeFraction >=volFraction_ && !removed_insb)
//        {
//            logger(INFO,"Stop inserting particles ...");
//            boundaryHandler.removeObject(insb->getIndex());
//            removed_insb = true;
//        }
//
//    }
//
//    //--------------------------------------------------
//    //Function to override the output file with specific parameters
//    void writeFstatHeader(std::ostream &os) const override
//    {
//        Mdouble volSystem = constants::pi*(std::pow(getXMax()/2,2))*getZMax();
//
//        Mdouble volParticle = (4.0/3.0)*constants::pi*(std::pow(particleHandler.getMeanRadius(),3.0));
//        Mdouble densityParticle = speciesHandler.getObject(0)->getDensity();
//
//        Mdouble volTotalParticles = particleHandler.getVolume();
//        Mdouble volumeFraction = volTotalParticles/volSystem;
//
//        Mdouble massParticle = densityParticle*volParticle;
//        Mdouble massTotalParticles = particleHandler.getMass();
//
//        //Mdouble volParticlesPlusVoidInSystem = constants::pi*(std::pow(getXMax()/2,2.0))*getSurfaceHeight();
//        //Mdouble bulkDensity = massTotalParticles/volParticlesPlusVoidInSystem;
//        //Mdouble TheoDensity = 550.0;
//
//        os << getTime()*units::tUnit/1.0e-9 //[1] Current time
//           << " " << getTimeMax()*units::tUnit/1.0e-9 //[2] Max time
//           << " " << volSystem*pow(units::lUnit,3) //[3] Volume of the system
//           << " " << getZMax()*units::lUnit //[4] Height of the system
//           << " " << particleHandler.getNumberOfObjects() //[5] Number of particles inserted
//           << " " << particleHandler.getMeanRadius()*units::lUnit // [6] Mean radius of a particle
//           << " " << volParticle*pow(units::lUnit,3) //[7] Mean volume of a particle
//           << " " << densityParticle*units::rhoUnit() //[8] Mean density of a particle
//           << " " << volTotalParticles*pow(units::lUnit,3) //[9] Volume of particles inserted
//           << " " << volumeFraction //[10] Volume fraction
//           << " " << massTotalParticles*units::mUnit //[11] total mass of particle inserted
//           << " " << getSurfaceHeight()*units::lUnit // [12] Max height reached by particles.
//           << " " << getKineticEnergy()*units::energyUnit() // [13] kinetic energy
//           << " " << getElasticEnergy()*units::energyUnit() //[14] Elastic energy
//           << " " << scalarNormalForce*units::fUnit() //[15] Normal force
//           << " " << meanCoordinationNumber // [16] meanCoordinationNumber
//           << std::endl;
//    }
//
//    //--------------------------------------------------
//    // Print parameters of the system.
//    void printTime() const override
//    {
//
//        Mdouble volSystem = constants::pi*(std::pow(getXMax()/2,2))*getZMax();
//        Mdouble volTotalParticles = particleHandler.getVolume();
//        Mdouble volParticle = (4.0/3.0)*constants::pi*(std::pow(particleHandler.getMeanRadius(),3.0));
//        Mdouble densityParticle = particleSpecies_->getDensity();
//
//        Mdouble massParticle = densityParticle*volParticle;
//        Mdouble massTotalParticles = massParticle*particleHandler.getNumberOfObjects();
//        Mdouble massTotalParticles2 =  particleHandler.getMass();
//
//        Mdouble volParticlesPlusVoidInSystem = constants::pi*(std::pow(getXMax()/2/2,2.0))*newHeight;
//
//        Mdouble bulkDensity = massTotalParticles/(volParticlesPlusVoidInSystem);
//        Mdouble relativeDensity = bulkDensity/450.0;
//
//        Mdouble volumeFraction = volParticlesPlusVoidInSystem/volSystem;
//
//        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
//        << ", tmax=" << std::setprecision(3) << std::left << std::setw(3) << getTimeMax()
//        << ", No. Particles=" << std::setprecision(3) << std::left << std::setw(6) << particleHandler.getNumberOfObjects()
//        << ", MassTotalParticles=" << std::setprecision(3) << std::left<< std::setw(6)<< massTotalParticles2
//        << ", KineticEnergy=" << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()
//        << ", Volume inserted=" << particleHandler.getVolume()
//        << ", Volume Fraction=" << particleHandler.getVolume()/volSystem
//        << ", MaxHeightParticle pos=" << getSurfaceHeight()
//        << std::endl;
//    }
//};
//
//// Main function:
//int main(int argc UNUSED, char *argv[] UNUSED)
//{
//    //Set problem parameters:
//    std::string setFilename = "S1_PA12_500R_500L_Scaled";
//    Mdouble tMax = 1e-4/units::tUnit;
//    Mdouble gravity = -9.81/units::accUnit();
//
//    double RVESize = 60.0e-6/units::lUnit;
//    Mdouble meanRadius = 0.5*RVESize; //[mm]//
//
//    //define domain size
//    Mdouble domainLength = 2*RVESize; //[m]
//    Mdouble domainWidth = 2*RVESize; //[m]
//    Mdouble domainDepth = 2*RVESize; //[m]
//    Mdouble volumeFraction = 0.59;
//
//    //----------------------------------------------
//    //define properties of particle-particle contacts
//    Mdouble density = 1000.0/units::rhoUnit(); //[ug/mm^3]
//    Mdouble mass = density*(4.0/3.0)*pi*cubic(meanRadius); //mass of the particle
//
//    //--------------------------------------------------
//    // Setup parameters:
//    const Mdouble E_Modulus = 1650e12/units::youngUnit(); //Young's modulus
//    const Mdouble eta = 0.34; //Poisson ratio
//    const Mdouble E_reduced = 1.0/((1.0-eta*eta)/E_Modulus + (1.0-eta*eta)/E_Modulus);
//    const Mdouble r_effective = (meanRadius*meanRadius)/(meanRadius+meanRadius);
//    const Mdouble k1 = 4.0/3.0 * E_reduced*sqrt(r_effective);
//
//    const Mdouble pDepth = 1.45; //Penetration depth
//    const Mdouble k2 = 10.0*k1; //Unloading stiffness
//    const Mdouble kc = 0.0; //Cohesive stiffness
//    const Mdouble restitutionCoefficient = 0.1; //
//    //**********************************************************************
//    //----------------------------------------------
//    //create a contact law for particle-particle interactions
//    ThermalSinterLinFrictionReversibleAdhesiveSpecies particleSpecies;
//
//    //-------->[Set Plastic Contribution]
//    particleSpecies.setRestitutionCoefficient(restitutionCoefficient,mass);
//    particleSpecies.setPlasticParameters(k1, k2, kc, pDepth);
//    particleSpecies.setDensity(density);
////    particleSpecies.setStiffnessAndRestitutionCoefficient(k1, restitutionCoefficient, mass);
//
//    //Create a solver and run the commands in the constructor PowderBed::PowderBed
//    S1_Insertion pb (meanRadius,particleSpecies,volumeFraction);
//
//    Mdouble TimeStep = (particleSpecies.getCollisionTime(mass));
//    //----------------------------------------------
//    pb.setGravity({0,0,gravity}); //mm/mu s^2 [Take care]
//    pb.setMin({0,0,0});
//    pb.setMax({domainLength,domainWidth,domainDepth});
//    //----------------------------------------------
//
//    pb.setSaveCount(100);
//    pb.setTimeStep(TimeStep);
//    pb.setTimeMax(tMax);
//    //----------------------------------------------
//
//    logger(INFO,"Time step: %", pb.getTimeStep());
//    //----------------------------------------------
//    //    pb.wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);
//    pb.setParticlesWriteVTK(false);
//    pb.wallHandler.setWriteVTK(false);
//
//    pb.setName(setFilename);
//    pb.setXBallsAdditionalArguments("-solidf -v0");
//    pb.setFileType(FileType::ONE_FILE);
//
//    pb.removeOldFiles();
//    pb.solve();
//
//
//    //This helper is to plot results in gnuplot.
//    std::cout << "Execute gnuplot, 'plot S1_Insertion.fstat' u1:2 to view output" << std::endl;
//
//    //plot "S1_Insertion.fstat" u 5:10 with lines lw 1.5 lt 1 dt 5 lc "black" title "Cylinder-500L"
//
//    //plot "S1_Insertion.fstat" u 1:10 with linespoints ls 6 lc "black" pi 6 title "Cylinder"
//    //set yrange [0:0.8]
//    //set xlabel "time [s]"
//    //set ylabel "Volume Fraction"
//
//    //plot "S1_Insertion.fstat" u 5:10 with lines lw 1.5 lt 1 dt 5 lc "black" title "Cylinder-500L"
//    //replot "S1_Insertion_428L.fstat" u 5:10 with lines lw 1.5 lt 1 dt 3 lc "black" title "Cylinder-428L"
//    //replot "S1_Insertion_545L.fstat" u 5:10 with lines lw 1.5 lc "black" title "Cylinder-545L"
//
//    return 0;
//}
