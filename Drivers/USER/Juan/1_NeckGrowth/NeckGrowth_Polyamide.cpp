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
#include <Walls/AxisymmetricIntersectionOfWalls.h>

using constants::pi;
using mathsFunc::cubic;

class NeckGrowth : public Mercury3D
{
private:

    Mdouble radius_ = 0.0;
    Mdouble deltaR = 1.0e-4; //[check] Thermal expansion coefficient

    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies_; //pointer to the particle properties

public:

    NeckGrowth (Mdouble radius, ParticleSpecies& particleSpecies)
    {
        radius_ = radius;
        //-------------------------------------------------
        speciesHandler.copyAndAddObject(particleSpecies); //particle-particle interactions
        //-------------------------------------------------
        particleSpecies_ = dynamic_cast<ThermalSinterLinFrictionReversibleAdhesiveSpecies*>(speciesHandler.getObject(0));
    }

    void setupInitialConditions() override {

        //++++++++++To store particle information++++++++++++++
        logger(INFO,"Adding particles ...");
        /* Introduce the InsertionBoundary */
//        ThermalParticle insertionBoundaryParticle;
//        insertionBoundaryParticle.setSpecies(speciesHandler.getObject(0));

        //--------------------------------------------------
        hGridRebuild();
        ThermalParticle P;
        P.setRadius(radius_);
        P.setSpecies(speciesHandler.getObject(0));

        unsigned nLayers = 1;
        std::array<Vec3D, 2> layerPos = {Vec3D(1, 0, 1), Vec3D(3, 0, 1)};

        double posFirst = 11 * radius_;
        for (unsigned i = 0; i < nLayers; i++) {
            Vec3D center = Vec3D(0.5 * (getXMax() - getXMin()),
                                 0.5 * (getYMax() - getYMin())+ i * 2 * radius_ + posFirst,
                                 0.5 * (getZMax() -getZMin()));
            for (auto pos : layerPos) {
                P.setPosition(center + radius_ * pos*0.9999 );
//                P.setTemperature(bedTemperature_);
//                P.setVelocity(Vec3D(0.0, 0.0, 0.0));
                particleHandler.copyAndAddObject(P);
                /// fix the last particles in space
//                if (i == nLayers - 1) {
//                    particleHandler.getLastObject()->fixParticle();
//                    /// give an initial kick to the last layer of particles
////                    particleHandler.getLastObject()->setPosition(P.getPosition() + (Vec3D(0.0, -0.05 * radius_, 0.0)));
//                }
            }
        }
//        hGridRebuild();
//        ThermalParticle P;
//        P.setRadius(radius_);
//        P.setSpecies(speciesHandler.getObject(0));
//
//        unsigned nLayers = 1;
//        std::array<Vec3D, 2> layerPos = {Vec3D(1, 0, 1), Vec3D(3, 0, 1)};
//
//        double posFirst = 11 * radius_;
//        for (unsigned i = 0; i < nLayers; i++) {
//            Vec3D center = Vec3D(0.5 * (getXMax() - getXMin()),
//                                 0.5 * (getYMax() - getYMin())+ i * 2 * radius_ + posFirst,
//                                 0.5 * (getZMax() -getZMin()));
//            for (auto pos : layerPos) {
//                P.setPosition(center + radius_ * pos*0.999);
////                P.setVelocity(Vec3D(0.0, 0.0, 0.0));
//                particleHandler.copyAndAddObject(P);
//                /// fix the last particles in space
////                if (i == nLayers - 1) {
////                    particleHandler.getLastObject()->fixParticle();
////                    /// give an initial kick to the last layer of particles
//////                    particleHandler.getLastObject()->setPosition(P.getPosition() + (Vec3D(0.0, -0.05 * radius_, 0.0)));
////                }
//            }
//        }

        std::cout << std::endl;
        logger(INFO,"Inserted % particles",particleHandler.getNumberOfObjects());
        std::cout << std::endl;
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
        Mdouble Force = dynamic_cast<const ThermalParticle*>(particleHandler.getObject(0))->getForce().X;
        Mdouble meanRadius = particleHandler.getMeanRadius();
        Mdouble meanContactRadius = sqrt(meanOverlap/meanRadius)*meanRadius;

        Mdouble K_s = 2.0*meanRadius/(std::pow(meanContactRadius,2.0));
        Mdouble Shrinkage = (1.0/5.0)*std::pow(getMeanRelativeContactRadius(),2.0);

        Mdouble stress = K_s*particleSpecies_->getSurfTension();

//        Mdouble rho_p = 1.0/(1 + std::pow(1 - X,3.0));

        os << getTime() //1
           << " " << 1.0 //2 getTemperature()
           << " " << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(0))->getRadius()//3
           << " " << particleSpecies_->getDensity()//4
           << " " << Force//5
           << " " << stress //6
           << " " << particleHandler.getMass()//7
           << " " << particleSpecies_->getLoadingStiffness()//8
           << " " << meanOverlap//9
           << " " << meanContactRadius //10
           << " " <<  particleHandler.getMeanRadius() //11
           << " " <<  1.0*getMeanRelativeContactRadius()//12 // Neckgrowth
           << " " << 1.0//13 [J]//getThermalEnergy()
           << " " << 1.0 //14 [W/s]getThermalEnergy()/(getTime()) //14 [W/s]
           << " " << getKineticEnergy() //15
           << " " << getElasticEnergy() //16
           << " " << K_s * particleSpecies_->getSurfTension() //17 Stress at the neckTip
           << " " << Shrinkage //18 shrinkage
           << " " << K_s//19
           << std::endl;
    }

    void printTime() const override
    {
        Mdouble volSystem = (getXMax()-getXMin())* (getYMax()-getYMin())*(getZMax()-getZMin());
        std::cout << "t: " << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << ", tmax: " <<  std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
                  << ", square(delta/R): " << std::setprecision(3) << std::left << std::setw(6) << getMeanRelativeContactRadius()
                  << ", Force_i: " << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(0))->getForce().X
                  << ", Force_j: " << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(1))->getForce().X
                  << std::endl;
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    //Set problem parameters:
    std::string setFilename = "NeckGrowthPA12";

    Mdouble meanRadius = 3.5e-5; //[mm]//
    //define domain size
    Mdouble domainLength = 6.0*(2*meanRadius); //[mm]
    Mdouble domainWidth =  4.0*(2*meanRadius); //[mm]
    Mdouble domainDepth = 4.0*(2*meanRadius); //[mm]

    //----------------------------------------------
    //define properties of particle-particle contacts
    Mdouble density = 930.0; //[mg/mm^3]
    Mdouble mass = density*(4.0/3.0)*pi*cubic(meanRadius); //mass of the particle

    //--------------------------------------------------
    // Setup parameters:
    const Mdouble E_Modulus = 1.650e9; //Young's modulus
    const Mdouble eta = 0.34; //Poisson ratio
    const Mdouble pDepth = 1.45; //Penetration depth
//    const Mdouble E_reduced = 1.0/(2.0*((1.0 - eta*eta)/E_Modulus + (1.0-eta*eta)/E_Modulus));
    const Mdouble r_effective = (meanRadius*meanRadius)/(meanRadius+meanRadius);

    const Mdouble E_reduced = 2.3e6; //[Pa]

    const Mdouble k1 = 0.001*E_reduced*r_effective;
    const Mdouble k2 = 10.0*k1; //Unloading stiffness
    const Mdouble kc = 0.0; //Cohesive stiffness
    const Mdouble restitutionCoefficient = 0.1; //

    //----------------------------------------------
    // sintering parameters
    //Virgin. deltaC = 2.8e-7, C1 = 1.71
    //1stAged: 2.9e-7, C1= 0.82
    //2ndAgeed: deltaC = 3.00e-7, C1 = 0.497
    Mdouble setDeltaC = 2.8e-7 ;
    Mdouble setC1 = 1.71;
    Mdouble compliance0 = 1.0/(2.0*E_Modulus);
    Mdouble surfaceTension = 0.041;

    //create a contact law for particle-particle interactions
    ThermalSinterLinFrictionReversibleAdhesiveSpecies particleSpecies;
    //----------------------------------------------

    //-------->[Set Plastic Contribution]
    particleSpecies.setRestitutionCoefficient(restitutionCoefficient,mass);
    particleSpecies.setPlasticParameters(k1, k2, kc, pDepth);
    particleSpecies.setDensity(density);
    particleSpecies.setSinterAdhesion(0.0016*k1);

    //-------------------
    //Sinter parameters
    particleSpecies.setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);  //FRENKEL OR VISCOELASTIC_CONTACT
    particleSpecies.setComplianceZero(compliance0); //Book: Adhesive Particles
    particleSpecies.setSurfTension(surfaceTension); // [N/m] or [J/m^2] surface energy density or surface tension.

    //Parameters to calibrate:
    particleSpecies.setFluidity(setC1);
    particleSpecies.setSeparationDis(setDeltaC);

    //----------------------------------------------
    //Create a solver and run the commands in the constructor PowderBed::PowderBed
    NeckGrowth pb (meanRadius,particleSpecies);

    //----------------------------------------------
    pb.setGravity({0,0,0}); //mm/mu s^2 [Take care]
    pb.setMin({0.0,0.0,0.0});
    pb.setMax({domainLength,domainWidth,domainDepth});
    //----------------------------------------------

    Mdouble TimeStep = particleSpecies.getCollisionTime(mass)/50;

    //----------------------------------------------
    pb.setSaveCount(10000);
    pb.setTimeStep(TimeStep);
    pb.setTimeMax(8.0);

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

    std::cout << "Execute gnuplot 'load 'NeckGrowthPA12.gnu' ' to view output" << std::endl;
    helpers::writeToFile(setFilename + ".gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'a(t)/R'\n"
                         "plot 'NeckGrowthPA12.fstat' u 1:12 title 'DEM simulation'  with lines linestyle 2\n"
                         //                         "replot 'PA12_Aged1_3.370e-5.txt' u ($1):($2) w p ls 1 title 'Experimental data 3.37 um' \n"
                         "replot 1.0 title 'a_{o} Limit' with lines linestyle 3 \n"
                         "replot 'PA12_2.584e-5.txt' u ($1):($2) w p ls 2 title 'Experimental data 2.584' \n"
                         "replot 'PA12_2.928e-5.txt' u ($1):($2) w p ls 3 title 'Experimental data 2.982 ' \n"
                         "replot 'PA12_3.044e-5.txt' u ($1):($2) w p ls 4 title 'Experimental data 3.044 '\n "
                         "replot 'PA12_3.210e-5.txt' u ($1):($2) w p ls 4 title 'Experimental data 3.210 ' \n"
                         "replot 'PA12_3.216e-5.txt' u ($1):($2) w p ls 4 title 'Experimental data 3.210 '\n "
                         "replot 'PA12_2.93e-5_Zhao.txt' u ($1):($2) w p ls 4 title 'Experimental data 2.93 Zao '\n "
//                         "replot '2AgedPA12_2.3e-05.txt' u ($1):($2) w p ls 4 title 'Experimental data First aged 2.330 '\n"
//                         "replot '2AgedPA12_3.8e-05.txt' u ($1):($2) w p ls 4 title 'Experimental data First aged 3.050 '\n"
//                         "replot '2AgedPA12_3.63e-05.txt' u ($1):($2) w p ls 4 title 'Experimental data First aged 3.37'\n "

    );


    std::cout << "Execute gnuplot 'load 'NeckGrowthPA122.gnu' ' to view output" << std::endl;
    helpers::writeToFile("NeckGrowthPA122.gnu",
                         "set xlabel 'a(t)/R'\n"
                         "set ylabel 'Q [J]'\n"
                         "plot 'NeckGrowthPA12.fstat' u 12:13 title 'Heat energy'  with lines linestyle 2");


    return 0;
}


//#include "DPMBase.h"
//#include <iostream>
//#include <vector>
//#include <Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h>
//#include <Walls/InfiniteWall.h>
//#include <Logger.h>
//
/////This code tests our plastic force model, as published in Luding 2008.
//class powdersAndgrains : public DPMBase {
//public:
//
//    explicit powdersAndgrains(Mdouble radius) {
//        //-----------------
//        //Global parameters
//        std::string r = helpers::to_string(radius);
//        setName("NeckGrowthPolyamide12_" + r);
//
//        setFileType(FileType::ONE_FILE);
////        setParticlesWriteVTK(true);
////        wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);
//
//        setParticleDimensions(3);
//        setSystemDimensions(3);
//        setGravity(Vec3D(0.0, 0.0,0.0));
//
//        //-------------------
//        //Boundary parameters
//        setXMax(2.0 * radius);
//        setYMax(radius);
//        setZMax(radius);
//        setXMin(-getXMax());
//        setYMin(-getYMax());
//        setZMin(-getZMax());
//
//        //-------------------
//        //Setup parameters:
//        const Mdouble density = 1000;
//
//        const Mdouble eta = 0.35; //Poisson ratio
//        const Mdouble E_Modulus = 1650e6; //[Pa] http://www.matweb.com/search/DataSheet.aspx?MatGUID=0e37a459c4eb452faa9d92659f9a0ccc&ckck=1
//
//        const Mdouble E_reduced = 1.0/(2.0*((1.0 - eta*eta)/E_Modulus + (1.0-eta*eta)/E_Modulus));
//        const Mdouble r_effective = (radius*radius)/(radius+radius);
//
//        const Mdouble k1 = E_reduced*r_effective;
////        const Mdouble k1 = 10.0 * radius; //[Stiffness depends on particle radius]
////        const Mdouble k1 = Stiffness;
//        const Mdouble restitutionCoefficient = 0.1;
//        const Mdouble pDepth = 1.45;
//
//        const Mdouble YoungM = 1.65e9; //[Pa] Young's Modulus for polyamide
//
//        //-------------------
//        //Species:
//        ThermalSinterLinFrictionReversibleAdhesiveSpecies sf;
//        //SinterSpecies sf;
//        sf.setHandler(&speciesHandler);
//        sf.setDensity(density);
//
//        const Mdouble mass = sf.getMassFromRadius(radius);
//        sf.setStiffnessAndRestitutionCoefficient(k1, restitutionCoefficient, mass);
//        sf.setPlasticParameters(k1, 2.0 * k1, 0.0, pDepth);
//
//        const Mdouble collisionTime = sf.getCollisionTime(mass);
//        logger(INFO, "Collision time %", collisionTime);
//
//        //-------------------
//        //Sinter parameters
//        sf.setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);  //FRENKEL OR VISCOELASTIC_CONTACT
//        sf.setSinterAdhesion(0.0016 * k1);
//
//        //-------------------
//        //Parameters to be used in contact growth driven by adhesive intersurface forces and accomodated by
//        //viscoelastic deformation.
//        sf.setComplianceZero(1/(2*YoungM)); //Book: Adhesive Particles
////        sf.setSurfTension(0.047); // [N/m] or [J/m^2] surface energy density or surface tension. Paper:Particle size limits for sintering polymer colloids. S Mazur.
//
//        //***********************************
//        // Control parameters:
////        sf.setSinterRate(1.7e-4);
//        sf.setFluidity(1.256); //$up$
//        sf.setSeparationDis(9.55e-07); //$down$
//        //***********************************
//
//        sf.setSurfTension(0.041); // [N/m] or [J/m^2] surface energy density or surface tension.
//        //http://www.surface-tension.de/solid-surface-energy.htm
//        auto species = speciesHandler.copyAndAddObject(sf);
//
//        setTimeMax(10);
//        setTimeStep(collisionTime/50.0);
//        setSaveCount(1000);
//
//        //-------------------
//        //Particle properties:
//        ThermalParticle P0, P1;
//        //SphericalParticle P0, P1;
//        P0.setSpecies(species);
//        P1.setSpecies(species);
//
//        P0.setRadius(radius);
//        P1.setRadius(radius);
//
//        P0.setPosition(Vec3D(-(1 - 1.0e-5) * radius, 0, 0));
//        P1.setPosition(-P0.getPosition());
//
//        particleHandler.copyAndAddObject(P0);
//        particleHandler.copyAndAddObject(P1);
//
//    }
//
//    Mdouble getMeanRelativeContactRadius() const
//    {
//        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
//        Mdouble meanRadius = particleHandler.getMeanRadius();
////        logger(INFO," a/r %",sqrt(meanOverlap/meanRadius));
//        return sqrt(meanOverlap/meanRadius);
//
//    }
//
//    void printTime() const override {
//        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
//        Mdouble meanRadius = particleHandler.getMeanRadius();
//
//        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
//                  << "t_step ="<< getTimeStep()
//                  << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
//                  << ", r=" << std::setprecision(3) << std::left << std::setw(6)
//                  << " Ene " << std::setw(12) << getKineticEnergy()/getElasticEnergy()
//                  //                  << " X/a " << std::setw(12) << getMeanRelativeContactRadius()
//                  << " square(delta/R) " << std::setw(12) << sqrt(meanOverlap/meanRadius)
//                  << std::endl;
//        //std::cout.flush();
//    }
//};
//
//int main(int argc UNUSED, char *argv[] UNUSED)
//{
//
//    std::ostringstream streamObj;
//
//    Mdouble materialC1 = 7.0;
//
//    std::string strObj = streamObj.str();
//
//    streamObj << materialC1;
//
//    std::cout << "materialC1";
//
//
//    powdersAndgrains sp0(3.39e-5);
//    sp0.solve();
//
//    //This helper is to see the Fn vs Overlap, rate of overlap.
//    std::cout << "Execute gnuplot 'load 'NeckGrowthPolyamide12.gnu' ' to view output" << std::endl;
//    helpers::writeToFile("NeckGrowthPolyamide12.gnu",
//                         "set xlabel 'displacement'\n"
//                         "set ylabel 'force'\n"
//                         "plot 'NeckGrowthPolyamide12_3.39e-05.fstat' u 7:9 w lp\n"
//                         "replot 'PA12_2.928e-5.txt' u ($1):($2) w p ls 3 title 'Experimental data 2.982' "
//    );
//
//    // time, i, j, x, y, z, delta, deltat, fn, ft, nx, ny, nz, tx, ty, tz
//
//    //This helper is to see the neck growth
//    std::cout << "Execute gnuplot 'load 'NeckGrowthPolyamide122.gnu' ' to view output" << std::endl;
//    helpers::writeToFile("PowdersAndGrains122.gnu",
//                         "set xlabel 'time [s]'\n"
//                         "set ylabel 'a(t)/R'\n"
//                         "plot 'NeckGrowthPolyamide12_3.39e-05.fstat' u ($1):(sqrt($7/3.39e-5)) title 'DEM simulation'  with lines linestyle 2\n"
//                         "replot 1.0 title 'a_{o} Limit' with lines linestyle 3 \n"
//                         "replot 'PA12_Aged1_3.39e-5.txt' ($1):($2) w p ls 1 title 'Experimental data'"
//
//    );
//
//    return 0;
//
//}