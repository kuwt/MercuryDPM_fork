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

using constants::pi;
using mathsFunc::cubic;

class NeckGrowth : public Mercury3D
{
private:

    Mdouble radius_ = 0.0;
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
        ThermalParticle insertionBoundaryParticle;
        insertionBoundaryParticle.setSpecies(speciesHandler.getObject(0));

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
                P.setPosition(center + radius_ * pos*0.999);
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
           << " " << 1.0 //2
           << " " << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(0))->getRadius()//3
           << " " << particleSpecies_->getDensity()//4
           << " " << Force//5
           << " " << stress //6
           << " " << particleHandler.getMass()//7
           << " " << particleSpecies_->getLoadingStiffness()//8
           << " " << meanOverlap//9
           << " " << meanContactRadius //10
           << " " <<  particleHandler.getMeanRadius() //11
           << " " <<  getMeanRelativeContactRadius()//12 // Neckgrowth
           << " " << 1.0//13 [J]
           << " " << 1.0 //14 [W/s]
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
        std::cout << "t " << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << " tmax " <<  std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
                  << " square(delta/R) " << std::setprecision(3) << std::left << std::setw(6) << getMeanRelativeContactRadius()
                  << std::endl;
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    //Set problem parameters:
    std::string setFilename = "NeckGrowth_Polystyrene6um";

    Mdouble meanRadius = 6.0e-5; //[mm]//
    //define domain size
    Mdouble domainLength = 6*(2*meanRadius); //[mm]
    Mdouble domainWidth =  6*(2*meanRadius); //[mm]
    Mdouble domainDepth = 4*(2*meanRadius); //[mm]

    //----------------------------------------------
    //define properties of particle-particle contacts
    Mdouble density = 1000.0; //[mg/mm^3]
    Mdouble mass = density*(4.0/3.0)*pi*cubic(meanRadius); //mass of the particle

    //--------------------------------------------------
    // Setup parameters:
    const Mdouble E_Modulus = 3.250e9; //Young's modulus
    const Mdouble eta = 0.34; //Poisson ratio
    const Mdouble pDepth = 1.05; //Penetration depth
    const Mdouble E_reduced = 1.0/(2.0*((1.0 - eta*eta)/E_Modulus + (1.0-eta*eta)/E_Modulus));
    const Mdouble r_effective = (meanRadius*meanRadius)/(meanRadius+meanRadius);

    const Mdouble k1 =0.001*E_reduced*r_effective;
    const Mdouble k2 = 10.0*k1; //Unloading stiffness
    const Mdouble kc = 0.0; //Cohesive stiffness
    const Mdouble restitutionCoefficient = 0.1; //

    //----------------------------------------------
    Mdouble setDeltaC = 4.1e-07 ;
    Mdouble setC1 = 27.2;

    //create a contact law for particle-particle interactions
    ThermalSinterLinFrictionReversibleAdhesiveSpecies particleSpecies;
    //----------------------------------------------

    //-------->[Set Plastic Contribution]
    particleSpecies.setRestitutionCoefficient(restitutionCoefficient,mass);
    particleSpecies.setPlasticParameters(k1, k2, 0.0, pDepth);
    particleSpecies.setDensity(density);
    particleSpecies.setSinterAdhesion(0.0001*k1);

    //-------------------
    //Sinter parameters
    particleSpecies.setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);  //FRENKEL OR VISCOELASTIC_CONTACT
    particleSpecies.setComplianceZero(1/(2*E_Modulus)); //Book: Adhesive Particles
    particleSpecies.setSurfTension(0.035); // [N/m] or [J/m^2] surface energy density or surface tension.

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
    pb.setSaveCount(100);
    pb.setTimeStep(TimeStep);
    pb.setTimeMax(0.3);

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

    //This helper is to see the neck growth
    std::cout << "Execute 'gnuplot load 'NeckGrowth_Polystyrene.gnu' to view output" << std::endl;
    helpers::writeToFile(setFilename + ".gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'a(t)/R'\n"
                         "set xrange [0.0:0.3]\n"
                         "set yrange [0.0:1.3] \n"
                         "plot 'PS_R6e-5.txt' u ($1):($2) w p ls 1 lw 2.5 title 'Experimental data 30 um - PS'\n"
                         "a = 3.5\n"
                         "b = 0.125 \n"
                         "c = 2.62\n"
                         "f(x) = 3.2*sqrt(x)\n"
                         "g(x) = b * a * (x**(1.0/7.0))\n"
                         "theta(x) = c*sqrt(x)\n"
                         "h(x) = sin(theta(x))\n "
                         "t(x) = h(x) *(4.0/((1.0+cos(theta(x))**2.0)*(2.0-cos(theta(x)))))**(1.0/3.0)\n"
                         "replot 0.09 lc 'blue' t '1st mechanism'\n"
                         "replot g(x) lc 'green' t '2nd mechanism'\n"
                         "replot t(x) lc 'red' t '3rd mechanism Frenkel-Pokluda'\n"
                         "replot 'NeckGrowth_Polystyrene6um.fstat' u 1:12 with lines lw 1.5 lt 1 dt 5 lc \"black\" title 'DEM simulation' \n"

    );

    return 0;
}
//4/((1+cos(theta(x)))**2) * (2-cos(theta(x)))
//#include "DPMBase.h"
//#include <iostream>
//#include <vector>
//#include <Species/ThermalSinterLinFrictionReversibleAdhesiveSpecies.h>
//#include <Walls/InfiniteWall.h>
//#include <Logger.h>
//
/////This code tests our plastic force model, as published in Luding 2008.
//class Calibration_Polystyrene : public DPMBase {
//public:
//
//    explicit Calibration_Polystyrene(Mdouble radius) {
//        //-----------------
//        //Global parameters
//        std::string r = helpers::toString(radius);
//        setName("NeckGrowth_Polystyrene");
//
//        setFileType(FileType::ONE_FILE);
//        setParticlesWriteVTK(true);
//
//        //setParticlesWriteVTK(true);
//        //wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);
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

//
//        auto species = speciesHandler.copyAndAddObject(sf);
//
//        setTimeStep(collisionTime/50);
//        setSaveCount(1000);
//        setTimeMax(0.6);
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
//        P0.setPosition(Vec3D(-(1 - 1e-15) * radius, 0, 0));
//        P1.setPosition(-P0.getPosition());
//
//        particleHandler.copyAndAddObject(P0);
//        particleHandler.copyAndAddObject(P1);
//
//        //------------------
////        InfiniteWall w0;
////        w0.setSpecies(species);
////        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0, 0, -radius/4));
////        wallHandler.copyAndAddObject(w0);
////        const Mdouble scale = 1e12;
//
//    }
//
//    Mdouble getMeanRelativeContactRadius() const
//    {
//        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
//        Mdouble meanRadius = particleHandler.getMeanRadius();
//        logger(INFO," a/r %",sqrt(meanOverlap/meanRadius));
//        return sqrt(meanOverlap/meanRadius);
//
//    }
//
//    void writeFstatHeader(std::ostream &os) const override
//    {
////        Mdouble theta = atan(getMeanRelativeContactRadius());
//        Mdouble theta = getMeanRelativeContactRadius();
//        os << getTime() //1
//           << " " << theta//2
//
//        << std::endl;
//    }
//
//    void printTime() const override {
//
//        Mdouble theta = getMeanRelativeContactRadius();
////        Mdouble extra = std::pow(8.0/(pow((1.0+ cos(theta)),2.0) *(2.0-cos(theta))),1.0/3.0);
//
//        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
//                  << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
//                  << ", r=" << std::setprecision(3) << std::left << std::setw(6)
//                  << " Ene " << std::setw(12) << getKineticEnergy()/getElasticEnergy()
//                  << " X/a " << theta
//                  //<< ", Ekin=" << std::setprecision(3) << std::left << std::setw(6) << getKineticEnergy()/getElasticEnergy()
//                  << std::endl;
//        //std::cout.flush();
//    }
//};
//
//int main(int argc UNUSED, char *argv[] UNUSED)
//{
//    Calibration_Polystyrene sp0(60e-6);
//    sp0.solve();
////    calibrationSinterForceUnitTest sp1(1.5e-6);
////    sp1.solve();
////    calibrationSinterForceUnitTest sp2(1e-6);
////    sp2.solve();
//
//    //This helper is to see the Fn vs Overlap, rate of overlap.
////    std::cout << "Execute 'gnuplot NeckGrowth_Polystyrene.gnu' to view output" << std::endl;
////    helpers::writeToFile("NeckGrowth_Polystyrene.gnu",
////                         "set xlabel 'displacement'\n"
////                         "set ylabel 'force'\n"
////                         "plot 'NeckGrowth_Polystyrene.fstat' u 7:9 w lp\n"
////    );
//
//    //This helper is to see the neck growth
//    std::cout << "Execute 'gnuplot NeckGrowth_Polystyrene.gnu' to view output" << std::endl;
//    helpers::writeToFile("NeckGrowth_Polystyrene.gnu",
//                         "set xlabel 'time [s]'\n"
//                         "set ylabel 'a(t)/R'\n"
//                         "plot 'NeckGrowth_Polystyrene.fstat' u 1:2 title 'DEM simulation'  with lines linestyle 2\n"
//                         "replot 'PS_R6e-5.txt' u ($1):($2) w p ls 1 title 'Experimental data 60 um' \n"
//                         "replot 'PS_R60_E27.txt' u ($1):($2) w p ls 1 title 'Experimental data 60 um 27'\n"
//                         "replot 'PS_R60_E25.txt' u ($1):($2) w p ls 1 title 'Experimental data 60 um 25'\n"
//                         "replot 'PS_R60_E23.txt' u ($1):($2) w p ls 1 title 'Experimental data 60 um 23'"
//                         //"replot 0.37 lc 'red' "
////                         "replot 'SinterUnitTest1.5e-06.fstat' u ($1):(sqrt($7/1.5e-06)) w lp lt rgb 'light-red'\n"
////                         "replot 'SinterUnitTest1e-06.fstat' u ($1):(sqrt($7/1e-06)) w lp lt rgb 'light-green'"
//    );
//
//    return 0;
//
//}

//                         "plot 'NeckGrowth_Polystyrene.fstat' u ($1):(sin(atan(sqrt($7/30.5e-06)))) title 'DEM simulation'  with lines linestyle 2\n"