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
                P.setPosition(center + radius_ * pos*0.9999);
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
        std::cout << "t " << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << " tmax " <<  std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
                  << " square(delta/R) " << std::setprecision(3) << std::left << std::setw(6) << getMeanRelativeContactRadius()
                  << std::endl;
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    //Set problem parameters:
    std::string setFilename = "NeckGrowth_PEEK450PF";

    Mdouble meanRadius = 2.5e-5; //[mm]//
    //define domain size
    Mdouble domainLength = 6*(2*meanRadius); //[mm]
    Mdouble domainWidth =  4*(2*meanRadius); //[mm]
    Mdouble domainDepth = 4*(2*meanRadius); //[mm]

    //----------------------------------------------
    //define properties of particle-particle contacts
    Mdouble density = 1000.0; //[mg/mm^3]
    Mdouble mass = density*(4.0/3.0)*pi*cubic(meanRadius); //mass of the particle

    //--------------------------------------------------
    // Setup parameters:
    const Mdouble E_Modulus = 9.9e9; //Young's modulus
    const Mdouble eta = 0.3; //Poisson ratio
    const Mdouble pDepth = 1.45; //Penetration depth
    const Mdouble E_reduced = 1.0/(2.0*((1.0 - eta*eta)/E_Modulus + (1.0-eta*eta)/E_Modulus));
    const Mdouble r_effective = (meanRadius*meanRadius)/(meanRadius+meanRadius);

    const Mdouble k1 =0.0000001*E_reduced*r_effective;
    const Mdouble k2 = 10.0*k1; //Unloading stiffness
    const Mdouble kc = 0.0; //Cohesive stiffness
    const Mdouble restitutionCoefficient = 0.1; //

    //----------------------------------------------
    Mdouble setDeltaC = 3.65e-07 ;
    Mdouble setC1 = 0.275;

    //create a contact law for particle-particle interactions
    ThermalSinterLinFrictionReversibleAdhesiveSpecies particleSpecies;
    //----------------------------------------------

    //-------->[Set Plastic Contribution]
    particleSpecies.setRestitutionCoefficient(restitutionCoefficient,mass);
    particleSpecies.setPlasticParameters(k1, k2, kc, pDepth);
    particleSpecies.setDensity(density);
    particleSpecies.setSinterAdhesion(0.001*k1);

    //-------------------
    //Sinter parameters
    particleSpecies.setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);  //FRENKEL OR VISCOELASTIC_CONTACT
    particleSpecies.setComplianceZero(1.0/(2.0*E_Modulus)); //Book: Adhesive Particles
    particleSpecies.setSurfTension(0.03); // [N/m] or [J/m^2] surface energy density or surface tension.

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
    pb.setSaveCount(1000);
    pb.setTimeStep(TimeStep);
    pb.setTimeMax(80.0);

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
    std::cout << "Execute gnuplot 'load 'NeckGrowth_PEEK450PF.gnu' ' to view output" << std::endl;
    helpers::writeToFile(setFilename + ".gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'a(t)/R'\n"
                         "set title 'Experimental information from S. Berretta et. al'\n"
                         "plot 'NeckGrowth_PEEK450PF.fstat' u 1:12 title 'DEM simulation'  with lines linestyle 2\n"
                         "replot 1.0 title 'a_{o} Limit' with lines linestyle 3 \n"
                         "replot 'PEEK450PF.txt' u ($1):($2) w p ls 1 title 'Experimental data'\n"
//                         "a = 3.5\n"
//                         "b = 0.125 \n"
//                         "c = 2.62\n"
//                         "f(x) = 3.2*sqrt(x)\n"
//                         "g(x) = b * a * (x**(1.0/7.0))\n"
//                         "theta(x) = c*sqrt(x)\n"
//                         "h(x) = sin(theta(x))\n "
//                         "t(x) = h(x) *(4.0/((1.0+cos(theta(x))**2.0)*(2.0-cos(theta(x)))))**(1.0/3.0)\n"
//                         "replot 0.09 lc 'blue' t '1st mechanism'\n"
//                         "replot g(x) lc 'green' t '2nd mechanism'\n"
//                         "replot t(x) lc 'red' t '3rd mechanism Frenkel-Pokluda'\n"

    );

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
//class Calibration_PEEK450PF : public DPMBase {
//public:
//
//    explicit Calibration_PEEK450PF(Mdouble radius) {
//        //-----------------
//        //Global parameters
//        std::string r = helpers::to_string(radius);
//        setName("NeckGrowth_PEEK450_" + r);
//
//        setFileType(FileType::ONE_FILE);
//        //setParticlesWriteVTK(true);
//
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
//        //Setup parameters:
//        const Mdouble density = 1000;
//        const Mdouble k1 = 10.0 * radius; //[Stiffness depends on particle radius]
//        const Mdouble restitutionCoefficient = 0.1;
//        const Mdouble pDepth = 1.5;
//
//        const Mdouble YoungM = 1.9e9; //[Pa] Young's Modulus for PEEK450PF
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
//        sf.setPlasticParameters(k1, 10.0 * k1, k1, pDepth);
//
//        const Mdouble collisionTime = sf.getCollisionTime(mass);
//        logger(INFO, "Collision time %", collisionTime);
//
//        //-------------------
//        //Sinter parameters
//        sf.setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);  //FRENKEL OR VISCOELASTIC_CONTACT
//        sf.setSinterAdhesion(0.00016 * k1);
//        sf.setComplianceZero(1/(2*YoungM)); //Book: Adhesive Particles
//        sf.setSurfTension(0.04); // [N/m] or [J/m^2] surface energy density or surface tension. Paper:Particle size limits for sintering polymer colloids. S Mazur.
//
//        //-------------------
//        //Parameters to be used in contact growth driven by adhesive intersurface forces and accomodated by
//        //viscoelastic deformation.
//        // Control parameters:
//        sf.setFluidity(1.5e-1);
//        sf.setSeparationDis(9.8e-7);
//
//        auto species = speciesHandler.copyAndAddObject(sf);
//
//        setTimeMax(80);
//        setTimeStep(0.006 * collisionTime);
//        setSaveCount(100);
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
//    Calibration_PEEK450PF sp0(25.0e-6);
//    sp0.solve();
//
//    //This helper is to see the Fn vs Overlap, rate of overlap.
//    std::cout << "Execute gnuplot 'load 'NeckGrowth_PEEK450.gnu' ' to view output" << std::endl;
//    helpers::writeToFile("NeckGrowth_PEEK450.gnu",
//                         "set xlabel 'displacement'\n"
//                         "set ylabel 'force'\n"
//                         "plot 'NeckGrowth_PEEK450_2.5e-05.fstat' u 1:9 w lp\n"
//    );
//
//    // time, i, j, x, y, z, delta, deltat, fn, ft, nx, ny, nz, tx, ty, tz
//
//    //This helper is to see the neck growth
//    std::cout << "Execute gnuplot 'load 'NeckGrowth_PEEK450_2.gnu' ' to view output" << std::endl;
//    helpers::writeToFile("NeckGrowth_PEEK450_2.gnu",
//                         "set xlabel 'time [s]'\n"
//                         "set ylabel 'a(t)/R'\n"
//                         "set title 'Experimental information from S. Berretta et. al'\n"
//                         "plot 'NeckGrowth_PEEK450_2.5e-05.fstat' u ($1):(sqrt($7/50.0e-6)) title 'DEM simulation'  with lines linestyle 2\n"
//                         "replot 1.0 title 'a_{o} Limit' with lines linestyle 3 \n"
//                         "replot 'PEEK450PF.txt' ($1):($2) w p ls 1 title 'Experimental data'"
//
//    );
//
//    return 0;
//
//}