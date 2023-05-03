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

class Densification : public Mercury3D
{
private:

    Mdouble radius_ = 0.0;
    Mdouble deltaR_ = 0.0; //[check] Thermal expansion coefficient
    Mdouble numLayers_ = 0.0;

    Mdouble startingTemp_ = 0.0; //[C]
    Mdouble gradientTemp_ = 0.0; //[C/s] Heating rate
    Mdouble meltTemp_ = 0.0; //[C]
    Mdouble maxTemp_ = 0.0;

    Mdouble deltaTemperature = 0.0;

    Mdouble timeCooling_ = 0.0; //[s]


    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies_; //pointer to the particle properties
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpeciesAtMeltingPoint_;

public:

    Densification (Mdouble radius, Mdouble deltaR, Mdouble numLayers, ParticleSpecies& particleSpecies, Mdouble  startingTemp, Mdouble gradientTemp,
                   Mdouble meltTemp, Mdouble maxTemp, Mdouble timeCooling)
    {
        radius_ = radius;
        deltaR_ = deltaR;
        numLayers_ = numLayers;

        startingTemp_ =  startingTemp;
        gradientTemp_ = gradientTemp;
        meltTemp_ = meltTemp;
        maxTemp_ = maxTemp;
        timeCooling_ = timeCooling;

        //-------------------------------------------------
        speciesHandler.copyAndAddObject(particleSpecies); //particle-particle interactions
        //-------------------------------------------------
        particleSpecies_ = dynamic_cast<ThermalSinterLinFrictionReversibleAdhesiveSpecies*>(speciesHandler.getObject(0));
        particleSpeciesAtMeltingPoint_ = new ThermalSinterLinFrictionReversibleAdhesiveSpecies(*particleSpecies_);
    }

    void setupInitialConditions() override {

        particleSpeciesAtMeltingPoint_->setLoadingStiffness(particleSpecies_->getLoadingStiffness());

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

        unsigned nLayers = numLayers_;
        std::array<Vec3D, 9> layerPos = {Vec3D(1, 1, 0), Vec3D(1, -1, 0), Vec3D(-1, -1, 0), Vec3D(-1, 1, 0),
                                         Vec3D(3, 1, 0), Vec3D(3, -1, 0), Vec3D(3, -3, 0), Vec3D(1, -3, -0), Vec3D(-1, -3, 0)};
//        std::array<Vec3D, 36> layerPos = {Vec3D(1, 1, 0), Vec3D(1, -1, 0), Vec3D(1, -3, -0), Vec3D(-1, -3, 0),
//                                          Vec3D(1, -5, 0), Vec3D(-1, -5, 0),
//                                          Vec3D(-1, -1, 0), Vec3D(-1, 1, 0),
//                                          Vec3D(3, 1, 0), Vec3D(3, -1, 0),
//                                          Vec3D(3, -3, 0), Vec3D(3, -5, 0),
//                                          Vec3D(5, 1, 0), Vec3D(5, -1, -0),
//                                          Vec3D(5, -3, 0), Vec3D(5, -5, 0),
//                                          Vec3D(7, 1, 0), Vec3D(7, -1, 0),
//                                          Vec3D(7, -3, 0), Vec3D(7, -5, 0),
//                                          Vec3D(7, -7, 0), Vec3D(5, -7, 0),
//                                          Vec3D(3, -7, 0), Vec3D(1, -7, 0),
//                                          Vec3D(-1, -7, 0),Vec3D(9, 1, 0),
//                                          Vec3D(9, -1, 0),Vec3D(9, -3, 0),
//                                          Vec3D(9, -5, 0),Vec3D(9, -7, 0),
//                                          Vec3D(9, -9, 0),
//                                          Vec3D(7, -9, 0),Vec3D(5, -9, 0),
//                                          Vec3D(3, -9, 0),Vec3D(1, -9, 0),
//                                          Vec3D(-1, -9, 0)};

//        std::array<Vec3D, 36> layerPos = {Vec3D(1, 0, 1), Vec3D(1, 0, -1), Vec3D(1, 0, -3), Vec3D(-1, 0, -3),
//                                         Vec3D(1, 0, -5), Vec3D(-1, 0, -5),
//                                         Vec3D(-1, 0, -1), Vec3D(-1, 0, 1),
//                                         Vec3D(3, 0, 1), Vec3D(3, 0, -1),
//                                          Vec3D(3, 0, -3), Vec3D(3, 0, -5),
//                                          Vec3D(5, 0, 1), Vec3D(5, 0, -1),
//                                          Vec3D(5, 0, -3), Vec3D(5, 0, -5),
//                                          Vec3D(7, 0, 1), Vec3D(7, 0, -1),
//                                          Vec3D(7, 0, -3), Vec3D(7, 0, -5),
//                                          Vec3D(7, 0, -7), Vec3D(5, 0, -7),
//                                          Vec3D(3, 0, -7), Vec3D(1, 0, -7),
//                                          Vec3D(-1, 0, -7),Vec3D(9, 0, 1),
//                                          Vec3D(9, 0, -1),Vec3D(9, 0, -3),
//                                          Vec3D(9, 0, -5),Vec3D(9, 0, -7),
//                                          Vec3D(9, 0, -9),
//                                          Vec3D(7, 0, -9),Vec3D(5, 0, -9),
//                                          Vec3D(3, 0, -9),Vec3D(1, 0, -9),
//                                          Vec3D(-1, 0, -9)};



//        double posFirst = 11 * radius_;
//        for (unsigned i = 0; i < nLayers; i++) {
//            Vec3D center = Vec3D(0.5 * (getXMax() - getXMin()),
//                                 0.5 * (getYMax() - getYMin())+ i * 2 * radius_ + posFirst,
//                                 0.5 * (getZMax() -getZMin()));
//            for (auto pos : layerPos) {
//                P.setPosition(center + radius_ * pos*0.9999 );
//                P.setTemperature(bedTemperature_);
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

        double posFirst = radius_;
//        double posFirst = 1.0;
        for (unsigned i = 0; i < numLayers_; i++) {
            Vec3D center = Vec3D(getXMin()+ 2.0*radius_,
                                 getYMin()+ 2.0*radius_,
                                 getZMin()+ i * 2 * radius_ + posFirst);
            for (auto pos : layerPos) {
                P.setPosition(center + radius_ * pos*0.999);
//                P.fixParticle();
                P.setTemperature(startingTemp_);
                P.setVelocity(Vec3D(0.0, 0.0, 0.0));
                particleHandler.copyAndAddObject(P);
//                / fix the last particles in space
//                if (i == numLayers_ - 1) {
//                    particleHandler.getLastObject()->fixParticle();
//                    /// give an initial kick to the last layer of particles
////                    particleHandler.getLastObject()->setPosition(P.getPosition() + (Vec3D(0.0, -0.05 * radius_, 0.0)));
//                }

//                if(P.getPosition().Z <= 2* radius_){
//                    P.fixParticle();
//                }
            }
        }

        std::cout << std::endl;
        logger(INFO,"Inserted % particles",particleHandler.getNumberOfObjects());
        std::cout << std::endl;
    }

    //--------------------------------------------------
    Mdouble getTemperatureProfile()const
    {
        Mdouble heatingTemperature = startingTemp_ + gradientTemp_*getTime();
        Mdouble coolingTemperature = maxTemp_;

        if(getTime()>=timeCooling_){
            coolingTemperature = startingTemp_ + gradientTemp_*(getTimeMax() - getTime())*0.7;
        }
        return std::min(std::min(heatingTemperature,coolingTemperature),maxTemp_);
    }


    //Change particle temperature according to temperature profile
    void setTemperature(Mdouble tempStep)
    {
        for (BaseParticle* p : particleHandler)
        {
            ThermalParticle* tp = dynamic_cast<ThermalParticle*>(p);
            tp->setTemperature(tempStep);
        }

        static Mdouble oldTemperature = startingTemp_;
        deltaTemperature = tempStep-oldTemperature;
        Mdouble deltaCooling = tempStep - meltTemp_;

        Mdouble factorRadius = 1.0 - (deltaR_ * deltaTemperature);

        for (BaseParticle* p : particleHandler)
        {
            if (p->getSpecies()==particleSpecies_)
            {
                p->setRadius((factorRadius * p->getRadius()));
            }
        }
        //change species properties
        Mdouble stiffnessFactor = 0.5*(1.0+tanh((meltTemp_-tempStep)/gradientTemp_));
        Mdouble oldLoadingStiffness = particleSpecies_->getLoadingStiffness();
        particleSpecies_->setLoadingStiffness(stiffnessFactor*particleSpeciesAtMeltingPoint_->getLoadingStiffness());

        //for decreasing temperature, change the maxOverlap
        if (deltaCooling<0.0)
        {
            for (BaseInteraction* cBase : interactionHandler)
            {
                auto c =  dynamic_cast<SinterLinInteraction*>(cBase);
                Mdouble unloadingStiffness = c->getUnloadingStiffness();
                c->setMaxOverlap(c->getMaxOverlap()
                                 *(unloadingStiffness-oldLoadingStiffness)
                                 /(unloadingStiffness-particleSpecies_->getLoadingStiffness())
                );
            }
        }
        oldTemperature = tempStep;
    }

    //--------------------------------------------------
    void actionsAfterTimeStep() override
    {
        setTemperature(getTemperatureProfile()); //With Temperature profile. i.e., dilatometer

    }

    //--------------------------------------------------
    double getInfo(const BaseParticle& p) const override
    {
        return (dynamic_cast<const ThermalParticle &>(p).getTemperature()-startingTemp_)/(maxTemp_-startingTemp_);
    }

    //System
    Mdouble getThermalEnergy() const //ToDo: Check units.
    {
        Mdouble eneEnergy;

        for (BaseParticle* p : particleHandler)
        {
            if (p->isFixed()) continue;

            ThermalParticle* tp = dynamic_cast<ThermalParticle*>(p);
            eneEnergy += (tp->getTemperature()-startingTemp_) * tp->getMass() * particleSpecies_->getHeatCapacity();
//            eneEnergy += deltaTemperature * tp->getMass() * particleSpecies->getHeatCapacity();
        }
        return eneEnergy;
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
        Mdouble Force = abs(dynamic_cast<const ThermalParticle*>(particleHandler.getObject(0))->getForce().X);
        Mdouble meanRadius = particleHandler.getMeanRadius();
        Mdouble meanContactRadius = sqrt(meanOverlap/meanRadius)*meanRadius;
        Mdouble stress = Force/meanContactRadius;

        Mdouble X = getMeanRelativeContactRadius()/2;
        Mdouble rho_p = 1.0/(1 + std::pow(1 - X,3.0));

        os << getTime() //1
           << " " << getTemperatureProfile() //2
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
           << " " << getThermalEnergy()/2//13 [J]
           << " " << getThermalEnergy()/(getTime()) //14 [W/s]
           << " " << getKineticEnergy() //15
           << " " << getElasticEnergy() //16
           << " " << rho_p //17
           << std::endl;
    }

    void printTime() const override
    {
        Mdouble volSystem = (getXMax()-getXMin())* (getYMax()-getYMin())*(getZMax()-getZMin());
        std::cout << "t " << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << ", tmax " <<  std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
                  << ", Temperature " << std::setprecision(3) << std::left << std::setw(6) << getTemperatureProfile()
                  << ", pRadius" << std::setprecision(3) << std::left << std::setw(6) << particleHandler.getMeanRadius()
                  << ", K1" << std::setprecision(3) << std::left << std::setw(6) << particleSpecies_->getLoadingStiffness()
                  << ", Heat Energy " << std::setprecision(3) << std::left << std::setw(6) << getThermalEnergy()
                  << ", square(delta/R) " << std::setprecision(3) << std::left << std::setw(6) << getMeanRelativeContactRadius()
                  << std::endl;
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    //Set problem parameters:
    std::string setFilename = "DensificationTest";

    Mdouble meanRadius = 33.0e-6; //[mm]//
    Mdouble deltaR = 1.0e-4; //Relate to thermal expasion coefficient
    //define domain size
    Mdouble numLayers = 1;
    Mdouble domainLength = 3*(2*meanRadius); //[mm]
    Mdouble domainWidth = 3*(2*meanRadius); //[mm]
    Mdouble domainDepth = numLayers*(2*meanRadius); //[mm]

    //----------------------------------------------
    //define properties of particle-particle contacts
    Mdouble density = 930.0; //[mg/mm^3]
    Mdouble mass = density*(4.0/3.0)*pi*cubic(meanRadius); //mass of the particle

    //--------------------------------------------------
    // Setup parameters:
    const Mdouble E_Modulus = 1.650e9; //Young's modulus
    const Mdouble eta = 0.34; //Poisson ratio
    const Mdouble pDepth = 1.45; //Penetration depth
    const Mdouble E_reduced = 1.0/(2.0*((1.0 - eta*eta)/E_Modulus + (1.0-eta*eta)/E_Modulus));
    const Mdouble r_effective = (meanRadius*meanRadius)/(meanRadius+meanRadius);

    const Mdouble k1 =0.01*E_reduced*r_effective;
    const Mdouble k2 = 10.0*k1; //Unloading stiffness
    const Mdouble kc = 0.0; //Cohesive stiffness
    const Mdouble restitutionCoefficient = 0.1; //

    //----------------------------------------------
    Mdouble setDeltaC = 5.1e-07 ;
    Mdouble setC1 = 25.2;

    //Thermal variables:
    //PS
    Mdouble setHCapacity = 1185.0; //[J/(kgK)]
    Mdouble setTConductivity = 0.231;  //[W/(mK)]

    Mdouble startingTemp_ = 165.0;
    Mdouble gradientTemp = 240.0;
    Mdouble meltTemp = 180.0;
    Mdouble maxTemp = 181.0;
    Mdouble timeCooling = 0.18;

    //create a contact law for particle-particle interactions
    ThermalSinterLinFrictionReversibleAdhesiveSpecies particleSpecies;
    //----------------------------------------------

    //-------->[Set Plastic Contribution]
    particleSpecies.setRestitutionCoefficient(restitutionCoefficient,mass);
    particleSpecies.setPlasticParameters(k1, k2, kc, pDepth);
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

    particleSpecies.setThermalConductivity(setTConductivity);
    particleSpecies.setHeatCapacity(setHCapacity);

    //----------------------------------------------
    //Create a solver and run the commands in the constructor PowderBed::PowderBed
    Densification pb (meanRadius,deltaR, numLayers,particleSpecies,startingTemp_,gradientTemp,
                      meltTemp, maxTemp, timeCooling);

    //----------------------------------------------
    pb.setGravity({0,0,0}); //mm/mu s^2 [Take care]
    pb.setMin({0.0,0.0,0.0});
    pb.setMax({domainLength,domainWidth,domainDepth});
    //----------------------------------------------

    Mdouble TimeStep = particleSpecies.getCollisionTime(mass)/50;

    //----------------------------------------------
    pb.setSaveCount(1000);
    pb.setTimeStep(TimeStep);
    pb.setTimeMax(1.0);

    logger(INFO,"Time step: %", pb.getTimeStep());
    //----------------------------------------------
//    pb.wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);
    pb.setParticlesWriteVTK(false);

    pb.setName(setFilename);
    pb.setXBallsAdditionalArguments("-solidf -v0");
    pb.setFileType(FileType::ONE_FILE);

    pb.removeOldFiles();
    pb.setXBallsAdditionalArguments("-solidf -v0 -cmode 8 -cmaxset 100 ");

    pb.solve();

    std::cout << "Execute gnuplot 'load 'DensificationTest.gnu' ' to view output" << std::endl;
    helpers::writeToFile(setFilename + ".gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'a(t)/R'\n"
                         "plot 'DensificationTest.fstat' u 1:12 title 'DEM simulation'  with lines linestyle 2\n"
                         //                         "replot 'PA12_Aged1_3.370e-5.txt' u ($1):($2) w p ls 1 title 'Experimental data 3.37 um' \n"
                         //                         "replot 1.0 title 'a_{o} Limit' with lines linestyle 3 \n"
                         //                         "replot 'PA12_2.584e-5.txt' u ($1):($2) w p ls 2 title 'Experimental data 2.584' \n"
                         //                         "replot 'PA12_2.928e-5.txt' u ($1):($2) w p ls 3 title 'Experimental data 2.982 ' \n"
                         "replot 'PA12_3.044e-5.txt' u ($1):($2) w p ls 4 title 'Experimental data 3.044 '\n "
                         "replot 'PA12_3.210e-5.txt' u ($1):($2) w p ls 4 title 'Experimental data 3.210 ' \n"
                         "replot 'PA12_3.216e-5.txt' u ($1):($2) w p ls 4 title 'Experimental data 3.210 '\n "
//                         "replot 'PA12_2.93e-5_Zhao.txt' u ($1):($2) w p ls 4 title 'Experimental data 2.93 Zhao '\n "

    );

    std::cout << "Execute gnuplot 'load 'DensificationTest2.gnu' ' to view output" << std::endl;
    helpers::writeToFile(setFilename + "2.gnu",
                         "set xlabel 'a(t)/R'\n"
                         "set ylabel 'Q [J]'\n"
                         "plot 'DensificationTest.fstat' u 12:13 title 'Heat energy'  with lines linestyle 2"
    );

    return 0;
}