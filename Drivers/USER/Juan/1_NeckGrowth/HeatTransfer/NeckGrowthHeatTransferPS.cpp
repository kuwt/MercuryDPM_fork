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

using constants::pi;
using mathsFunc::cubic;

///This test models the neck growth kinetics + heat transfer of two polystyrene particles.
class NeckGrowthHeatTransfer : public Mercury3D {

private:

    Mdouble pRadius_ = 0.0;
    Mdouble deltaR_ = 0.0; //[check] Thermal expansion coefficient

    Mdouble startingTemp_ = 0.0; //[C]
    Mdouble gradientTemp_ = 0.0; //[C/s] Heating rate
    Mdouble meltTemp_ = 0.0; //[C]
    Mdouble maxTemp_ = 0.0;

    Mdouble deltaTemperature = 0.0;

    Mdouble timeCooling_ = 0.0; //[s]

//    unsigned index0_, index1_ ;
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies_; //pointer to the particle properties
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpeciesAtMeltingPoint;
//    Mdouble stress;
public:

    explicit NeckGrowthHeatTransfer(Mdouble radius, Mdouble deltaR, ParticleSpecies& particleSpecies, Mdouble  startingTemp, Mdouble gradientTemp,
                                    Mdouble meltTemp, Mdouble maxTemp, Mdouble timeCooling){

        //-----------------
        //Global parameters
        pRadius_ = radius;
        deltaR_ = deltaR;

        startingTemp_ =  startingTemp;
        gradientTemp_ = gradientTemp;
        meltTemp_ = meltTemp;
        maxTemp_ = maxTemp;
        timeCooling_ = timeCooling;

        //-------------------------------------------------
        speciesHandler.copyAndAddObject(particleSpecies); //particle-particle interactions
        //-------------------------------------------------
        particleSpecies_ = dynamic_cast<ThermalSinterLinFrictionReversibleAdhesiveSpecies*>(speciesHandler.getObject(0));
        particleSpeciesAtMeltingPoint = new ThermalSinterLinFrictionReversibleAdhesiveSpecies(*particleSpecies_);

        //-------------------
    }

    void setupInitialConditions() override
    {
        particleSpeciesAtMeltingPoint->setLoadingStiffness(particleSpecies_->getLoadingStiffness());
        //-------------------
//        //Particle properties:
//        ThermalParticle P0, P1;
//        //SphericalParticle P0, P1;
//        P0.setSpecies(particleSpecies_);
//        P1.setSpecies(particleSpecies_);
//
//        P0.setRadius(pRadius_);
//        P0.setTemperature(bedTemperature_);
//        P1.setRadius(pRadius_);
//        P1.setTemperature(bedTemperature_);
//
//        P0.setPosition(Vec3D(-(1 - 1.0e-5) * pRadius_, 0, 0));
//        P1.setPosition(-P0.getPosition());

//
//        //++++++++++To store particle information++++++++++++++
        logger(INFO,"Adding particles ...");
        /* Introduce the InsertionBoundary */
//        ThermalParticle insertionBoundaryParticle;
//        insertionBoundaryParticle.setSpecies(speciesHandler.getObject(0));

        //--------------------------------------------------
        hGridRebuild();
        ThermalParticle P;
        P.setRadius(pRadius_);
        P.setSpecies(speciesHandler.getObject(0));

        unsigned nLayers = 1;

        std::array<Vec3D, 2> layerPos = {Vec3D(1, 0, 1), Vec3D(3, 0, 1)};


        double posFirst = pRadius_;
//        double posFirst = 1.0;
        for (unsigned i = 0; i < nLayers; i++) {
            Vec3D center = Vec3D(getXMin()+ 2.0*pRadius_,
                                 getYMin()+ 2.0*pRadius_,
                                 getZMin()+ i * 1.999 * pRadius_ + posFirst);
            for (auto pos : layerPos) {
                P.setPosition(center + pRadius_ * pos*0.999);
//                P.fixParticle();
                P.setTemperature(startingTemp_);
//                P.setVelocity(Vec3D(0.0, 0.0, 0.0));
                particleHandler.copyAndAddObject(P);
            }
        }







//        std::array<Vec3D, 2> layerPos = {Vec3D(1, 0, 1), Vec3D(3, 0, 1)};
//
//        double posFirst = 11 * pRadius_;
//        for (unsigned i = 0; i < nLayers; i++) {
//            Vec3D center = Vec3D(0.5 * (getXMax() - getXMin()),
//                                 0.5 * (getYMax() - getYMin())+ i * 2 * pRadius_ + posFirst,
//                                 0.5 * (getZMax() -getZMin()));
//            for (auto pos : layerPos) {
//                P.setPosition(center + pRadius_ * pos*0.999);
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

//        particleHandler.copyAndAddObject(P0);
////        index0_ = P0.getIndex();
//        particleHandler.copyAndAddObject(P1);
//        index1_ = P1.getIndex();
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
    //--------------------------------------------------

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

        Mdouble factorRadius = 1.0 + (deltaR_ * deltaTemperature);

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
        particleSpecies_->setLoadingStiffness(stiffnessFactor*particleSpeciesAtMeltingPoint->getLoadingStiffness());

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
    Mdouble getMeanRelativeContactRadius() const
    {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanRadius = particleHandler.getMeanRadius();
        return sqrt(meanOverlap/meanRadius);
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

//    Function to override the output file with specific parameters
    void writeFstatHeader(std::ostream &os) const override
    {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble Force = dynamic_cast<const ThermalParticle*>(particleHandler.getObject(0))->getForce().X;
        Mdouble meanRadius = particleHandler.getMeanRadius();
        Mdouble meanContactRadius = sqrt(meanOverlap/meanRadius)*meanRadius;
        Mdouble stress = Force/meanContactRadius;

        os << getTime() //1
           << " " << getTemperatureProfile() //2
           << " " << particleHandler.getMeanRadius() // [6] Mean radius of a particle//3
           << " " << particleSpecies_->getDensity()//4
           << " " << Force//5
           << " " << stress //6
           << " " << particleHandler.getMass()//7
           << " " << particleSpecies_->getLoadingStiffness()//8
           << " " << meanOverlap//9
           << " " << meanContactRadius //10
           << " " <<  particleHandler.getMeanRadius() //11
           << " " <<  getMeanRelativeContactRadius()//12 // Neckgrowth
           << " " << getThermalEnergy()//13 [J]
           << " " << getThermalEnergy()/(getTime()) //14 [W/s]
           << " " << getKineticEnergy() //15
           << " " << getElasticEnergy() //16
           << std::endl;
    }

    void printTime() const override {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanRadius = particleHandler.getMeanRadius();

        std::cout << "t=" << getTime()
                  << ", tmax=" << getTimeMax()
                  << ", Ene " << getKineticEnergy()/getElasticEnergy()
                  << ", Temperature " << getTemperatureProfile()
                  << ", Heat Energy " << getThermalEnergy()
                  << ", Heat Energy/s [W] " << getThermalEnergy()/getTime()
                  << ", square(delta/R) " << sqrt(meanOverlap/meanRadius)
                  //                  << ", T_P0 " << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(index0_))->getTemperature()
                  //                  << ", T_P1 " << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(index1_))->getTemperature()
                  //                  << ", ThermalEnergy"<< getThermalEnergy()
                  //                  << ", MeanRadius" << particleHandler.getMeanRadius()
                  //                  << ", loading stiffness " << particleSpecies->getLoadingStiffness()
                  << std::endl;
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{

    std::string setFilename = "NeckGrowth_Polystyrene6um_29E";

    //Particle radius:
    Mdouble meanRadius = 3.0e-5; //[mm]//
    Mdouble deltaR = 8.0e-5; //[check] Thermal expansion coefficient

    Mdouble domainLength = 4*(2*meanRadius); //[mm]
    Mdouble domainWidth =  2*(2*meanRadius); //[mm]
    Mdouble domainDepth = 2*(2*meanRadius); //[mm]

    //----------------------------------------------
    //define properties of particle-particle contacts
    Mdouble density = 1000.0; //[mg/mm^3]
    Mdouble mass = density*(4.0/3.0)*pi*cubic(meanRadius); //mass of the particle

    //--------------------------------------------------
    // Setup parameters:
    const Mdouble E_Modulus = 3.250e9; //Young's modulus
    const Mdouble eta = 0.34; //Poisson ratio
    const Mdouble pDepth = 1.45; //Penetration depth
    const Mdouble E_reduced = 1.0/(2.0*((1.0 - eta*eta)/E_Modulus + (1.0-eta*eta)/E_Modulus));
    const Mdouble r_effective = (meanRadius*meanRadius)/(meanRadius+meanRadius);

    const Mdouble k1 =0.001*E_reduced*r_effective;
    const Mdouble k2 = 10.0*k1; //Unloading stiffness
    const Mdouble kc = 0.0; //Cohesive stiffness
    const Mdouble restitutionCoefficient = 0.1; //

    Mdouble setDeltaC = 1.1e-07 ;
    Mdouble setC1 = 0.8;
    Mdouble surfaceTension =0.04;

    //Thermal variables:
    Mdouble setHCapacity = 1320.0; //[J/(kgK)]
    Mdouble setTConductivity = 0.167;  //[W/(mK)]

    Mdouble startingTemp_ = 50.0;
    Mdouble gradientTemp = 240.0;
    Mdouble meltTemp = 72.0;
    Mdouble maxTemp = 73.0;
    Mdouble timeCooling = 0.18;

    //create a contact law for particle-particle interactions
    ThermalSinterLinFrictionReversibleAdhesiveSpecies particleSpecies;
    //----------------------------------------------

    //-------->[Set Plastic Contribution]
    particleSpecies.setRestitutionCoefficient(restitutionCoefficient,mass);
    particleSpecies.setPlasticParameters(k1, k2, 0.0, pDepth);
    particleSpecies.setDensity(density);
    particleSpecies.setSinterAdhesion(0.0001*k1);
    particleSpecies.setAdhesionForceMax(0.00000005*k1);

    //-------------------
    //Sinter parameters
    particleSpecies.setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);  //FRENKEL OR VISCOELASTIC_CONTACT
    particleSpecies.setComplianceZero(1/(2*E_Modulus)); //Book: Adhesive Particles
    particleSpecies.setSurfTension(surfaceTension); // [N/m] or [J/m^2] surface energy density or surface tension.

    //Parameters to calibrate:
    particleSpecies.setFluidity(setC1);
    particleSpecies.setSeparationDis(setDeltaC);

    particleSpecies.setThermalConductivity(setTConductivity);
    particleSpecies.setHeatCapacity(setHCapacity);


    NeckGrowthHeatTransfer test(meanRadius,deltaR,particleSpecies,startingTemp_,gradientTemp,
                                meltTemp, maxTemp, timeCooling);

    //----------------------------------------------
    test.setGravity({0,0,0}); //mm/mu s^2 [Take care]
    test.setMin({0.0,0.0,0.0});
    test.setMax({domainLength,domainWidth,domainDepth});
    //----------------------------------------------
    Mdouble TimeStep = particleSpecies.getCollisionTime(mass)/50;

    //----------------------------------------------
    test.setSaveCount(1000);
    test.setTimeStep(TimeStep);
    test.setTimeMax(0.45);

    logger(INFO,"Time step: %", test.getTimeStep());

    //----------------------------------------------
    //pb.wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);
    test.setParticlesWriteVTK(false);

    test.setName(setFilename);
    test.setXBallsAdditionalArguments("-solidf -v0 -cmode 8 -cmaxset 100 ");
    test.setFileType(FileType::ONE_FILE);

    test.removeOldFiles();
    test.solve();

    std::cout << "Execute gnuplot 'load 'NeckGrowth_Polystyrene6um_29E.gnu' ' to view output" << std::endl;
    helpers::writeToFile("NeckGrowth_Polystyrene6um_29E.gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'a(t)/R'\n"
                         "set xrange [0.0:0.45]\n"
                         "plot 'NeckGrowth_Polystyrene6um_29E.fstat' u 1:12 lc 'green 'title 'DEM simulation'\n"
                         //                         "replot 'NeckGrowth_Polystyrene6um_28E.fstat' u 1:12 title 'DEM simulation'  with lines linestyle 2\n"
                         //                         "replot 'NeckGrowth_Polystyrene6um_27E.fstat' u 1:12 title 'DEM simulation'  with lines linestyle 2\n"
                         "replot 1.0 title 'a_{o} Limit' with lines linestyle 3 \n"
                         //                         "replot 'PS_R60_E23.txt' u ($1):($2) w p ls 2 title 'Experimental data 60 um E23' \n"
                         //                         "replot 'PS_R60_E25.txt' u ($1):($2) w p ls 3 title 'Experimental data 60 um E25' \n"
                         "replot 'PS_R60_E27.txt' u ($1):($2) w p ls 2 lc 'red' title 'Experimental data 60 um E27' \n"
                         "replot 'PS_R60_E28.txt' u ($1):($2) w p ls 3 title 'Experimental data 60 um E28' \n"
                         "replot 'PS_R60_E29.txt' u ($1):($2) w p ls 3 title 'Experimental data 60 um E29'"
    );

    std::cout << "Execute gnuplot 'load 'NeckGrowth_Polystyrene6um_29E_2.gnu' ' to view output" << std::endl;
    helpers::writeToFile("NeckGrowth_Polystyrene6um_29E_2.gnu",
                         "set xlabel 'a(t)/R'\n"
                         "set ylabel 'Q [J]'\n"
                         "plot 'NeckGrowth_Polystyrene6um_29E.fstat' u 12:13 title 'Heat energy'  with lines linestyle 2"
    );
    return 0;

}