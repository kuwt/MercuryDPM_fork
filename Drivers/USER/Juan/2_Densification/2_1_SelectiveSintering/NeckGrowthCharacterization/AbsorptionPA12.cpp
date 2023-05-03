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

using constants::pi;
using mathsFunc::cubic;

class NeckGrowth : public Mercury3D
{
private:

    Mdouble radius_ = 0.0;
    Mdouble deltaR_ = 0.0; //[check] Thermal expansion coefficient
    Mdouble heatTransferCoefficient_ = 0.0;

    Mdouble bedTemp_ = 0.0; //[C]

    Mdouble meltingTemp_ = 0.0; //[C]

    Mdouble deltaTemperature = 0.0;

    Mdouble holdingTime_ = 0.0; //[s]

    Mdouble scalarNormalForce = 0.0;

    Mdouble laserEnergy_ = 0.0;
    Mdouble laserPower = 0.0; //[J/s]
    Mdouble pulseDuration_ = 0.0; //[s]

    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies_; //pointer to the particle properties
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpeciesAtMeltingPoint_;

public:

    NeckGrowth (Mdouble radius, Mdouble deltaR, Mdouble heatTransferCoefficient, ParticleSpecies& particleSpecies, Mdouble  startingTemp,
                Mdouble meltTemp,Mdouble holdingTime, Mdouble laserEnergy, Mdouble pulseDuration)
    {
        radius_ = radius;
        deltaR_ = deltaR;
        heatTransferCoefficient_ = heatTransferCoefficient;

        bedTemp_ =  startingTemp;
        meltingTemp_ = meltTemp;
        holdingTime_ = holdingTime;

        laserEnergy_ = laserEnergy;
        pulseDuration_ = pulseDuration;

        //-------------------------------------------------
        speciesHandler.copyAndAddObject(particleSpecies); //particle-particle interactions
        //-------------------------------------------------
        particleSpecies_ = dynamic_cast<ThermalSinterLinFrictionReversibleAdhesiveSpecies*>(speciesHandler.getObject(0));
        particleSpeciesAtMeltingPoint_ = new ThermalSinterLinFrictionReversibleAdhesiveSpecies(*particleSpecies_);
    }

    void setupInitialConditions() override {

        particleSpeciesAtMeltingPoint_->setLoadingStiffness(particleSpecies_->getLoadingStiffness());
        particleSpeciesAtMeltingPoint_->setSinterAdhesion(particleSpecies_->getSinterAdhesion());
        particleSpeciesAtMeltingPoint_->setSurfTension(particleSpecies_->getSurfTension());
        particleSpeciesAtMeltingPoint_->setFluidity(particleSpecies_->getFluidity());

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
                P.setPosition(center + radius_ * pos*0.9999 );
                P.setTemperature(bedTemp_);
                particleHandler.copyAndAddObject(P);
            }
        }
        std::cout << std::endl;
        logger(INFO,"Inserted % particles",particleHandler.getNumberOfObjects());
        std::cout << std::endl;

        Mdouble massParticles = 0.5*particleHandler.getMass(); // [kg]
        Mdouble heatCapacity = particleSpecies_->getHeatCapacity(); //[J/KgK]

        laserPower = laserEnergy_/pulseDuration_;
    }

    Mdouble getHeatingRate() const
    {
        Mdouble energyAbsorbed = 0.0; //[%]

        //Coefficients of Absorption for R = 60.0e-6 m.
        Mdouble a1; Mdouble a2; Mdouble a3; Mdouble a4; Mdouble a5; Mdouble a6; Mdouble a7;

        if(getMeanRelativeContactRadius()<=0.312)
        {
            a1 = 4681.786;
            a2 = -4115.34240411;
            a3 = 1331.69895091;
            a4 = -193.62799263;
            a5 = 12.2508503;
            a6 = -0.31293957;
            a7 = 0.49883571;
        }else{
            a1 = 0.0;
            a2 = 0.0;
            a3 = 0.0;
            a4 = 0.0;
            a5 = 0.0;
            a6 = 0.0;
            a7 = 0.47833691;
        }

        energyAbsorbed = a1*pow(getMeanRelativeContactRadius(),6.0);
        energyAbsorbed += a2*pow(getMeanRelativeContactRadius(),5.0);
        energyAbsorbed += a3*pow(getMeanRelativeContactRadius(),4.0);
        energyAbsorbed += a4*pow(getMeanRelativeContactRadius(),3.0);
        energyAbsorbed += a5*pow(getMeanRelativeContactRadius(),2.0);
        energyAbsorbed += a6*getMeanRelativeContactRadius();
        energyAbsorbed += a7;

        return 4.0*(laserPower * energyAbsorbed) /(0.5*particleHandler.getMass()* particleSpecies_->getHeatCapacity()); //[k/s];;
    }
    //--------------------------------------------------
    //--------------------------------------------------
    Mdouble getTemperatureEvolution()const
    {

        Mdouble heatingTemperature = bedTemp_ + getHeatingRate()*getTime();

        Mdouble coolingTemperature = meltingTemp_;

        if(getTime() >= pulseDuration_){

            const static Mdouble maxTemp = heatingTemperature;

//            Mdouble m = (bedTemp_- maxTemp /holdingTime_);
//            coolingTemperature =  m*(getTime() - pulseDuration_) + maxTemp;
            double surfaceArea = 4.0*pi*std::pow(particleHandler.getMeanRadius(),2.0);
            double m = - ((heatTransferCoefficient_*surfaceArea*getTime())/(0.5*particleHandler.getMass()* particleSpecies_->getHeatCapacity()));
            coolingTemperature = (maxTemp-bedTemp_)*exp(m) + bedTemp_;

            if (coolingTemperature<bedTemp_) coolingTemperature = bedTemp_;
        }

//        }
        return std::min(heatingTemperature,coolingTemperature);
    }

    //--------------------------------------------------
    //Change particle temperature according to temperature profile
    void setTemperature(Mdouble tempStep)
    {
        //Variables to update
        static Mdouble oldTemperature = bedTemp_;
        Mdouble deltaTemperature = tempStep-oldTemperature;
        Mdouble deltaCooling = tempStep - meltingTemp_;

        for (BaseParticle* p : particleHandler)
        {
            auto *tp = dynamic_cast<ThermalParticle *>(p);

            if (p->getSpecies()==particleSpecies_)
            {
                tp->setTemperature(tempStep); //Change the temperature of the particle
            }
        }

        //Set changes of particles heated close to the melting poi
        setHeatEffect(tempStep,deltaTemperature,deltaCooling);

        oldTemperature = tempStep;
    }

    //--------------------------------------------------
    //Function to change particle properties based on the temperature
    void setHeatEffect(Mdouble tempStep,Mdouble deltaTemperature,Mdouble deltaCooling)
    {
        Mdouble factorRadius = 1.0 + (deltaR_ * deltaTemperature);

        for (BaseParticle* p : particleHandler)
        {
            if (p->getSpecies()==particleSpecies_)
            {
                auto* tp = dynamic_cast<ThermalParticle*>(p);
                Mdouble pT = tp->getTemperature();
                p->setRadius((factorRadius * p->getRadius()));

            }
        }
        //change species properties
        Mdouble stiffnessFactor = 0.5*(1.0+tanh((meltingTemp_-tempStep)/getHeatingRate()));

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
    }

    //--------------------------------------------------
    void actionsAfterTimeStep() override
    {
        setTemperature(getTemperatureEvolution());

        //To compute the normal force
        for (auto i : interactionHandler)
        {
            scalarNormalForce += Vec3D::dot(i->getForce(),i->getNormal());
        }

    }

    //--------------------------------------------------
    double getInfo(const BaseParticle& p) const override
    {
        return (dynamic_cast<const ThermalParticle &>(p).getTemperature()-bedTemp_)/(meltingTemp_-bedTemp_);  //        return getTemperatureEvolution();
    }

    //Thermal energy of the particle
    Mdouble getThermalEnergy() const
    {
        Mdouble eneEnergy;

        for (BaseParticle* p : particleHandler)
        {
            if (p->isFixed()) continue;

            ThermalParticle* tp = dynamic_cast<ThermalParticle*>(p);
            eneEnergy +=  (tp->getMass()/2.0) * particleSpecies_->getHeatCapacity()*(tp->getTemperature()-bedTemp_);
//            eneEnergy +=  (tp->getMass()) * particleSpecies_->getHeatCapacity()*(tp->getTemperature()-startingTemp_);

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
        Mdouble Force = scalarNormalForce;
        Mdouble meanRadius = particleHandler.getMeanRadius();
        Mdouble meanContactRadius = sqrt(meanOverlap/meanRadius)*meanRadius;
        Mdouble areaContact = pi* std::pow(meanContactRadius,2.0);
        Mdouble stress = Force/(areaContact);

        Mdouble X = getMeanRelativeContactRadius()/2;
        Mdouble rho_p = 1.0/(1 + std::pow(1 - X,3.0));

        Mdouble K_s = 2.0*meanRadius/(std::pow(meanContactRadius,2.0));
        Mdouble Shrinkage = (1.0/5.0)*std::pow(getMeanRelativeContactRadius(),2.0);

        Mdouble tangentialStress = K_s * particleSpecies_->getSurfTension();

        Mdouble rateContactGrowth = getMeanRelativeContactRadius()/getTime();

        //Characteristic time for heat diffusion tc for this system:
        //Laser sintering of PA12 particles studied by in-situ optical, thermal and X-ray characterization
        double tc = std::pow(2.0*radius_,2.0);
        tc *= particleSpecies_->getDensity()*particleSpecies_->getHeatCapacity();
        tc /= particleSpecies_->getThermalConductivity();

        os << getTime() //1
           << " " << getTemperatureEvolution() //2
           << " " << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(0))->getRadius()//3
           << " " << particleSpecies_->getDensity()//4
           << " " << Force//5
           << " " << stress //6
           << " " << particleHandler.getMass()//7
           << " " << particleSpecies_->getLoadingStiffness()//8
           << " " << meanOverlap//9
           << " " << meanContactRadius //10
           << " " << particleHandler.getMeanRadius() //11
           << " " << getMeanRelativeContactRadius()//12 // Neckgrowth
           << " " << getThermalEnergy()//13 [J]
           << " " << getThermalEnergy()/(getTime()) //14 [W]
           << " " << getKineticEnergy() //15
           << " " << getElasticEnergy() //16
           << " " << tangentialStress //17 Stress at the neckTip
           << " " << Shrinkage //18 shrinkage
           << " " << K_s//19
           << " " << rho_p //20
           << " " << scalarNormalForce //[21] Normal force
           << " " << particleSpecies_->getSinterAdhesion()// [22] sintering
           << " " << particleSpecies_->getFluidity()// [23] fluidity
           << " " << particleSpecies_->getSurfTension()// [24] Surface tension
           << " " << rateContactGrowth // [25] Rate of contact radius
           << " " << (getThermalEnergy())/(particleHandler.getMass()) // [26] J/kg
           << " " << tc //[27] characteristic heat diffusion time
           << " " << tangentialStress*(getTemperatureEvolution() - bedTemp_) //[28] Marangoni stress
           << " " << tangentialStress * (getTemperatureEvolution() - bedTemp_) * meanContactRadius //[29] Marangoni Force
           << std::endl;
    }

    void printTime() const override
    {
        Mdouble volSystem = (getXMax()-getXMin())* (getYMax()-getYMin())*(getZMax()-getZMin());
        std::cout << "t " << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << " tmax " << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
                  << " Temperature " << std::setprecision(3) << std::left << std::setw(6) << getTemperatureEvolution()
                  << " Fluidity" << std::setprecision(3) << std::left << std::setw(6) << particleSpecies_->getFluidity()
                  //                  << " pRadius" << std::setprecision(3) << std::left << std::setw(6) << particleHandler.getMeanRadius()
                  //                  << " Loading Stiffness" << std::setprecision(3) << std::left << std::setw(6) << particleSpecies_->getLoadingStiffness()
                  //                  << " Heat Energy " << std::setprecision(3) << std::left << std::setw(6) << getThermalEnergy()
                  << " square(delta/R) " << std::setprecision(3) << std::left << std::setw(6) << getMeanRelativeContactRadius()
                  << std::endl;
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    //File name
    std::string filename = "AbsorptionPA12_E192";

    //Particle information:
    Mdouble meanRadius = 125.0e-6; //[m]//
    Mdouble deltaRadius = 0.0001; // Todo: Thermal expansion coefficient (how large particle growth during the process)
    Mdouble heatTransferCoefficient = 200.0; //[W/(m2 K)]
    //Thermal variables:
    Mdouble setHCapacity = 1200.0; //[J/(kgK)]
    Mdouble setTConductivity = 0.231;  //[W/(mK)]

    //define domain size
    Mdouble domainLength = 6.0*(2*meanRadius); //[mm]
    Mdouble domainWidth =  4.0*(2*meanRadius); //[mm]
    Mdouble domainDepth = 4.0*(2*meanRadius); //[mm]

    //----------------------------------------------
    //define properties of particle
    //----------------------------------------------
    Mdouble density = 1250.0; //[kg/m^3]
    Mdouble mass = density*(4.0/3.0)*pi*cubic(meanRadius); //[kg] mass of the particle

    //--------------------------------------------------
    // Set-up parameters:
    const Mdouble E_Modulus = 1700e6; //[Pa] Young's modulus for PS particles
    const Mdouble eta = 0.34; //Poisson ratio
    const Mdouble pDepth = 1.35; //Penetration depth
    const Mdouble E_reduced = 1.0/((1.0-eta*eta)/E_Modulus + (1.0-eta*eta)/E_Modulus);
    const Mdouble r_effective = (meanRadius*meanRadius)/(meanRadius+meanRadius);
//    const Mdouble E_reduced = 2.3e6; //[Pa]

    Mdouble k1 = 4.0/3.0 * E_reduced*sqrt(r_effective)*1e-6;
    const Mdouble k2 = 5.0*k1; //Unloading stiffness
    const Mdouble kc = 0.0; //Cohesive stiffness
    const Mdouble restitutionCoefficient = 0.1; //

    Mdouble sinterAdhesion = 0.0001*k1;//Change while decreasing temperature

    //----------------------------------------------
    //Thermal parameters:
    Mdouble bedTemperature = 155.0; //[C] Initial temperature
    Mdouble meltingTemp = 178.0; //[C] Melting temperature of the material

    //LaserParameters:
    Mdouble laserEnergy = 2.0*192.0e-6; //[J]
    Mdouble pulseDuration = 0.001; //[s]
    Mdouble holdingTime = 0.9; // [s] This defines the cooling rate, change
    //----------------------------------------------

    // sintering parameters
    Mdouble setDeltaC = 0.6e-04 ; //[m] Adhesive traction in term of separation distance
    Mdouble setC1 = 24.3; // [1/Pa s] Fluidity

    Mdouble compliance0 = 1.0/(2.0*E_Modulus);
    Mdouble surfaceTension = 0.0323;

    //----------------------------------------------
    Mdouble simulationTime = 5.0;

    //create a contact law for particle-particle interactions
    ThermalSinterLinFrictionReversibleAdhesiveSpecies particleSpecies;
    //----------------------------------------------
    //-------->[Set Plastic Contribution]
    particleSpecies.setRestitutionCoefficient(restitutionCoefficient,mass);
    particleSpecies.setPlasticParameters(k1, k2, kc, pDepth);
    particleSpecies.setDensity(density);
    particleSpecies.setSinterAdhesion(sinterAdhesion);

    //-------------------
    //Sinter parameters
    particleSpecies.setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);  //FRENKEL OR VISCOELASTIC_CONTACT
    particleSpecies.setComplianceZero(compliance0); //Book: Adhesive Particles
    particleSpecies.setSurfTension(surfaceTension); // [N/m] or [J/m^2] surface energy density or surface tension.

    //Parameters to calibrate:
    particleSpecies.setFluidity(setC1);
    particleSpecies.setSeparationDis(setDeltaC);

    particleSpecies.setThermalConductivity(setTConductivity);
    particleSpecies.setHeatCapacity(setHCapacity);

    //----------------------------------------------
    //Create a solver and run the commands in the constructor PowderBed::PowderBed
    NeckGrowth pb (meanRadius, deltaRadius,heatTransferCoefficient, particleSpecies, bedTemperature,
                   meltingTemp, holdingTime, laserEnergy, pulseDuration);

    //----------------------------------------------
    pb.setGravity({0,0,0}); //mm/mu s^2 [Take care]
    pb.setMin({0.0,0.0,0.0});
    pb.setMax({domainLength,domainWidth,domainDepth});
    //----------------------------------------------

    Mdouble TimeStep = particleSpecies.getCollisionTime(mass)/50;

    //----------------------------------------------
    pb.setSaveCount(1000);
    pb.setTimeStep(TimeStep);
    pb.setTimeMax(simulationTime);

    logger(INFO,"Time step: %", pb.getTimeStep());
    //----------------------------------------------
//    pb.wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);
    pb.setParticlesWriteVTK(false);

    pb.setName(filename);
//    pb.setXBallsAdditionalArguments("-solidf -v0");
    pb.setFileType(FileType::ONE_FILE);

    pb.removeOldFiles();
    pb.setXBallsAdditionalArguments("-solidf -v0 -cmode 8 -cmaxset 100 ");

    pb.solve();

    //Neck growth kinetics
    std::cout << "Execute gnuplot 'load 'AbsorptionPA12.gnu' ' to view output" << std::endl;
    helpers::writeToFile("AbsorptionPA12.gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'a(t)/R'\n"
                         "set xrange [0:5.0]\n"
                         "plot 'AbsorptionPA12_E192.fstat' u 1:12 title 'E=192 $u$ J'\n"
                         "replot 1.0 title 'a_{o} Limit' with lines linestyle 3 \n"
                         "replot 'PA12_R125_E192.txt' u ($1):($2) w p ls 4 title 'Experimental data 3.210 ' \n"
    );

    //Thermal diffusion
    std::cout << "Execute gnuplot 'load 'AbsorptionPA12_HeatEnergyAndTime.gnu' ' to view output" << std::endl;
    helpers::writeToFile("AbsorptionPA12_HeatEnergyAndTime.gnu",
                         "set log x2 \n"
                         "set logscale xy \n"
                         "unset log x \n"
                         "unset log y \n"
                         "set xlabel 'a/R [-] '\n"
                         "set x2label 'Time for heat diffusion [s]'\n"
                         "set ylabel 'Q [J]'\n"
                         "set xtics nomirror\n"
                         "set x2tics\n"
                         "set tics out\n"
                         "set autoscale  x \n"
                         //                         "set autoscale x2 \n"
                         "plot 'AbsorptionPA12_E192.fstat' u 12:13 axes x1y1, 'AbsorptionPA12_E192.fstat' u (0.4):13 lw 0 axes x2y1");

    //Stress
    std::cout << "Execute gnuplot 'load 'AbsorptionPA12_Stress.gnu' ' to view output" << std::endl;
    helpers::writeToFile("AbsorptionPA12_Stress.gnu",
                         "set xlabel 'a(t)/R'\n"
                         "set ylabel 'Thermal stress at neck [Pa K]'\n"
                         "plot 'AbsorptionPA12_E192.fstat' u 2:28 title 'Thermal stress'  with lines linestyle 2");


    return 0;
}