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
//
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

    Mdouble bedTemperature_ = 0.0; //[C]
    Mdouble heatingRate_ = 0.0; //[C/s] Heating rate
    Mdouble meltingTemp_ = 0.0; //[C]
    Mdouble maxTemp_ = 0.0;

    Mdouble deltaTemperature = 0.0;

    Mdouble holdingTime_ = 0.0; //[s]

    Mdouble scalarNormalForce = 0.0;

    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies_; //pointer to the particle properties
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpeciesAtMeltingPoint_;

public:

    NeckGrowth (Mdouble radius, Mdouble deltaR, ParticleSpecies& particleSpecies, Mdouble  startingTemp, Mdouble gradientTemp,
                   Mdouble meltTemp, Mdouble maxTemp, Mdouble holdingTime)
    {
        radius_ = radius;
        deltaR_ = deltaR;

        bedTemperature_ =  startingTemp;
        heatingRate_ = gradientTemp;
        meltingTemp_ = meltTemp;
        maxTemp_ = maxTemp;
        holdingTime_ = holdingTime;

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
                P.setTemperature(bedTemperature_);
                particleHandler.copyAndAddObject(P);
            }
        }
        std::cout << std::endl;
        logger(INFO,"Inserted % particles",particleHandler.getNumberOfObjects());
        std::cout << std::endl;
    }

    //--------------------------------------------------
    Mdouble getTemperatureEvolution()const
    {
        Mdouble heatingTemperature = bedTemperature_ + heatingRate_*getTime();
        Mdouble coolingTemperature = maxTemp_;
        Mdouble isTimeForColling =  getTimeMax()- holdingTime_;

        if(heatingTemperature >= maxTemp_)
        {
            heatingTemperature = maxTemp_;

            if(getTime() >= isTimeForColling){
                coolingTemperature = bedTemperature_ + (maxTemp_-bedTemperature_)/(getTimeMax()-isTimeForColling)*(getTimeMax() - getTime());
            }
        }
        return std::min(std::min(heatingTemperature,coolingTemperature),maxTemp_);
    }

    //--------------------------------------------------
    //Change particle temperature according to temperature profile
    void setTemperature(Mdouble tempStep)
    {
        //Variables to update
        static Mdouble oldTemperature = bedTemperature_;
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
                //Check which particles are 95% closer to the melting point, then the particle radius can decrease
                if(pT >= 0.95*meltingTemp_)
                {
                    p->setRadius((factorRadius * p->getRadius()));
                }
            }
        }
        //change species properties
        Mdouble stiffnessFactor = 0.5*(1.0+tanh((meltingTemp_-tempStep)/heatingRate_));
        Mdouble adhesionFactor = 0.5*(1.0-tanh((meltingTemp_-tempStep)/heatingRate_));
        Mdouble complianceFactor = (2.0-tanh((maxTemp_-tempStep)/heatingRate_));
        Mdouble oldLoadingStiffness = particleSpecies_->getLoadingStiffness();

        particleSpecies_->setLoadingStiffness(stiffnessFactor*particleSpeciesAtMeltingPoint_->getLoadingStiffness());
        particleSpecies_->setSinterAdhesion(adhesionFactor*particleSpeciesAtMeltingPoint_->getSinterAdhesion());
        particleSpecies_->setFluidity(complianceFactor * particleSpeciesAtMeltingPoint_->getFluidity());
        particleSpecies_->setSurfTension(stiffnessFactor*particleSpeciesAtMeltingPoint_->getSurfTension());

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
        return (dynamic_cast<const ThermalParticle &>(p).getTemperature()-bedTemperature_)/(maxTemp_-bedTemperature_);  //        return getTemperatureEvolution();
    }

    //System
    Mdouble getThermalEnergy() const //ToDo: Check units.
    {
        Mdouble eneEnergy;

        for (BaseParticle* p : particleHandler)
        {
            if (p->isFixed()) continue;

            ThermalParticle* tp = dynamic_cast<ThermalParticle*>(p);
            eneEnergy += (tp->getTemperature()-bedTemperature_) * tp->getMass() * particleSpecies_->getHeatCapacity();
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
           << " " << getThermalEnergy()/(getTime()) //14 [W/s]
           << " " << getKineticEnergy() //15
           << " " << getElasticEnergy() //16
           << " " << tangentialStress //17 Stress at the neckTip
           << " " << Shrinkage //18 shrinkage
           << " " << K_s//19
           << " " << rho_p //20
           << " " << scalarNormalForce //[21] Normal force
           << " " << particleSpecies_->getSinterAdhesion()// [22] sintering
           << " " << particleSpecies_->getFluidity()// [23] fluidity
           << " " << particleSpecies_->getSurfTension()// [24] fluidity
           << " " << rateContactGrowth // [25] Rate of contact radius
           << std::endl;
    }

    void printTime() const override
    {
        Mdouble volSystem = (getXMax()-getXMin())* (getYMax()-getYMin())*(getZMax()-getZMin());
        std::cout << "t " << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << " tmax " << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
                  << " Temperature " << std::setprecision(3) << std::left << std::setw(6) << getTemperatureEvolution()
                  << " pRadius" << std::setprecision(3) << std::left << std::setw(6) << particleHandler.getMeanRadius()
                  << " Loading Stiffness" << std::setprecision(3) << std::left << std::setw(6) << particleSpecies_->getLoadingStiffness()
                  << " Heat Energy " << std::setprecision(3) << std::left << std::setw(6) << getThermalEnergy()
                  << " square(delta/R) " << std::setprecision(3) << std::left << std::setw(6) << getMeanRelativeContactRadius()
                  << std::endl;
    }
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    //Set problem parameters:
    std::string setFilename = "NeckGrowthHeatTransferPA12_T190";

    Mdouble meanRadius = 3.0e-5; //[mm]//
    Mdouble deltaRadius = 0.0002; //Thermal expansion coefficient 1/C
    //Thermal variables:
    Mdouble setHCapacity = 1185.0; //[J/(kgK)]
    Mdouble setTConductivity = 0.231;  //[W/(mK)]

    //define domain size
    Mdouble domainLength = 6*(2*meanRadius); //[mm]
    Mdouble domainWidth =  4*(2*meanRadius); //[mm]
    Mdouble domainDepth = 4*(2*meanRadius); //[mm]

    //----------------------------------------------
    //define properties of particle-particle contacts
    Mdouble density = 930.0; //[mg/mm^3]
    Mdouble mass = density*(4.0/3.0)*pi*cubic(meanRadius); //mass of the particle

    //--------------------------------------------------
    // Setup parameters:
    const Mdouble E_Modulus = 1700.0e6; //Young's modulus
    const Mdouble eta = 0.34; //Poisson ratio
    const Mdouble pDepth = 1.45; //Penetration depth
    const Mdouble E_reduced = 1.0/((1.0-eta*eta)/E_Modulus + (1.0-eta*eta)/E_Modulus);
    const Mdouble r_effective = (meanRadius*meanRadius)/(meanRadius+meanRadius);
//    const Mdouble E_reduced = 2.3e6; //[Pa]

    Mdouble k1 = 4.0/3.0 * E_reduced*sqrt(r_effective)*1e-8;
    const Mdouble k2 = 5.0*k1; //Unloading stiffness
    const Mdouble kc = 0.0; //Cohesive stiffness
    const Mdouble restitutionCoefficient = 0.1; //

    Mdouble sinterAdhesion = 0.0001*k1;
    //----------------------------------------------
    // define temperature profile
    Mdouble startingTemp_ = 160.0;
    Mdouble maxTemp = 190.0;
    Mdouble gradientTemp = 10.0;//[C/s]
    Mdouble meltingTemp = 175.0;

    Mdouble holdindingTime = 1.0;
    Mdouble maxTime = 10.0;

    // sintering parameters
    Mdouble setDeltaC = 2.8e-07 ;
    Mdouble setC1 = 1.7;
    Mdouble compliance0 = 1.0/(2.0*E_Modulus);
    Mdouble surfaceTension = 0.04;

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
    NeckGrowth pb (meanRadius,deltaRadius,particleSpecies,startingTemp_,gradientTemp,
                      meltingTemp, maxTemp, holdindingTime);

    //----------------------------------------------
    pb.setGravity({0,0,0}); //mm/mu s^2 [Take care]
    pb.setMin({0.0,0.0,0.0});
    pb.setMax({domainLength,domainWidth,domainDepth});
    //----------------------------------------------

    Mdouble TimeStep = particleSpecies.getCollisionTime(mass)/50;

    //----------------------------------------------
    pb.setSaveCount(1000);
    pb.setTimeStep(TimeStep);
    pb.setTimeMax(maxTime);

    logger(INFO,"Time step: %", pb.getTimeStep());
    //----------------------------------------------
//    pb.wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);
    pb.setParticlesWriteVTK(false);

    pb.setName(setFilename);
//    pb.setXBallsAdditionalArguments("-solidf -v0");
    pb.setFileType(FileType::ONE_FILE);

    pb.removeOldFiles();
    pb.setXBallsAdditionalArguments("-solidf -v0 -cmode 8 -cmaxset 100 ");

    pb.solve();

    std::cout << "Execute gnuplot 'load 'NeckGrowthHeatTransferPA12.gnu' ' to view output" << std::endl;
    helpers::writeToFile("NeckGrowthHeatTransferPA12.gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'a(t)/R'\n"
                         "set xrange [0:8.0]\n"
//                         "plot 'NeckGrowthHeatTransferPA12.fstat' u 1:12 title 'DEM simulation - Tmelting' \n"
                         "plot 'NeckGrowthHeatTransferPA12_T160.fstat' u 1:12 title 'DEM simulation - Tmelting - (160Tmax)'\n"
                         "replot 'NeckGrowthHeatTransferPA12_T170.fstat' u 1:12 title 'DEM simulation - Tmelting - (170Tmax)'\n"
                         "replot 'NeckGrowthHeatTransferPA12_T180.fstat' u 1:12 title 'DEM simulation - Tmelting - (180Tmax)'\n"
                         "replot 'NeckGrowthHeatTransferPA12_T190.fstat' u 1:12 title 'DEM simulation - Tmelting - (190Tmax)'\n"

//                         "replot 'NeckGrowthHeatTransferPA12_6.fstat' u 1:12 title 'DEM simulation - Tmelting - (185Tmax)'\n"
//                         "replot 'NeckGrowthHeatTransferPA12_7.fstat' u 1:12 title 'DEM simulation - Tmelting - (175Tmax)'\n"

//                         "replot 'PA12_Aged1_3.370e-5.txt' u ($1):($2) w p ls 1 title 'Experimental data 3.37 um' \n"
//                         "replot 1.0 title 'a_{o} Limit' with lines linestyle 3 \n"
//                         "replot 'PA12_2.584e-5.txt' u ($1):($2) w p ls 2 title 'Experimental data 2.584' \n"
//                         "replot 'PA12_2.928e-5.txt' u ($1):($2) w p ls 3 title 'Experimental data 2.982 ' \n"
//                         "replot 'PA12_3.044e-5.txt' u ($1):($2) w p ls 4 title 'Experimental data 3.044 '\n "
//                         "replot 'PA12_3.210e-5.txt' u ($1):($2) w p ls 4 title 'Experimental data 3.210 ' \n"
//                         "replot 'PA12_3.216e-5.txt' u ($1):($2) w p ls 4 title 'Experimental data 3.210 '\n "
                         "replot 'PA12_2.93e-5_Zhao.txt' u ($1):($2) w p ls 4 title 'Experimental data 2.93 Zhao '\n "
//                         "replot 'PA12_Aged1_2.330e-5.txt' u ($1):($2) w p ls 4 title 'Experimental data First aged 2.330 '\n"
//                         "replot 'PA12_Aged1_3.050e-5.txt' u ($1):($2) w p ls 4 title 'Experimental data First aged 3.050 '\n"

    );


    std::cout << "Execute gnuplot 'load 'NeckGrowthHeatTransferPA122.gnu' ' to view output" << std::endl;
    helpers::writeToFile("NeckGrowthHeatTransferPA122.gnu",
                         "set xlabel 'a(t)/R'\n"
                         "set ylabel 'Q [J]'\n"
                         "plot 'NeckGrowthHeatTransferPA12.fstat' u 12:13 title 'Heat energy'  with lines linestyle 2");


    return 0;
}