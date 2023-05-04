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

///This test models the neck growth kinetics + heat transfer of two polystyrene particles.
class NeckGrowthHeatTransfer : public Mercury3D {
public:

    explicit NeckGrowthHeatTransfer(Mdouble radius,Mdouble setC1,Mdouble setDeltaC,
            Mdouble setHCapacity, Mdouble setTConductivity) {

        //-----------------
        //Global parameters
        pRadius_ = radius;

        setFileType(FileType::ONE_FILE);
        setParticlesWriteVTK(true);
//        wallHandler.setWriteVTK(FileType::MULTIPLE_FILES);

        std::string r = helpers::to_string(pRadius_);
        setName("NeckGrowthHeatTransfer");
        setXBallsAdditionalArguments("-solidf -v0 -cmode 8 -cmaxset 100 ");
//        setXBallsAdditionalArguments("-v0 -solidf -p 1 -cmode 8");

        setGravity(Vec3D(0.0, 0.0,0.0));

        //-------------------
        //Boundary parameters
        setXMax(2.0 * pRadius_);
        setYMax(pRadius_);
        setZMax(pRadius_);
        setXMin(-getXMax());
        setYMin(-getYMax());
        setZMin(-getZMax());

        //-------------------
        //Setup parameters:
        // Polystyrene
        const Mdouble YoungM = 1.6e9; //[Pa] Young's Modulus
        const Mdouble density = 1040.0;
        const Mdouble pDepth = 1.26;
        const Mdouble restitutionCoefficient = 0.7;
//        const Mdouble k1 = 100.0 * pRadius_;//4*YoungM*pRadius_;//0.07;

        const Mdouble E_Modulus = 3250.0e6; //[Pa] http://www.matweb.com/search/DataSheet.aspx?MatGUID=0e37a459c4eb452faa9d92659f9a0ccc&ckck=1
        const Mdouble eta = 0.34; //Poisson ratio

        const Mdouble E_reduced = 1.0/(2.0*((1.0 - eta*eta)/E_Modulus + (1.0-eta*eta)/E_Modulus));
        const Mdouble r_effective = (pRadius_*pRadius_)/(pRadius_+pRadius_);

        const Mdouble k1 = 0.0001*E_reduced*r_effective;

        //--------------------------------------------------
        //Particle species
        particleSpecies = speciesHandler.copyAndAddObject(ThermalSinterLinFrictionReversibleAdhesiveSpecies());

        particleSpecies->setDensity(density);
        const Mdouble mass = particleSpecies->getMassFromRadius(pRadius_);
        particleSpecies->setRestitutionCoefficient(restitutionCoefficient,mass);
        particleSpecies->setPlasticParameters(k1, 10.0 * k1, 0.0, pDepth);
        particleSpecies->setSinterAdhesion(0.0001*k1);

        //-------------------
        //Sinter parameters
        particleSpecies->setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);  //FRENKEL OR VISCOELASTIC_CONTACT
        particleSpecies->setComplianceZero(1/(2*YoungM)); //Book: Adhesive Particles
        particleSpecies->setSurfTension(0.0356); // [N/m] or [J/m^2] surface energy density or surface tension.

        //Parameters to calibrate:
        particleSpecies->setFluidity(setC1);
        particleSpecies->setSeparationDis(setDeltaC);

        //Thermal parameters:
        particleSpecies->setThermalConductivity(setTConductivity);
        particleSpecies->setHeatCapacity(setHCapacity);

        //---Integration variables:
        const Mdouble collisionTime = particleSpecies->getCollisionTime(mass);
        logger(INFO, "Collision time %", collisionTime);
        setTimeStep(collisionTime/50.0);
        setSaveCount(1000);
        //--------------------------
        particleSpeciesAtMeltingPoint = new ThermalSinterLinFrictionReversibleAdhesiveSpecies(*particleSpecies);
    }

    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies; //pointer to the particle properties
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpeciesAtMeltingPoint;

    void setupInitialConditions() override
    {
        particleSpeciesAtMeltingPoint->setLoadingStiffness(particleSpecies->getLoadingStiffness());
        //-------------------
        //Particle properties:
        ThermalParticle P0, P1;
        //SphericalParticle P0, P1;
        P0.setSpecies(particleSpecies);
        P1.setSpecies(particleSpecies);

        P0.setRadius(pRadius_);
        P0.setTemperature(startingTemp);
        P1.setRadius(pRadius_);
        P1.setTemperature(startingTemp);

        P0.setPosition(Vec3D(-(1 - 1.0e-5) * pRadius_, 0, 0));
        P1.setPosition(-P0.getPosition());

        particleHandler.copyAndAddObject(P0);
        index0_ = P0.getIndex();
        particleHandler.copyAndAddObject(P1);
        index1_ = P1.getIndex();
    }
    //--------------------------------------------------
    Mdouble getTemperatureProfile()const
    {
        Mdouble heatingTemperature = startingTemp + gradientTemp*getTime();
        Mdouble coolingTemperature = maxTemp;

        if(getTime()>=timeCooling){
            coolingTemperature = startingTemp + gradientTemp*(getTimeMax() - getTime())*0.7;
        }
        return std::min(std::min(heatingTemperature,coolingTemperature),maxTemp);
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

        static Mdouble oldTemperature = startingTemp;
        Mdouble deltaTemperature = tempStep-oldTemperature;
        Mdouble deltaCooling = tempStep - meltTemp;

        Mdouble factorRadius = 1.0 - (deltaR * deltaTemperature);

        for (BaseParticle* p : particleHandler)
        {
            if (p->getSpecies()==particleSpecies)
            {
                p->setRadius((factorRadius * p->getRadius()));
            }
        }
        //change species properties
        Mdouble stiffnessFactor = 0.5*(1.0+tanh((meltTemp-tempStep)/gradientTemp));
        Mdouble oldLoadingStiffness = particleSpecies->getLoadingStiffness();
        particleSpecies->setLoadingStiffness(stiffnessFactor*particleSpeciesAtMeltingPoint->getLoadingStiffness());

        //for decreasing temperature, change the maxOverlap
        if (deltaCooling<0.0)
        {
            for (BaseInteraction* cBase : interactionHandler)
            {
                auto c =  dynamic_cast<SinterLinInteraction*>(cBase);
                Mdouble unloadingStiffness = c->getUnloadingStiffness();
                c->setMaxOverlap(c->getMaxOverlap()
                                 *(unloadingStiffness-oldLoadingStiffness)
                                 /(unloadingStiffness-particleSpecies->getLoadingStiffness())
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
        return (dynamic_cast<const ThermalParticle &>(p).getTemperature()-startingTemp)/(maxTemp-startingTemp);
    }

    Mdouble getThermalEnergy() const //ToDo: Check units.
    {
        Mdouble eneEnergy;

        for (BaseParticle* p : particleHandler)
        {
            if (p->isFixed()) continue;

            ThermalParticle* tp = dynamic_cast<ThermalParticle*>(p);
            eneEnergy += tp->getTemperature() * tp->getMass() * particleSpecies->getHeatCapacity();
        }
        return eneEnergy;
    }

//    Function to override the output file with specific parameters
    void writeFstatHeader(std::ostream &os) const override
    {
    Mdouble meanOverlap = interactionHandler.getMeanOverlap();
    Mdouble meanRadius = particleHandler.getMeanRadius();
    Mdouble meanContactRadius = sqrt(meanOverlap/meanRadius)*meanRadius;

        os << getTime() //1
           << " " << getTemperatureProfile() //2
           << " " << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(index0_))->getRadius()//3
           << " " << particleSpecies->getDensity()//4
           << " " << particleHandler.getMass()//5
           << " " << particleSpecies->getLoadingStiffness()//6
           << " " << meanOverlap//7
           << " " << meanContactRadius //8
           << " " <<  particleHandler.getMeanRadius() //9
           << " " <<  getMeanRelativeContactRadius()//10
           << " " << getThermalEnergy()//11
           << " " << getThermalEnergy()/getTime() //12
           << " " << ((getThermalEnergy()/getTime())/particleHandler.getMass())/1000.0 //13
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
                  << ", T_P0 " << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(index0_))->getTemperature()
                  << ", T_P1 " << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(index1_))->getTemperature()
                  << ", volume_P0 " << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(index0_))->getVolume()
//                  << ", ThermalEnergy"<< getThermalEnergy()
//                  << ", MeanRadius" << particleHandler.getMeanRadius()
//                  << ", loading stiffness " << particleSpecies->getLoadingStiffness()
                  << std::endl;
    }
private:

    unsigned index0_, index1_ ;

    Mdouble deltaR = 1.0e-4; //[check] Thermal expansion coefficient
    Mdouble gradientTemp = (170.0); //[C/min]
    Mdouble timeCooling = 0.3; //[s]

    const Mdouble startingTemp = 215.0; //[C]
    const Mdouble meltTemp = 240.0; //[C]
    const Mdouble maxTemp = 241.0;

    Mdouble sinterEnergy_ = 0.19e-6; //[J] toDo: Check the amount of energy provided by a laser
    Mdouble pRadius_ = 0.0;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    //Particle radius:
    Mdouble pRadius = 3.0e-5;
    Mdouble setC1 = 22.2;
    Mdouble setDeltaC = 5.1e-07 ;

    //Thermal variables:
    https://designerdata.nl/materials/plastics/thermo-plastics/polyamide-12
    Mdouble setHCapacity = 1320.0; //[J/(kgK)]
    Mdouble setTConductivity = 0.167;  //[W/(mK)]

    NeckGrowthHeatTransfer test(pRadius,setC1,setDeltaC,setHCapacity,setTConductivity);

    test.setTimeMax(0.3); //[s]
    test.removeOldFiles();
//    test.restartFile.setFileType(FileType::MULTIPLE_FILES);

    test.solve();

    //This helper is to see the Fn vs Overlap
//    std::cout << "Execute  gnuplot 'load 'NeckGrowthHeatTransfer1.gnu' to view output" << std::endl;
//    helpers::writeToFile("NeckGrowthHeatTransfer1.gnu",
//                         "set xlabel 'displacement'\n"
//                         "set ylabel 'force'\n"
//                         "plot 'NeckGrowthHeatTransfer.fstat' u 7:9 w lp\n"
//    );
//
//    std::cout << "Execute gnuplot 'load 'NeckGrowthHeatTransfer2.gnu' ' to view output" << std::endl;
//    helpers::writeToFile("NeckGrowthHeatTransfer2.gnu",
//                         "set xlabel 'time [s]'\n"
//                         "set ylabel 'a(t)/R'\n"
//                         "plot 'NeckGrowthHeatTransfer.fstat' u ($1):(sqrt($7/3.39e-5)) title 'DEM simulation'  with lines linestyle 2\n"
//                         "replot 1.0 title 'a_{o} Limit' with lines linestyle 3 \n"
//                         "replot 'PA12_Aged1_3.39e-5.txt' ($1):($2) w p ls 1 title 'Experimental data'"

//    );

    std::cout << "Execute gnuplot 'load 'NeckGrowthHeatTransfer.gnu' ' to view output" << std::endl;
    helpers::writeToFile("NeckGrowthHeatTransfer.gnu",
                         "set xlabel 'time [s]'\n"
                         "set ylabel 'a(t)/R'\n"
                         "plot 'NeckGrowthHeatTransfer.fstat' u 1:10 title 'DEM simulation'  with lines linestyle 2\n"
                         "replot 1.0 title 'a_{o} Limit' with lines linestyle 3 \n"
                         "replot 'PS_R3e-5.txt' ($1):($2) w p ls 1 title 'Experimental data'"
                         );

    std::cout << "Execute gnuplot 'load 'NeckGrowthHeatTransfer2.gnu' ' to view output" << std::endl;
    helpers::writeToFile("NeckGrowthHeatTransfer2.gnu",
                         "set xlabel 'a(t)/R'\n"
                         "set ylabel 'Q [J]'\n"
                         "plot 'NeckGrowthHeatTransfer.fstat' u 8:11 title 'Heat energy'  with lines linestyle 2"
                         );
    return 0;

}