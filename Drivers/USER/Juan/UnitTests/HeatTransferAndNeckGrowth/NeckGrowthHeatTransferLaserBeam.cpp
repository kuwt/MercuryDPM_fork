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

///This code tests our plastic force model, as published in Luding 2008.
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

        std::string r = helpers::toString(pRadius_);
        setName("NeckGrowthHeatTransfer2");
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
        const Mdouble YoungM = 1.65e9; //[Pa] Young's Modulus for polyamide
        const Mdouble density = 930.0;
        const Mdouble pDepth = 1.35;
        const Mdouble restitutionCoefficient = 0.7;
        const Mdouble k1 = 10.0 * pRadius_;//4*YoungM*pRadius_;//0.07;

        //--------------------------------------------------
        //Particle species
        particleSpecies = speciesHandler.copyAndAddObject(ThermalSinterLinFrictionReversibleAdhesiveSpecies());

        particleSpecies->setDensity(density);
        const Mdouble mass = particleSpecies->getMassFromRadius(pRadius_);
        particleSpecies->setRestitutionCoefficient(restitutionCoefficient,mass);
        particleSpecies->setPlasticParameters(k1, 2.0 * k1, 0.0, pDepth);

        //-------------------
        //Sinter parameters
        particleSpecies->setSinterType(SINTER_APPROACH::VISCOELASTIC_CONTACT);  //FRENKEL OR VISCOELASTIC_CONTACT
        particleSpecies->setComplianceZero(1/(2*YoungM)); //Book: Adhesive Particles
        particleSpecies->setSurfTension(0.047); // [N/m] or [J/m^2] surface energy density or surface tension.
        particleSpecies->setSinterAdhesion(0.001*k1); //ToDo: It acts as a sintering force

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
        setSaveCount(100);
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
            coolingTemperature = startingTemp + gradientTemp*(getTimeMax() - getTime());
        }
        return std::min(std::min(heatingTemperature,coolingTemperature),maxTemp);
    }
    //--------------------------------------------------
    void getTempFromLaser()
    {
        //Count particles
        Mdouble massHeatedParticles = 0;
        Mdouble tempLaser = 0.0;
        for (BaseParticle* p : particleHandler)
        {
            massHeatedParticles += p->getMass();
        }

        //heat particles
        for (BaseParticle* p : particleHandler)
        {
            ThermalParticle* tp = dynamic_cast<ThermalParticle*>(p);
            tempLaser = sinterEnergy_ / massHeatedParticles / particleSpecies->getHeatCapacity();
            tp->addTemperature(tempLaser);
        }
        tempStep_ = tempLaser - 273.0;
    }

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
//        getTempFromLaser(); // 2_1_LaserBeam
//        setTemperature(tempStep_);

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
        Mdouble eneThermal = 0;
        for (BaseParticle* p : particleHandler)
        {
            if (p->isFixed()) continue;
            ThermalParticle* tp = dynamic_cast<ThermalParticle*>(p);
            eneThermal += tp->getTemperature() * tp->getMass() * particleSpecies->getHeatCapacity();
        }
        return eneThermal;
    }

    //Function to override the output file with specific parameters
    void writeFstatHeader(std::ostream &os) const override
    {
    Mdouble meanOverlap = interactionHandler.getMeanOverlap();
    Mdouble meanRadius = particleHandler.getMeanRadius();
    Mdouble meanContactRadius = sqrt(meanOverlap/meanRadius)*meanRadius;

        os << getTime() //1
           << " " << getTemperatureProfile() //2
           << " " << particleSpecies->getDensity()//3
           << " " << particleSpecies->getLoadingStiffness()//4
           << " " << meanOverlap//5
           << " " << meanContactRadius //6
           << " " <<  getMeanRelativeContactRadius()//7
           << " " << getThermalEnergy() //8
           << std::endl;
    }

    void printTime() const override {
        Mdouble meanOverlap = interactionHandler.getMeanOverlap();
        Mdouble meanRadius = particleHandler.getMeanRadius();

        std::cout << "t=" << getTime()
                  << ", tmax=" << getTimeMax()
                  << ", Ene " << getKineticEnergy()/getElasticEnergy()
                  << ", square(delta/R) " << sqrt(meanOverlap/meanRadius)
                  << ", T_P0 " << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(index0_))->getTemperature()
                  << ", T_P1 " << dynamic_cast<const ThermalParticle*>(particleHandler.getObject(index1_))->getTemperature()
                  << ", ThermalEnergy"<< getThermalEnergy()
                  << ", MeanRadius" << particleHandler.getMeanRadius()
                  << ", loading stiffness " << particleSpecies->getLoadingStiffness()
                  << std::endl;
    }
private:

    Mdouble tempStep_ = 0.0;
    unsigned index0_, index1_ ;

    Mdouble deltaR = 1.0e-4; //[check] Thermal expansion coefficient
    Mdouble gradientTemp = 7.0; //[C/s]
    Mdouble timeCooling = 0.05; //[s]

    const Mdouble startingTemp = 170.0; //[C]
    const Mdouble meltTemp = 175.0; //[C]
    const Mdouble maxTemp = 176.0;

    Mdouble sinterEnergy_ = 0.16e-3; //[J] toDo: Check the amount of energy provided by a laser
    Mdouble pRadius_ = 0.0;
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    //Particle radius:
    Mdouble pRadius = 3.0e-5;
    Mdouble setC1 = 0.4;
    Mdouble setDeltaC = 2.7e-07 ;

    //Thermal variables:
    https://designerdata.nl/materials/plastics/thermo-plastics/polyamide-12
    Mdouble setHCapacity = 1185.0; //[J/kg/K]
    Mdouble setTConductivity = 0.231;  //W/m/K

    NeckGrowthHeatTransfer test(pRadius,setC1,setDeltaC,setHCapacity,setTConductivity);

    test.setTimeMax(2.0); //[s]
    test.removeOldFiles();
//    test.restartFile.setFileType(FileType::MULTIPLE_FILES);

    test.solve();

    return 0;

}