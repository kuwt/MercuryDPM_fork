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
#include <Walls/InfiniteWall.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <iomanip>

/** This code sinters particles and compress them at the same time.
*/
class Densification_PolyamideS1 : public Mercury3D{
public:

    //set default values
    Densification_PolyamideS1()
    {
        setName("DensificationPolyamideS1");
        readRestartFile();
        setRestarted(false);
//        fStatFile.setFileType(FileType::MULTIPLE_FILES);
        setName("DensificationPolyamideS2");
        particleSpecies = dynamic_cast<ThermalSinterLinFrictionReversibleAdhesiveSpecies*>(speciesHandler.getObject(0));
        particleSpeciesAtMeltingPoint = new ThermalSinterLinFrictionReversibleAdhesiveSpecies(*particleSpecies);
//        lid = dynamic_cast<InfiniteWall*>(wallHandler.getLastObject());
    }

    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpecies;
    ThermalSinterLinFrictionReversibleAdhesiveSpecies* particleSpeciesAtMeltingPoint;
    Mdouble meltTemperature;
    Mdouble varTemperature;
    Mdouble temperatureGradient;
    Mdouble deltaRadius;
    Mdouble maxTemperature;

    ~Densification_PolyamideS1()
    {
        delete particleSpeciesAtMeltingPoint;
    }

    //set slave variables
    void setupInitialConditions() override
    {
        Mdouble stiffnessFactor = 0.5*(1.0+tanh((meltTemperature- getTemperature())/varTemperature));
        particleSpeciesAtMeltingPoint->setLoadingStiffness(particleSpecies->getLoadingStiffness()/stiffnessFactor);
        particleSpeciesAtMeltingPoint->setCohesionStiffness(particleSpecies->getCohesionStiffness()/stiffnessFactor);
    }

    Mdouble getTemperature() const
    {
        Mdouble heatingTemperature = 20.0 + temperatureGradient*getTime();
        Mdouble coolingTemperature = 20.0 + temperatureGradient*(getTimeMax() - getTime());
        return std::min(std::min(heatingTemperature,coolingTemperature),maxTemperature);
    }

    void setTemperature(Mdouble temperature)
    {
//        static int counter = 0;
//        counter++;
        static Mdouble oldTemperature = 20.0;
        Mdouble deltaTemperature = temperature-oldTemperature;
        Mdouble factorRadius = 1.0-deltaRadius*deltaTemperature;
//
//        //change the density and particle radius
        particleSpecies->setDensity(mathsFunc::cubic(factorRadius)*particleSpecies->getDensity());
        for (BaseParticle* p : particleHandler)
        {
            if (p->getSpecies()==particleSpecies)
                p->setRadius(factorRadius*p->getRadius());
        }

        //change species properties
        Mdouble stiffnessFactor = 0.5*(1+tanh((meltTemperature-temperature)/varTemperature));
        Mdouble oldLoadingStiffness = particleSpecies->getLoadingStiffness();
        particleSpecies->setLoadingStiffness(stiffnessFactor*particleSpeciesAtMeltingPoint->getLoadingStiffness());
        //particleSpecies->setUnloadingStiffnessMax(stiffnessFactor*particleSpeciesAtMeltingPoint->getUnloadingStiffnessMax());
        particleSpecies->setCohesionStiffness(stiffnessFactor*particleSpeciesAtMeltingPoint->getCohesionStiffness());
//
        //for decreasing temperature, change the maxOverlap
        if (deltaTemperature<0.0)
        {
//            for (BaseInteraction* cBase : interactionHandler)
//            {
//                auto c = dynamic_cast<ThermalSinterLinFrictionReversibleAdhesiveSpecies*>(cBase);
//                Mdouble unloadingStiffness = c->getUnloadingStiffness();
//                c->setMaxOverlap(c->getMaxOverlap()
//                                 *(unloadingStiffness-oldLoadingStiffness)
//                                 /(unloadingStiffness-particleSpecies->getLoadingStiffness())
//                );
//            }
            exit(-1);
        }

//        change old temperature

    }


    void actionsAfterTimeStep() override
    {
        setTemperature(getTemperature());
//
//        static unsigned int counter = 0;
//        if (++counter>100)
//        {
//            counter=0;
//
//            static Mdouble lidArea = 0.25 * constants::pi * mathsFunc::square(getXMax()- getXMin());
//            static Mdouble particleArea = constants::pi * mathsFunc::square(particleHandler.getObject(0)->getRadius());
//            static Mdouble stiffness = dynamic_cast<const LinearPlasticViscoelasticSpecies*>(speciesHandler.getObject(0))->getLoadingStiffness();
//            // amount by which the pressure has to be increased
//            Mdouble dPressure = lid->getForce().Z/lidArea - pressure;
//            // amount by which position should be changed to achieve the right pressure
//            Mdouble dZ = dPressure * particleArea / stiffness;
//            lid->setVelocity(Vec3D(0.0,0.0,dZ * 100.0));
//        }
    }

    void printTime() const override
    {
//        static Mdouble lidArea = 0.25 * constants::pi * mathsFunc::square(getXMax()- getXMin());
        std::cout << "t=" << std::left << std::setw(4) << getTime()
                  << " T " << std::left << std::setw(5) << getTemperature()
                  << " k " << std::left << std::setw(7) << particleSpecies->getLoadingStiffness()
                  << " meanOverlap " << std::left << std::setw(8) << 1000*std::sqrt(getElasticEnergy()/particleSpecies->getLoadingStiffness()/interactionHandler.getNumberOfObjects())
//                  //<< " d " << std::left << std::setw(8) << 1000*dynamic_cast<const LinearPlasticViscoelasticInteraction*>(interactionHandler.getObject(0))->getOverlap()
//                  //<< " dMax " << std::left << std::setw(8) << 1000*dynamic_cast<const LinearPlasticViscoelasticInteraction*>(interactionHandler.getObject(0))->getPlasticOverlap()
//                  //<< " Ene " << std::left << std::setw(6) << getKineticEnergy()/getElasticEnergy()
//                  << " Lid " << std::left << std::setw(6) << lid->getPosition().Z
                  << std::endl;
    }



//    InfiniteWall* lid;
//    Mdouble pressure;

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    Densification_PolyamideS1 s;
//    s.pressure = 200;
    s.meltTemperature = 100.0;
    s.maxTemperature = 180.0;
    s.varTemperature = 20.0;
    s.deltaRadius = 0.1e-3;
    s.setTimeMax(0.003);
    s.temperatureGradient = 25.0;
    s.solve();
}
///todo{ I need to change the amount of overlap according to the experimental calibration}