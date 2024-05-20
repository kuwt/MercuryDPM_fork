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

#ifndef MERCURY_MELTABLENORMALSPECIES_H
#define MERCURY_MELTABLENORMALSPECIES_H

#include "Species/BaseSpecies.h"
#include "Math/ExtendedMath.h"
#include "Interactions/NormalForceInteractions/MeltableInteraction.h"
#include "DPMBase.h"

class MeltableNormalSpecies : public BaseNormalForce {
public:
    typedef MeltableInteraction InteractionType;

    MeltableNormalSpecies() = default;

    MeltableNormalSpecies(const MeltableNormalSpecies &p) = default;

    ~MeltableNormalSpecies() = default;

    void read(std::istream &is);

    void write(std::ostream &os) const;

    std::string getBaseName() const;

    void mix(MeltableNormalSpecies *S, MeltableNormalSpecies *T);

    void setElasticModulus(Mdouble elasticModulus);

    Mdouble getElasticModulus() const;

    void setPoissonRatio(Mdouble poissonRatio);

    Mdouble getPoissonRatio() const;

    void setDissipation(Mdouble dissipation);

    Mdouble getDissipation() const;

    void setLatentHeat(Mdouble latentHeat);

    Mdouble getLatentHeat() const;

    void setThermalConductivityCoefficient(Mdouble thermalConductivityCoefficient);

    Mdouble getThermalConductivityCoefficient() const;

    void setThermalConvectionCoefficient(Mdouble thermalConvectionCoefficient);

    Mdouble getThermalConvectionCoefficient() const;

    void setMaterialEmissivity(Mdouble materialEmissivity);

    Mdouble getMaterialEmissivity() const;

    void setSolidHeatCapacity(Mdouble heatCapacity);

    Mdouble getSolidHeatCapacity() const;

    void setMeltingTemperature(Mdouble meltingTemperature);

    Mdouble getMeltingTemperature() const;

    void setAmbientTemperature(Mdouble ambientTemperature);

    Mdouble getAmbientTemperature() const;

    void setHeatInput(std::function<double(const BaseParticle*)>& heatInputFunction) {
        heatInputFunction_ = heatInputFunction;
    }

    void setHeatInput(Mdouble heatInput);

    Mdouble getHeatInput(const BaseParticle* p) const;

    //Mdouble computeMaxTimeStep(Mdouble minRadius, Mdouble minDensity, Mdouble minElasticModulus);

    void setMaterialAbsorptivity(Mdouble materialAbsorptivity);

    Mdouble getMaterialAbsorptivity() const;

    void setDeltaT(Mdouble deltaT);

    Mdouble getDeltaT() const;

    void setLiquidHeatCapacity(Mdouble liquidHeatCapacity);

    Mdouble getLiquidHeatCapacity() const;

    void setThermalExpansionCoefficient(Mdouble thermalExpansionCoefficient);

    Mdouble getThermalExpansionCoefficient() const;

    void setActivationEnergy(Mdouble activationEnergy);

    Mdouble getActivationEnergy() const;

    void setSurfaceTension(Mdouble surfaceTension);

    Mdouble getSurfaceTension() const;

    void setRefViscosity(Mdouble refViscosity);

    Mdouble getRefViscosity() const;

    void setMinRelativeSolidRadius(Mdouble minRelativeSolidRadius) {
        minRelativeSolidRadius_ = minRelativeSolidRadius;
    }

    Mdouble getMinRelativeSolidRadius() const {
        return minRelativeSolidRadius_;
    }

    void setWallTemperature(Mdouble wallTemperature) {
        wallTemperature_ = wallTemperature;
    }

    Mdouble getWallTemperature() const {
        return wallTemperature_;
    }

    Mdouble getEffectiveHeatCapacity(double temperature) const;

    Mdouble getEffectiveLatentHeat() const;

    Mdouble getRelativeSolidRadius(double temperature) const;

    Mdouble getEffectiveElasticModulus() const {
        return 0.5*getElasticModulus()/(1.0-mathsFunc::square(getPoissonRatio()));
    }

    void analyseTimeScales(double radius, double density, double temperature) const;

private:
    Mdouble elasticModulus_ = 0.0;
    Mdouble poissonRatio_ = 0.0;
    Mdouble dissipation_ = 0.0;
    Mdouble deltaT_ = 0.0;
    Mdouble solidHeatCapacity_ = 0.0;
    Mdouble liquidHeatCapacity_ = 0.0;
    Mdouble latentHeat_ = 0.0;
    Mdouble meltingTemperature_ = constants::inf;
    Mdouble thermalConductivityCoefficient_ = 0.0;
    Mdouble thermalConvectionCoefficient_ = 0.0;
    Mdouble materialEmissivity_ = 0.0; //rename
    Mdouble ambientTemperature_ = 0.0;
    Mdouble wallTemperature_ = -1.0;
    Mdouble heatInput_ = 0;
    std::function<double(const BaseParticle*)> heatInputFunction_ = nullptr;
    Mdouble materialAbsorptivity_ = 0.0;
    Mdouble thermalExpansionCoefficient_ = 0.0;
    Mdouble activationEnergy_ = 0.0;
    Mdouble surfaceTension_ = 0.0;
    Mdouble refViscosity_ = 0.0;
    Mdouble minRelativeSolidRadius_ = 0.1;
};

#endif
