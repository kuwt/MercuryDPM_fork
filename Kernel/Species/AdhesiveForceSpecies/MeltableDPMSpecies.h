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

#ifndef MERCURY_MELTABLEDPMSPECIES_H
#define MERCURY_MELTABLEDPMSPECIES_H

#include "Species/BaseSpecies.h"
#include "Math/ExtendedMath.h"
#include "Interactions/AdhesiveForceInteractions/MeltableDPMInteraction.h"

class MeltableDPMSpecies : public virtual BaseSpecies
{
public:
    typedef MeltableDPMInteraction InteractionType;

    MeltableDPMSpecies();

    MeltableDPMSpecies(const MeltableDPMSpecies &p);

    virtual ~MeltableDPMSpecies();

    void read(std::istream& is) override ;

    void write(std::ostream& os) const override;

    std::string getBaseName() const;

    void mix(MeltableDPMSpecies* const S, MeltableDPMSpecies* const T);
    // Species specific funstions:

    //void setEffectiveElasticModulus(Mdouble poisson_ratio_1, Mdouble poisson_ratio_2, Mdouble young_modulus_1, Mdouble young_modulus_2);
    void setEffectiveElasticModulus(Mdouble set_effective_elastic_modulus);
    Mdouble getEffectiveElasticModulus () const;

    //setter - getter:
    void setViscosityCoefficient (Mdouble set_viscosity_coefficient);
    Mdouble getViscosityCoefficient() const;

    void setMaterialStrength (Mdouble set_material_stregnth);
    Mdouble getMaterialStrength() const;
    void setLatentHeat (Mdouble set_latent_heat);
    Mdouble getLatentHeat() const;
    void setThermalConductivity (Mdouble set_thermal_conductivity);
    Mdouble getThermalConductivity() const;
    void setMaterialEmissivity (Mdouble set_material_emissivity);
    Mdouble getMaterialEmissivity() const;
    void setHeatCapacity(Mdouble set_heat_capacity);
    Mdouble getHeatCapacity() const;

    void setMeltingTemperature(Mdouble set_melting_temperature);
    Mdouble getMeltingTemperatures() const;
    void setAmbientTemperature(Mdouble set_ambient_temperature);
    Mdouble getAmbientTemperatures() const;


private:
    Mdouble effective_elastic_modulus_;
    Mdouble viscosity_coefficient_;
    Mdouble material_strength_;
    Mdouble latent_heat_;
    Mdouble thermal_conductivity_;
    Mdouble material_emissivity_;
    Mdouble heat_capacity_;
    Mdouble melting_temperature_;
    Mdouble ambient_temperature_;

};

#endif //MERCURY_MELTABLEDPMSPECIES_H
