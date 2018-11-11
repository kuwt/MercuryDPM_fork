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

#include "MeltableDPMSpecies.h"
//#include "Interactions/BaseInteraction.h"
//#include <cmath>
//#include <Interactions/AdhesiveForceInteractions/MeltableDPMInteraction.h>
//#include "Logger.h"

//#include "Particles/BaseParticle.h"
class BaseParticle;
class BaseInteractable;

MeltableDPMSpecies::MeltableDPMSpecies()
{
    effective_elastic_modulus_ = 0.0;
    viscosity_coefficient_ = 0.0;
    material_strength_ = 0.0;
    latent_heat_ = 0.0;
    thermal_conductivity_ = 0.0;
    material_emissivity_ = 0.0;
    heat_capacity_ = 0.0;
    melting_temperature_ = 0.0;
    ambient_temperature_ = 0.0;

#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"MeltableDPMSpecies::MeltableDPMSpecies() finished"<<std::endl;
#endif
}

MeltableDPMSpecies::MeltableDPMSpecies(const MeltableDPMSpecies &p)
{
    effective_elastic_modulus_ = p.effective_elastic_modulus_;
    viscosity_coefficient_ = p.viscosity_coefficient_;
    material_strength_= p.material_strength_;
    latent_heat_ = p.latent_heat_;
    thermal_conductivity_ = p.thermal_conductivity_;
    material_emissivity_ = p.material_emissivity_;
    heat_capacity_ = p.heat_capacity_;
    melting_temperature_ = p.melting_temperature_;
    ambient_temperature_ = p.ambient_temperature_;

#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"MeltableDPMSpecies::MeltableDPMSpecies(const MeltableDPMSpecies &p) finished"<<std::endl;
#endif
}

MeltableDPMSpecies::~MeltableDPMSpecies()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"MeltableDPMSPecies::~MeltableDPMSpecies() finished"<<std::endl;
#endif
}

void MeltableDPMSpecies::write(std::ostream &os) const
{
    os << "effective alstic modulus" << effective_elastic_modulus_
       << "dissipation" << viscosity_coefficient_
       << "material strength" << material_strength_
       << "latent heat" << latent_heat_
       << "thermal conductivity" << thermal_conductivity_
       << "material_emissivity" << material_emissivity_
       << "heat capacity" << heat_capacity_
       << "melting temperature" << melting_temperature_
       << "ambient temperature" << ambient_temperature_;
}

void MeltableDPMSpecies::read(std::istream &is)
{
    std::string dummy;
    is >> dummy >> effective_elastic_modulus_
       >> dummy >> viscosity_coefficient_
       >> dummy >> material_strength_
       >> dummy >> latent_heat_
       >> dummy >> thermal_conductivity_
       >> dummy >> material_emissivity_
       >> dummy >> heat_capacity_
       >> dummy >> melting_temperature_
       >> dummy >> ambient_temperature_;
}

std::string MeltableDPMSpecies::getBaseName() const
{
    return "MeltableDPM";
}

//specific funstions
void MeltableDPMSpecies::setEffectiveElasticModulus(Mdouble set_effective_elastic_modulus)
{
    effective_elastic_modulus_ = set_effective_elastic_modulus;
/*    if (poisson_ratio_1 >= 0 & poisson_ratio_2 >= 0 & young_modulus_1 >=0 & young_modulus_2 >=0)
    {
        effective_elastic_modulus_ = (4.0/3.0)*((1-mathsFunc::square(poisson_ratio_1)/young_modulus_1)+(1-mathsFunc::square(poisson_ratio_2)/young_modulus_2));
    }
    else
    {
        std::cerr << "Error in setEffectiveElasticModulus(" << effective_elastic_modulus_ << ")" << std::endl;
        exit(-1);
    }*/
}
Mdouble MeltableDPMSpecies::getEffectiveElasticModulus() const
{
    return effective_elastic_modulus_;
}
//
void MeltableDPMSpecies::setViscosityCoefficient(Mdouble set_viscosity_coefficient)
{
    viscosity_coefficient_ = set_viscosity_coefficient;
}
Mdouble MeltableDPMSpecies::getViscosityCoefficient() const
{
    return viscosity_coefficient_;
}

void MeltableDPMSpecies::setMaterialStrength(Mdouble set_material_stregnth)
{
    material_strength_ = set_material_stregnth;
}
Mdouble MeltableDPMSpecies::getMaterialStrength() const
{
    return material_strength_;
}
void MeltableDPMSpecies::setLatentHeat(Mdouble set_latent_heat)
{
    latent_heat_ = set_latent_heat;
}
Mdouble MeltableDPMSpecies::getLatentHeat() const
{
    return latent_heat_;
}
void MeltableDPMSpecies::setThermalConductivity(Mdouble set_thermal_conductivity)
{
    thermal_conductivity_ = set_thermal_conductivity;
}
Mdouble MeltableDPMSpecies::getThermalConductivity() const
{
    return thermal_conductivity_;
}
void MeltableDPMSpecies::setMaterialEmissivity(Mdouble set_material_emissivity)
{
    material_emissivity_ = set_material_emissivity;
}
Mdouble MeltableDPMSpecies::getMaterialEmissivity() const
{
    return material_emissivity_;
}
void MeltableDPMSpecies::setHeatCapacity(Mdouble set_heat_capacity)
{
    heat_capacity_ = set_heat_capacity;
}
Mdouble MeltableDPMSpecies::getHeatCapacity() const
{
    return heat_capacity_;
}
void MeltableDPMSpecies::setAmbientTemperature(Mdouble set_ambient_temperature)
{
    ambient_temperature_ = set_ambient_temperature;
}
Mdouble MeltableDPMSpecies::getAmbientTemperatures() const
{
    return ambient_temperature_;
}
void MeltableDPMSpecies::setMeltingTemperature(Mdouble set_melting_temperature)
{
    melting_temperature_ = set_melting_temperature;
}
Mdouble MeltableDPMSpecies::getMeltingTemperatures() const
{
    return melting_temperature_;
}
//
void MeltableDPMSpecies::mix(MeltableDPMSpecies *const S, MeltableDPMSpecies *const T)
{
    effective_elastic_modulus_ = average(S->getEffectiveElasticModulus(), T->getEffectiveElasticModulus());
}
