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

#include "MeltableDPMInteraction.h"
#include "Species/AdhesiveForceSpecies/MeltableDPMSpecies.h"
#include "Particles/BaseParticle.h"
#include "Particles/ThermalParticle.h"
//#include "BaseInteractable.h"
#include "InteractionHandler.h"
#include <iomanip>

MeltableDPMInteraction::MeltableDPMInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"MeltableDPMInteraction::MeltableDPMInteraction() finished"<<std::endl;
#endif
}

MeltableDPMInteraction::MeltableDPMInteraction(const MeltableDPMInteraction& p)
        : BaseInteraction(p)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"MeltableDPMInteraction::MeltableDPMInteraction(const MeltableDPMInteraction& p) finished"<<std::endl;
#endif
}

MeltableDPMInteraction::~MeltableDPMInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"MeltableDPMInteraction::~MeltableDPMInteraction() finished"<<std::endl;
#endif
}
//
void MeltableDPMInteraction::write(std::ostream& os UNUSED) const
{
    //BaseInteraction::write(os);
}

void MeltableDPMInteraction::read(std::istream& is UNUSED)
{
    //BaseInteraction::read(is);
}
//
std::string MeltableDPMInteraction::getBaseName() const
{
    return "MeltableDPM";
}

void MeltableDPMInteraction::computeNormalForce()
{
    ThermalParticle* i = dynamic_cast<ThermalParticle*>(getP());
    logger.assert(i!=null_ptr,"i is not a Thermal Particle");
    ThermalParticle* j = dynamic_cast<ThermalParticle*>(getI());
    logger.assert(j!=null_ptr,"j is not a Thermal Particle");
    //
    Mdouble molten_layer_thickness_i = 0.0;
    Mdouble molten_layer_thickness_j = 0.0;
    Mdouble solid_radius_i = i->getRadius() - molten_layer_thickness_i;
    Mdouble solid_radius_j = j->getRadius() - molten_layer_thickness_j;
    //
    //setRelativeVelocity(getP()->getVelocityAtContact(getContactPoint()))-getI()->getVelocityAtContact(getContactPoint());
    //setNormalRelativeVelocity(Vec3D::dot(getRelativeVelocity(),getNormal()));
    //
    const MeltableDPMSpecies* species = getSpecies();
    //
    Mdouble stefan_boltzman_constant = 5.67e-8;
    Mdouble particle_area = constants::pi*getEffectiveRadius();
    //
    Mdouble heat_source;
    /*if ()
        heat_source = ;
    else
        heat_source = 0;*/
    //
    if (getOverlap() >= 0) {
        // radiative term:
        Mdouble heat_transfer_radiative_i =
                -1 * species->getMaterialEmissivity() * stefan_boltzman_constant * particle_area *
                ((mathsFunc::cubic(i->getTemperature()) * i->getTemperature()) -
                 (mathsFunc::cubic(species->getAmbientTemperatures()) * species->getAmbientTemperatures()));
        // conductive term:
        Mdouble contact_radius = sqrt(getEffectiveRadius() * getOverlap());
        Mdouble heat_conduction_ij =
                2 * contact_radius * species->getThermalConductivity() * (j->getTemperature() - i->getTemperature());
        // Total heat transfer between two particles:
        Mdouble heat_transfer = heat_conduction_ij + heat_transfer_radiative_i;
        // dT/dt:
        Mdouble temperature_evolution = (heat_transfer + heat_source) / (i->getMass() * species->getHeatCapacity());
        // Ti:
        Mdouble temperature_i;
        temperature_i += temperature_evolution * getTimeStamp();
        i->setTemperature(temperature_i);
        //Tj:
        Mdouble temperature_j;
        temperature_j += temperature_evolution * getTimeStamp();
        j->setTemperature(temperature_j);
        //
        if (i->getTemperature() >= species->getMeltingTemperatures() & (0 <= molten_layer_thickness_i & molten_layer_thickness_i <= i->getRadius() ))
        {
            Mdouble molten_layer_thickness_evolution_i = 0;
        }
        //
    }
    else
    {

    }
}

Mdouble MeltableDPMInteraction::getElasticEnergy() const
{
    return 0.0;
}

/*!
 * \return a constant pointer to an instance of this class.
 */
const MeltableDPMSpecies* MeltableDPMInteraction::getSpecies() const
{
    return dynamic_cast<const SpeciesType*>(getBaseSpecies()); //downcast
}