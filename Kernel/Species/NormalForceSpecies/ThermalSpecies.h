//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
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

#ifndef THERMALSPECIES_H
#define THERMALSPECIES_H

#include "Interactions/NormalForceInteractions/ThermalInteraction.h"

class BaseInteraction;

template<class NormalForceSpecies>
class ThermalSpecies : public NormalForceSpecies
{
public:
    
    ///\brief The correct Interaction type for this FrictionForceSpecies
    typedef ThermalInteraction<typename NormalForceSpecies::InteractionType> InteractionType;
    
    ///\brief The default constructor.
    ThermalSpecies();
    
    ///\brief The default copy constructor.
    ThermalSpecies(const ThermalSpecies& s);
    
    ///\brief The default destructor.
    virtual ~ThermalSpecies();
    
    /// \brief Writes the species properties to an output stream.
    void write(std::ostream& os) const;
    
    /// \brief Reads the species properties from an input stream.
    void read(std::istream& is);
    
    /// \brief Used in Species::getName to obtain a unique name for each Species.
    std::string getBaseName() const;
    
    ///Allows heatCapacity_ to be accessed
    Mdouble getHeatCapacity() const;
    
    ///Allows heatCapacity_ to be changed
    void setHeatCapacity(Mdouble heatCapacity);
    
    ///Allows heatCapacity_ to be accessed
    Mdouble getThermalConductivity() const;
    
    ///Allows heatCapacity_ to be changed
    void setThermalConductivity(Mdouble thermalConductivity);

private:
    /*!
     * \brief The heat capacity.
     */
    Mdouble heatCapacity_;
    
    /*!
     * \brief The thermal conductivity.
     */
    Mdouble thermalConductivity_;
};

template<class NormalForceSpecies>
ThermalSpecies<NormalForceSpecies>::ThermalSpecies()
        : NormalForceSpecies()
{
    heatCapacity_ = 0.0;
    thermalConductivity_ = 0.0;
}

template<class NormalForceSpecies>
ThermalSpecies<NormalForceSpecies>::ThermalSpecies(const ThermalSpecies& s)
        : NormalForceSpecies(s)
{
    heatCapacity_ = s.heatCapacity_;
    thermalConductivity_ = s.thermalConductivity_;
}

template<class NormalForceSpecies>
ThermalSpecies<NormalForceSpecies>::~ThermalSpecies()
{}

template<class NormalForceSpecies>
void ThermalSpecies<NormalForceSpecies>::write(std::ostream& os) const
{
    NormalForceSpecies::write(os);
    os << " heatCapacity " << heatCapacity_;
    os << " thermalConductivity " << thermalConductivity_;
}

template<class NormalForceSpecies>
void ThermalSpecies<NormalForceSpecies>::read(std::istream& is)
{
    std::string dummy;
    NormalForceSpecies::read(is);
    is >> dummy >> heatCapacity_;
    is >> dummy >> thermalConductivity_;
}

template<class NormalForceSpecies>
std::string ThermalSpecies<NormalForceSpecies>::getBaseName() const
{
    return "Thermal" + NormalForceSpecies::getBaseName();
}


template<class NormalForceSpecies>
Mdouble ThermalSpecies<NormalForceSpecies>::getHeatCapacity() const
{
    return heatCapacity_;
}

template<class NormalForceSpecies>
void ThermalSpecies<NormalForceSpecies>::setHeatCapacity(Mdouble heatCapacity)
{
    logger.assert_always(heatCapacity > 0,
                         "[ThermalSpecies<>::setHeatCapacity(%)] value has to be positive",
                         heatCapacity);
    heatCapacity_ = heatCapacity;
}

template<class NormalForceSpecies>
Mdouble ThermalSpecies<NormalForceSpecies>::getThermalConductivity() const
{
    return thermalConductivity_;
}

template<class NormalForceSpecies>
void ThermalSpecies<NormalForceSpecies>::setThermalConductivity(Mdouble thermalConductivity)
{
    logger.assert_always(thermalConductivity >= 0,
                         "[ThermalSpecies<>::setThermalConductivity(%)] value has to be positive",
                         thermalConductivity);
    thermalConductivity_ = thermalConductivity;
}

#endif
