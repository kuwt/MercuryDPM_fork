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

#ifndef HEATFLUIDCOUPLEDSPECIES_H
#define HEATFLUIDCOUPLEDSPECIES_H

class BaseInteraction;

template<class NormalForceSpecies>
class HeatFluidCoupledSpecies : public NormalForceSpecies
{
public:
    
   // ///\brief The correct Interaction type for this FrictionForceSpecies
    //typename NormalForceSpecies::InteractionType InteractionType;
    
    ///\brief The default constructor.
    HeatFluidCoupledSpecies();
    
    ///\brief The default copy constructor.
    HeatFluidCoupledSpecies(const HeatFluidCoupledSpecies& s);
    
    ///\brief The default destructor.
    virtual ~HeatFluidCoupledSpecies();
    
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
    
    ///Allows thermalConductivity_ to be accessed
    Mdouble getThermalConductivity() const;
    
    ///Allows thermalConductivity_ to be changed
    void setThermalConductivity(Mdouble thermalConductivity);

    ///Allows massTransferCoefficient_ to be accessed
    Mdouble getMassTransferCoefficient() const;

    ///Allows massTransferCoefficient_ to be changed
    void setMassTransferCoefficient(Mdouble massTransferCoefficient);

    ///Allows latentHeatVaporization_ to be accessed
    Mdouble getLatentHeatVaporization() const;

    ///Allows latentHeatVaporization_ to be changed
    void setLatentHeatVaporization(Mdouble latentHeatVaporization);

    ///Allows liquidDensity_ to be accessed
    Mdouble getLiquidDensity() const;

    ///Allows liquidDensity_ to be changed
    void setLiquidDensity(Mdouble liquidDensity);

    ///Allows evaporationCoefficientA_ to be accessed
    Mdouble getEvaporationCoefficientA() const;

    ///Allows evaporationCoefficientA_ to be changed
    void setEvaporationCoefficientA(Mdouble evaporationCoefficientA);

    ///Allows evaporationCoefficientB_ to be accessed
    Mdouble getEvaporationCoefficientB() const;

    ///Allows evaporationCoefficientB_ to be changed
    void setEvaporationCoefficientB(Mdouble evaporationCoefficientB);

    ///Allows ambientHumidity_ to be accessed
    Mdouble getAmbientHumidity() const;

    ///Allows ambientHumidity_ to be changed
    void setAmbientHumidity(Mdouble ambientHumidity);

    ///Allows ambientEquilibriumMoistureContent_ to be accessed
    Mdouble getAmbientEquilibriumMoistureContent() const;

    ///Allows ambientEquilibriumMoistureContent_ to be changed
    void setAmbientEquilibriumMoistureContent(Mdouble ambientEquilibriumMoistureContent);

    ///Allows ambientVapourConcentration_ to be accessed
    Mdouble getAmbientVapourConcentration() const;

    ///Allows ambientVapourConcentration_ to be changed
    void setAmbientVapourConcentration(Mdouble ambientVapourConcentration);

    ///Allows ambientTemperature_ to be accessed
    Mdouble getAmbientTemperature() const;

    ///Allows ambientTemperature_ to be changed
    void setAmbientTemperature(Mdouble ambientTemperature);

private:
    /*!
     * \brief The heat capacity.
     */
    Mdouble heatCapacity_;
    
    /*!
     * \brief The thermal conductivity.
     */
    Mdouble thermalConductivity_;

    /*!
     * \brief The mass transfer coefficient.
     */
    Mdouble massTransferCoefficient_;

    /*!
     * \brief The latent heat of vaporization.
     */
    Mdouble latentHeatVaporization_;

    /*!
     * \brief The liquid density.
     */
    Mdouble liquidDensity_;

    /*!
     * \brief The evaporation coefficient a.
     */
    Mdouble evaporationCoefficientA_;

    /*!
     * \brief The evaporation coefficient b.
     */
    Mdouble evaporationCoefficientB_;

    /*!
     * \brief The ambient humidity.
     */
    Mdouble ambientHumidity_;

    /*!
     * \brief The ambient equilibrium moisture content.
     */
    Mdouble ambientEquilibriumMoistureContent_;

    /*!
     * \brief The ambient vapour concentration.
     */
    Mdouble ambientVapourConcentration_;

    /*!
     * \brief The ambient temperature.
     */
    Mdouble ambientTemperature_;

};

template<class NormalForceSpecies>
HeatFluidCoupledSpecies<NormalForceSpecies>::HeatFluidCoupledSpecies()
        : NormalForceSpecies()
{
    heatCapacity_ = 0.0;
    thermalConductivity_ = 0.0;
    massTransferCoefficient_=0.0;
    latentHeatVaporization_=0.0;
    liquidDensity_=0.0;
    evaporationCoefficientA_=0.0;
    evaporationCoefficientB_=0.0;
    ambientHumidity_=0.0;
    ambientEquilibriumMoistureContent_=0.0;
    ambientVapourConcentration_=0.0;
    ambientTemperature_=0.0;
}

template<class NormalForceSpecies>
HeatFluidCoupledSpecies<NormalForceSpecies>::HeatFluidCoupledSpecies(const HeatFluidCoupledSpecies& s)
        : NormalForceSpecies(s)
{
    heatCapacity_ = s.heatCapacity_;
    thermalConductivity_ = s.thermalConductivity_;
    massTransferCoefficient_= s.massTransferCoefficient_;
    latentHeatVaporization_= s.latentHeatVaporization_;
    liquidDensity_=s.liquidDensity_;
    evaporationCoefficientA_=s.evaporationCoefficientA_;
    evaporationCoefficientB_=s.evaporationCoefficientB_;
    ambientHumidity_=s.ambientHumidity_;
    ambientEquilibriumMoistureContent_=s.ambientEquilibriumMoistureContent_;
    ambientVapourConcentration_=s.ambientVapourConcentration_;
    ambientTemperature_=s.ambientTemperature_;
}

template<class NormalForceSpecies>
HeatFluidCoupledSpecies<NormalForceSpecies>::~HeatFluidCoupledSpecies()
{}

template<class NormalForceSpecies>
void HeatFluidCoupledSpecies<NormalForceSpecies>::write(std::ostream& os) const
{
    NormalForceSpecies::write(os);
    os << " heatCapacity " << heatCapacity_;
    os << " thermalConductivity " << thermalConductivity_;
    os << " massTransferCoefficient " << massTransferCoefficient_;
    os << " latentHeatVaporization " << latentHeatVaporization_;
    os << " liquidDensity " << liquidDensity_;
    os << " evaporationCoefficientA " << evaporationCoefficientA_;
    os << " evaporationCoefficientA " << evaporationCoefficientB_;
    os << " ambientHumidity " << ambientHumidity_;
    os << " ambientEquilibriumMoistureContent " << ambientEquilibriumMoistureContent_;
    os << " ambientVapourConcentration " << ambientVapourConcentration_;
    os << " ambientTemperature " << ambientTemperature_;
}

template<class NormalForceSpecies>
void HeatFluidCoupledSpecies<NormalForceSpecies>::read(std::istream& is)
{
    std::string dummy;
    NormalForceSpecies::read(is);
    is >> dummy >> heatCapacity_;
    is >> dummy >> thermalConductivity_;
    is >> dummy >> massTransferCoefficient_;
    is >> dummy >> latentHeatVaporization_;
    is >> dummy >> liquidDensity_;
    is >> dummy >> evaporationCoefficientA_;
    is >> dummy >> evaporationCoefficientB_;
    is >> dummy >> ambientHumidity_;
    is >> dummy >> ambientEquilibriumMoistureContent_;
    is >> dummy >> ambientVapourConcentration_;
    is >> dummy >> ambientTemperature_;
}

template<class NormalForceSpecies>
std::string HeatFluidCoupledSpecies<NormalForceSpecies>::getBaseName() const
{
    return "HeatFluidCoupled" + NormalForceSpecies::getBaseName();
}


template<class NormalForceSpecies>
Mdouble HeatFluidCoupledSpecies<NormalForceSpecies>::getHeatCapacity() const
{
    return heatCapacity_;
}

template<class NormalForceSpecies>
void HeatFluidCoupledSpecies<NormalForceSpecies>::setHeatCapacity(Mdouble heatCapacity)
{
    logger.assert_always(heatCapacity > 0,
                         "[HeatFluidCoupledSpecies<>::setHeatCapacity(%)] value has to be positive",
                         heatCapacity);
    heatCapacity_ = heatCapacity;
}

template<class NormalForceSpecies>
Mdouble HeatFluidCoupledSpecies<NormalForceSpecies>::getThermalConductivity() const
{
    return thermalConductivity_;
}

template<class NormalForceSpecies>
void HeatFluidCoupledSpecies<NormalForceSpecies>::setThermalConductivity(Mdouble thermalConductivity)
{
    logger.assert_always(thermalConductivity >= 0,
                         "[HeatFluidCoupledSpecies<>::setThermalConductivity(%)] value has to be positive",
                         thermalConductivity);
    thermalConductivity_ = thermalConductivity;
}

template<class NormalForceSpecies>
Mdouble HeatFluidCoupledSpecies<NormalForceSpecies>::getMassTransferCoefficient() const
{
    return massTransferCoefficient_;
}

template<class NormalForceSpecies>
void HeatFluidCoupledSpecies<NormalForceSpecies>::setMassTransferCoefficient(Mdouble massTransferCoefficient)
{
    logger.assert_always(massTransferCoefficient > 0,
                         "[HeatFluidCoupledSpecies<>::setMassTransferCoefficient(%)] value has to be positive",
                         massTransferCoefficient);
    massTransferCoefficient_ = massTransferCoefficient;
}

template<class NormalForceSpecies>
Mdouble HeatFluidCoupledSpecies<NormalForceSpecies>::getLatentHeatVaporization() const
{
    return latentHeatVaporization_;
}

template<class NormalForceSpecies>
void HeatFluidCoupledSpecies<NormalForceSpecies>::setLatentHeatVaporization(Mdouble latentHeatVaporization)
{
    logger.assert_always(latentHeatVaporization > 0,
                         "[HeatFluidCoupledSpecies<>::setLatentHeatVaporization(%)] value has to be positive",
                         latentHeatVaporization);
    latentHeatVaporization_ = latentHeatVaporization;
}

template<class NormalForceSpecies>
Mdouble HeatFluidCoupledSpecies<NormalForceSpecies>::getLiquidDensity() const
{
    return liquidDensity_;
}

template<class NormalForceSpecies>
void HeatFluidCoupledSpecies<NormalForceSpecies>::setLiquidDensity(Mdouble liquidDensity)
{
    logger.assert_always(liquidDensity > 0,
                         "[HeatFluidCoupledSpecies<>::setLiquidDensity(%)] value has to be positive",
                         liquidDensity);
    liquidDensity_ = liquidDensity;
}

template<class NormalForceSpecies>
Mdouble HeatFluidCoupledSpecies<NormalForceSpecies>::getEvaporationCoefficientA() const
{
    return evaporationCoefficientA_;
}

template<class NormalForceSpecies>
void HeatFluidCoupledSpecies<NormalForceSpecies>::setEvaporationCoefficientA(Mdouble evaporationCoefficientA)
{
    logger.assert_always(evaporationCoefficientA > 0,
                         "[HeatFluidCoupledSpecies<>::setEvaporationCoefficientA(%)] value has to be positive",
                         evaporationCoefficientA);
    evaporationCoefficientA_ = evaporationCoefficientA;
}

template<class NormalForceSpecies>
Mdouble HeatFluidCoupledSpecies<NormalForceSpecies>::getEvaporationCoefficientB() const
{
    return evaporationCoefficientB_;
}

template<class NormalForceSpecies>
void HeatFluidCoupledSpecies<NormalForceSpecies>::setEvaporationCoefficientB(Mdouble evaporationCoefficientB)
{
    logger.assert_always(evaporationCoefficientB < 0,
                         "[HeatFluidCoupledSpecies<>::setEvaporationCoefficientB(%)] value has to be negative",
                         evaporationCoefficientB);
    evaporationCoefficientB_ = evaporationCoefficientB;
}

template<class NormalForceSpecies>
Mdouble HeatFluidCoupledSpecies<NormalForceSpecies>::getAmbientHumidity() const
{
    return ambientHumidity_;
}

template<class NormalForceSpecies>
void HeatFluidCoupledSpecies<NormalForceSpecies>::setAmbientHumidity(Mdouble ambientHumidity)
{
    logger.assert_always(ambientHumidity > 0,
                         "[HeatFluidCoupledSpecies<>::setAmbientHumidity(%)] value has to be positive",
                         ambientHumidity);
    ambientHumidity_ = ambientHumidity;
}

template<class NormalForceSpecies>
Mdouble HeatFluidCoupledSpecies<NormalForceSpecies>::getAmbientEquilibriumMoistureContent() const
{
    return ambientEquilibriumMoistureContent_;
}

template<class NormalForceSpecies>
void HeatFluidCoupledSpecies<NormalForceSpecies>::setAmbientEquilibriumMoistureContent(Mdouble ambientEquilibriumMoistureContent)
{
    logger.assert_always(ambientEquilibriumMoistureContent > 0,
                         "[HeatFluidCoupledSpecies<>::setAmbientEquilibriumMoistureContent(%)] value has to be positive",
                         ambientEquilibriumMoistureContent);
    ambientEquilibriumMoistureContent_ = ambientEquilibriumMoistureContent;
}

template<class NormalForceSpecies>
Mdouble HeatFluidCoupledSpecies<NormalForceSpecies>::getAmbientVapourConcentration() const
{
    return ambientVapourConcentration_;
}

template<class NormalForceSpecies>
void HeatFluidCoupledSpecies<NormalForceSpecies>::setAmbientVapourConcentration(Mdouble ambientVapourConcentration)
{
    logger.assert_always(ambientVapourConcentration > 0,
                         "[HeatFluidCoupledSpecies<>::setAmbientVapourConcentration(%)] value has to be positive",
                         ambientVapourConcentration);
    ambientVapourConcentration_ = ambientVapourConcentration;
}

template<class NormalForceSpecies>
Mdouble HeatFluidCoupledSpecies<NormalForceSpecies>::getAmbientTemperature() const
{
    return ambientTemperature_;
}

template<class NormalForceSpecies>
void HeatFluidCoupledSpecies<NormalForceSpecies>::setAmbientTemperature(Mdouble ambientTemperature)
{
    logger.assert_always(ambientTemperature > 0,
                         "[HeatFluidCoupledSpecies<>::setAmbientTemperature(%)] value has to be positive",
                         ambientTemperature);
    ambientTemperature_ = ambientTemperature;
}

#endif
