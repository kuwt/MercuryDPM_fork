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

#ifndef HEATFLUIDCOUPLEDSPECIES_H
#define HEATFLUIDCOUPLEDSPECIES_H
#include "ThermalSpecies.h"
#include "Interactions/NormalForceInteractions/HeatFluidCoupledInteraction.h"
class BaseInteraction;

/// Species for the HeatFluidCoupledParticle
template<class NormalForceSpecies>
class HeatFluidCoupledSpecies : public ThermalSpecies<NormalForceSpecies>
{
public:

    typedef HeatFluidCoupledInteraction<typename NormalForceSpecies::InteractionType> InteractionType;

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

    void actionsAfterTimeStep(BaseParticle* particle) const override;

    /// f1 is used in Runge–Kutta method.
    std::array<double,2> f(double liquidVolume,double temperature, double mass, double surfaceArea) const;

private:

    /*!
     * \brief The mass transfer rate (m/s)
     */
    Mdouble massTransferCoefficient_;

    /*!
     * \brief The latent heat of vaporization (J/kg)
     */
    Mdouble latentHeatVaporization_;

    /*!
     * \brief The liquid density (kg/m^3)
     */
    Mdouble liquidDensity_;

    /*!
     * \brief The evaporation coefficient a (dimensionless)
     */
    Mdouble evaporationCoefficientA_;

    /*!
     * \brief The evaporation coefficient b (dimensionless)
     */
    Mdouble evaporationCoefficientB_;

    /*!
     * \brief The ambient humidity (dimensionless, between 0 and 1, but cannot be 0)
     */
    Mdouble ambientHumidity_;

    /*!
     * \brief The ambient equilibrium moisture content (dimensionless, between 0 and 1).
     */
    Mdouble ambientEquilibriumMoistureContent_;

    /*!
     * \brief The ambient vapour concentration (kg/m^3).
     */
    Mdouble ambientVapourConcentration_;

    /*!
     * \brief The ambient temperature (K).
     */
    Mdouble ambientTemperature_;

};

template<class NormalForceSpecies>
HeatFluidCoupledSpecies<NormalForceSpecies>::HeatFluidCoupledSpecies()
        : ThermalSpecies<NormalForceSpecies>()
{
    massTransferCoefficient_=0.0;
    latentHeatVaporization_=0.0;
    liquidDensity_=1.0;
    evaporationCoefficientA_=0.0;
    evaporationCoefficientB_=0.0;
    // note, this default value is unrealistically high; 0.3 is realistic.
    ambientHumidity_=1.0;
    ambientEquilibriumMoistureContent_=0.0;
    ambientVapourConcentration_=0.0;
    ambientTemperature_=0.0;
}

template<class NormalForceSpecies>
HeatFluidCoupledSpecies<NormalForceSpecies>::HeatFluidCoupledSpecies(const HeatFluidCoupledSpecies& s)
        : ThermalSpecies<NormalForceSpecies>(s)
{
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
    ThermalSpecies<NormalForceSpecies>::write(os);
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
    ThermalSpecies<NormalForceSpecies>::read(is);
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
    logger.assert_always(evaporationCoefficientA >= 0,
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
    logger.assert_always(evaporationCoefficientB <= 0,
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
    //note, this is not allowed to be zero!
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
    logger.assert_always(ambientEquilibriumMoistureContent >= 0,
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
    logger.assert_always(ambientVapourConcentration >= 0,
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

template<class NormalForceSpecies>
void HeatFluidCoupledSpecies<NormalForceSpecies>::actionsAfterTimeStep(BaseParticle* baseParticle) const
{
    double dt = baseParticle->getHandler()->getDPMBase()->getTimeStep();
    auto p = dynamic_cast<HeatFluidCoupledParticle*>(baseParticle);
    double mass = p->getMass();
    double surfaceArea = p->getSurfaceArea();
    //Runge–Kutta method
    std::array<double,2> k1 = f(p->getLiquidVolume(),p->getTemperature(), mass, surfaceArea);
    //logger(INFO, "dL % dT %", k1[0], k1[1]);
    std::array<double,2> k2 = f(p->getLiquidVolume()+dt*0.5*k1[0],p->getTemperature()+dt*0.5*k1[1], mass, surfaceArea);
    std::array<double,2> k3 = f(p->getLiquidVolume()+dt*0.5*k2[0],p->getTemperature()+dt*0.5*k2[1], mass, surfaceArea);
    std::array<double,2> k4 = f(p->getLiquidVolume()+dt*k3[0],p->getTemperature()+dt*k3[1], mass, surfaceArea);
    double dliquidVolume = dt*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6.0;
    double dTemperature = dt*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6.0;
    // if liquid film volume is larger than the volume that needs to be subtracted
    if (p->getLiquidVolume()+dliquidVolume>=0.0) {
        p->setLiquidVolume(p->getLiquidVolume()+dliquidVolume);
        p->setTemperature(std::max(0.0,p->getTemperature()+dTemperature));
        //logger(INFO,"% LF % T %", p->getIndex(), p->getLiquidVolume(), p->getTemperature());
    } else {
        // how much to remove from liquid bridges
        double liquidVolumeToDistribute = - (p->getLiquidVolume()+dliquidVolume);
        p->setLiquidVolume(0.0);
        // get liquid bridge volume
        double liquidBridgeVolume = 0.0;
        for (BaseInteraction* i : p-> getInteractions()) {
            auto j = dynamic_cast<LiquidMigrationWilletInteraction*>(i);
            liquidBridgeVolume += j->getLiquidBridgeVolume();
        }
        // if liquid bridge volume is larger than the volume that needs to be subtracted
        if (liquidVolumeToDistribute<=liquidBridgeVolume) {
            double factor = 1.0-liquidVolumeToDistribute/liquidBridgeVolume;
            for (BaseInteraction* i : p-> getInteractions()) {
                auto j = dynamic_cast<LiquidMigrationWilletInteraction*>(i);
                j->setLiquidBridgeVolume(factor*j->getLiquidBridgeVolume());
            }
            p->setTemperature(std::max(0.0,p->getTemperature()+dTemperature));
        } else {
            // if both liquid bridges and liquid films are empty after the subtraction
            if (liquidBridgeVolume!=0.0) {
                liquidVolumeToDistribute -= liquidBridgeVolume;
                for (BaseInteraction *i: p->getInteractions()) {
                    auto j = dynamic_cast<LiquidMigrationWilletInteraction *>(i);
                    j->setLiquidBridgeVolume(0.0);
                }
            }
            double factor = 1.0+liquidVolumeToDistribute/dliquidVolume;
            p->setTemperature(std::max(0.0,p->getTemperature()+factor*dTemperature));
        }
    }
}

/**
 * Computes the drying and cooling rate of a particle; based on equations in (Azmir et al., 2018), which are summarised in EvaporationModel.pdf
 * @param liquidVolume_  Liquid film volume
 * @param temperature_   Temperature of particle
 * @param mass           Mass of particle
 * @param surfaceArea    Surface area of particle
 * @return
 */
template<class NormalForceSpecies>
std::array<double,2> HeatFluidCoupledSpecies<NormalForceSpecies>::f(double liquidVolume_,double temperature_, double mass, double surfaceArea) const {
    // Saturated vapour concentration  (kg/m^3), eq(7)
    // Dependency on temperature is valid between 1 and 200 degree Celsius
    // (it is unphysically decreasing between 0 and 1 degree),
    // extrapolated as linear functions
    double dT = temperature_-273;
    double saturatedVapourConcentration = dT < 1 ? 8.319815774e-3*temperature_/274 :
                                          (((4.844e-9*dT-1.4807e-7)*dT+2.6572e-5)*dT-4.8613e-5)*dT+8.342e-3;

    // Relative activation energy of evaporation, eq(6)
    double equilibriumActivationEnergy=-constants::R*getAmbientTemperature()*log(getAmbientHumidity());
    // \todo what is the difference between the variable X and Y in Azmir?
    double moistureContent = liquidVolume_*getLiquidDensity()/mass;
    double activationEnergy=equilibriumActivationEnergy*(moistureContent-getAmbientEquilibriumMoistureContent());

    // Vapour concentrations at the particle-medium interface (kg/m^3), eq(5)
    // Implemented such that the code does not necessarily fail if temperature is 0
    double interfaceVapourConcentration = temperature_==0?0.0:
                                          exp(-activationEnergy/(constants::R*temperature_))*saturatedVapourConcentration;

    // Drying rate of liquid film mass (kg/s), eq (4)
    double dLiquidMass=-getMassTransferCoefficient()*surfaceArea
                       *(interfaceVapourConcentration-getAmbientVapourConcentration());

    // Drying rate of liquid film volume (kg/s), eq (3)
    double dLiquidVolume = dLiquidMass/getLiquidDensity();

    // Heat of evaporation (J/s), eq (2)
    double heatOfEvaporation = getLatentHeatVaporization()*(1.0+getEvaporationCoefficientA()*exp((getEvaporationCoefficientB()*getLiquidDensity()/mass)*liquidVolume_))*dLiquidMass;

    // Cooling rate of particle (K/s), eq (1)
    double dTemperature = heatOfEvaporation / (mass * this->getHeatCapacity());

    return {dLiquidVolume, dTemperature};
}
#endif
