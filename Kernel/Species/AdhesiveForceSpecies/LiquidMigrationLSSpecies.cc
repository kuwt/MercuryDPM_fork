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

#include "LiquidMigrationLSSpecies.h"
#include "Logger.h"

LiquidMigrationLSSpecies::LiquidMigrationLSSpecies()
{
    liquidBridgeVolumeMax_ = 0.0; //std::numeric_limits<double>::infinity();
    liquidBridgeVolumeMin_ = 0.0;
    distributionCoefficient_ = 1.0;
    surfaceTension_ = 0.0;
    contactAngle_ = 0.0;
    viscosity_ = 0.0;

#ifdef DEBUG_CONSTRUCTOR
    logger(INFO, "LiquidMigrationLSSpecies::LiquidMigrationLSSpecies() finished");
#endif
}

/*!
 * \param[in] the species that is copied
 */
LiquidMigrationLSSpecies::LiquidMigrationLSSpecies(const LiquidMigrationLSSpecies& s)
{
    liquidBridgeVolumeMax_ = s.liquidBridgeVolumeMax_;
    liquidBridgeVolumeMin_ = s.liquidBridgeVolumeMin_;
    distributionCoefficient_ = s.distributionCoefficient_;
    surfaceTension_ = s.surfaceTension_;
    contactAngle_ = s.contactAngle_;
    viscosity_ = s.viscosity_;

#ifdef DEBUG_CONSTRUCTOR
    logger(INFO, "LiquidMigrationLSSpecies::LiquidMigrationLSSpecies(const LiquidMigrationLSSpecies &p) finished");
#endif
}

LiquidMigrationLSSpecies::~LiquidMigrationLSSpecies()
{
#ifdef DEBUG_DESTRUCTOR
    logger(INFO, "LiquidMigrationLSSpecies::~LiquidMigrationLSSpecies() finished");
#endif
}

/*!
 * \param[out] output stream (typically the restart file)
 */
void LiquidMigrationLSSpecies::write(std::ostream& os) const
{
    os << " liquidBridgeVolumeMax " << liquidBridgeVolumeMax_;
    if (liquidBridgeVolumeMin_) os << " liquidBridgeVolumeMin " << liquidBridgeVolumeMin_;
    os << " distributionCoefficient " << distributionCoefficient_;
    os << " surfaceTension " << surfaceTension_;
    os << " contactAngle " << contactAngle_;
    os << " viscosity_ " << viscosity_;
}

/*!
 * \param[in] input stream (typically the restart file)
 */
void LiquidMigrationLSSpecies::read(std::istream& is)
{
    std::string dummy;
    is >> dummy >> liquidBridgeVolumeMax_;
    helpers::readOptionalVariable(is,"liquidBridgeVolumeMin",liquidBridgeVolumeMin_);
    is >> dummy >> distributionCoefficient_;
    is >> dummy >> surfaceTension_;
    is >> dummy >> contactAngle_;
    is >> dummy >> viscosity_;
}

/*!
 * \return a string containing the name of the species (minus the word "Species")
 */
std::string LiquidMigrationLSSpecies::getBaseName() const
{
    return "LiquidMigrationLS";
}

/*!
 * \details For all parameters we assume that the harmonic mean of the parameters of the 
 * original two species is a sensible default.
 * \param[in] S,T the two species whose properties are mixed to create the new species
 */
void LiquidMigrationLSSpecies::mix(LiquidMigrationLSSpecies* const S, LiquidMigrationLSSpecies* const T)
{
    setLiquidBridgeVolumeMax(BaseSpecies::average(S->getLiquidBridgeVolumeMax(), T->getLiquidBridgeVolumeMax()));
    setLiquidBridgeVolumeMin(BaseSpecies::average(S->getLiquidBridgeVolumeMin(), T->getLiquidBridgeVolumeMin()));
    distributionCoefficient_ = BaseSpecies::average(S->getDistributionCoefficient(), T->getDistributionCoefficient());
    surfaceTension_ = BaseSpecies::average(S->getSurfaceTension(), T->getSurfaceTension());
    contactAngle_ = BaseSpecies::average(S->getContactAngle(), T->getContactAngle());
    viscosity_ = BaseSpecies::average(S->getViscosity(), T->getViscosity());
}

///\return the maximum separation distance below which adhesive forces can occur (needed for contact detection)
void LiquidMigrationLSSpecies::setInteractionDistance()
{
    getBaseSpecies()->setInteractionDistance((1.0 + 0.5 * contactAngle_) * (cbrt(liquidBridgeVolumeMax_) + 0.1 * pow(liquidBridgeVolumeMax_, 2.0/3.0)));
}

/*!
 * \param[in] liquidBridgeVolume the maximum volume of the liquid bridge.
 */
void LiquidMigrationLSSpecies::setLiquidBridgeVolumeMax(Mdouble liquidBridgeVolumeMax)
{
    logger.assert_always(liquidBridgeVolumeMax>=0,
                         "Error in setLiquidBridgeVolumeMax: liquidBridgeVolumeMax=%", liquidBridgeVolumeMax);
    liquidBridgeVolumeMax_ = liquidBridgeVolumeMax;
    setInteractionDistance();
}

/*!
 * \param[in] liquidBridgeVolume the minimum volume of the liquid bridge, the criterion value of the summation of the liquid film volumes for the liquid bridge forming.
 */
void LiquidMigrationLSSpecies::setLiquidBridgeVolumeMin(Mdouble liquidBridgeVolumeMin)
{
    logger.assert_always(liquidBridgeVolumeMin>=0,
            "Error in setLiquidBridgeVolumeMin: liquidBridgeVolumeMin=%", liquidBridgeVolumeMin);
    liquidBridgeVolumeMin_ = liquidBridgeVolumeMin;
}

/*!
 * \return the maximum volume of the liquid bridge.
 */
Mdouble LiquidMigrationLSSpecies::getLiquidBridgeVolumeMax() const
{
    return liquidBridgeVolumeMax_;
}

/*!
 * \return the minimum volume of the liquid bridge.
 */
Mdouble LiquidMigrationLSSpecies::getLiquidBridgeVolumeMin() const
{
    return liquidBridgeVolumeMin_;
}

/*!
 * \param[in] distributionCoefficient the distribution coefficient of the liquid.
 */
void LiquidMigrationLSSpecies::setDistributionCoefficient(Mdouble distributionCoefficient)
{
    if (distributionCoefficient >= 0 && distributionCoefficient <= 1.0)
        distributionCoefficient_ = distributionCoefficient;
    else
    {
        logger(ERROR, "Error in setDistributionCoefficient, DistributionCoefficient should be between 0 and 1");
        exit(-1);
    }
}

/*!
 * \return the distribution coefficient of the liquid.
 */
Mdouble LiquidMigrationLSSpecies::getDistributionCoefficient() const
{
    return distributionCoefficient_;
}


/*!
 * \param[in] surfaceTension the surface tension of the liquid.
 */
void LiquidMigrationLSSpecies::setSurfaceTension(Mdouble surfaceTension)
{
    if (surfaceTension >= 0)
        surfaceTension_ = surfaceTension;
    else
    {
        logger(ERROR, "Error in setSurfaceTension, surfaceTension should be >= 0");
        exit(-1);
    }
}

/*!
 * \return the surface tension of the liquid.
 */
Mdouble LiquidMigrationLSSpecies::getSurfaceTension() const
{
    return surfaceTension_;
}

/*!
 * \param[in] contactAngle the contact angle between particle and liquid bridge surface.
 */
void LiquidMigrationLSSpecies::setContactAngle(Mdouble contactAngle)
{
    if (contactAngle >= 0)
    {
        contactAngle_ = contactAngle;
        setInteractionDistance();
    }
    else
    {
        logger(ERROR, "Error in setContactAngle, contactAngle should be >= 0");
        exit(-1);
    }
}

/*!
 * \return the contact angle between particle and liquid bridge surface.
 */
Mdouble LiquidMigrationLSSpecies::getContactAngle() const
{
    return contactAngle_;
}

/*!
 * \param[in] viscosity the viscosity of the liquid.
 */
void LiquidMigrationLSSpecies::setViscosity(Mdouble viscosity)
{
    if (viscosity >= 0)
    {
        viscosity_ = viscosity;
    }
    else
    {
        logger(ERROR, "Error in setViscosity, viscosity should be >= 0");
        exit(-1);
    }
}

/*!
 * \return the the viscosity of the liquid.
 */
Mdouble LiquidMigrationLSSpecies::getViscosity() const
{
    return viscosity_;
}


