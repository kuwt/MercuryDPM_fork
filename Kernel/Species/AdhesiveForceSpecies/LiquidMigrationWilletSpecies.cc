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

#include "LiquidMigrationWilletSpecies.h"
#include <Logger.h>

LiquidMigrationWilletSpecies::LiquidMigrationWilletSpecies()
{
    liquidBridgeVolumeMax_ = 0.0; //std::numeric_limits<double>::infinity();
    distributionCoefficient_ = 1.0;
    maxInteractionDistance_ = 0.0;
    surfaceTension_ = 0.0;
    contactAngle_ = 0.0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LiquidMigrationWilletSpecies::LiquidMigrationWilletSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[in] the species that is copied
 */
LiquidMigrationWilletSpecies::LiquidMigrationWilletSpecies(const LiquidMigrationWilletSpecies& s)
{
    liquidBridgeVolumeMax_ = s.liquidBridgeVolumeMax_;
    distributionCoefficient_ = s.distributionCoefficient_;
    maxInteractionDistance_ = s.maxInteractionDistance_;
    surfaceTension_ = s.surfaceTension_;
    contactAngle_ = s.contactAngle_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LiquidMigrationWilletSpecies::LiquidMigrationWilletSpecies(const LiquidMigrationWilletSpecies &p) finished"<<std::endl;
#endif
}

LiquidMigrationWilletSpecies::~LiquidMigrationWilletSpecies()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"LiquidMigrationWilletSpecies::~LiquidMigrationWilletSpecies() finished"<<std::endl;
#endif
}

/*!
 * \param[out] output stream (typically the restart file)
 */
void LiquidMigrationWilletSpecies::write(std::ostream& os) const
{
    os << " liquidBridgeVolume " << liquidBridgeVolumeMax_;
    os << " distributionCoefficient " << distributionCoefficient_;
    os << " surfaceTension " << surfaceTension_;
    os << " contactAngle " << contactAngle_;
}

/*!
 * \param[in] input stream (typically the restart file)
 */
void LiquidMigrationWilletSpecies::read(std::istream& is)
{
    std::string dummy;
    is >> dummy >> liquidBridgeVolumeMax_;
    is >> dummy >> distributionCoefficient_;
    is >> dummy >> surfaceTension_;
    is >> dummy >> contactAngle_;
}

/*!
 * \return a string containing the name of the species (minus the word "Species")
 */
std::string LiquidMigrationWilletSpecies::getBaseName() const
{
    return "LiquidMigrationWillet";
}

/*!
 * \details For all parameters we assume that the harmonic mean of the parameters of the 
 * original two species is a sensible default.
 * \param[in] S,T the two species whose properties are mixed to create the new species
 */
void LiquidMigrationWilletSpecies::mix(LiquidMigrationWilletSpecies* const S, LiquidMigrationWilletSpecies* const T)
{
    setLiquidBridgeVolumeMax(average(S->getLiquidBridgeVolumeMax(), T->getLiquidBridgeVolumeMax()));
    distributionCoefficient_ = average(S->getDistributionCoefficient(), T->getDistributionCoefficient());
    surfaceTension_ = average(S->getSurfaceTension(), T->getSurfaceTension());
    contactAngle_ = average(S->getContactAngle(), T->getContactAngle());
}

/*!
 * \return the maximum separation distance between particles below which the adhesive force is active.
 */
Mdouble LiquidMigrationWilletSpecies::getInteractionDistance() const
{
    return maxInteractionDistance_;
}

/*!
 * \param[in] liquidBridgeVolume the volume of the liquid bridge.
 */
void LiquidMigrationWilletSpecies::setLiquidBridgeVolumeMax(Mdouble liquidBridgeVolumeMax)
{
    if (liquidBridgeVolumeMax >= 0)
    {
        liquidBridgeVolumeMax_ = liquidBridgeVolumeMax;
        maxInteractionDistance_ = (1.0 + 0.5 * contactAngle_) * cbrt(liquidBridgeVolumeMax_);
    }
    else
    {
        std::cerr << "Error in setLiquidBridgeVolumeMax: liquidBridgeVolumeMax=" << liquidBridgeVolumeMax << std::endl;
        exit(-1);
    }
}

/*!
 * \return the volume of the liquid bridge.
 */
Mdouble LiquidMigrationWilletSpecies::getLiquidBridgeVolumeMax() const
{
    return liquidBridgeVolumeMax_;
}

/*!
 * \param[in] distributionCoefficient the distribution coefficient of the liquid.
 */
void LiquidMigrationWilletSpecies::setDistributionCoefficient(Mdouble distributionCoefficient)
{
    if (distributionCoefficient >= 0 && distributionCoefficient <= 1.0)
        distributionCoefficient_ = distributionCoefficient;
    else
    {
        std::cerr << "Error in setDistributionCoefficient" << std::endl;
        exit(-1);
    }
}

/*!
 * \return the distribution coefficient of the liquid.
 */
Mdouble LiquidMigrationWilletSpecies::getDistributionCoefficient() const
{
    return distributionCoefficient_;
}


/*!
 * \param[in] surfaceTension the surface tension of the liquid.
 */
void LiquidMigrationWilletSpecies::setSurfaceTension(Mdouble surfaceTension)
{
    if (surfaceTension >= 0)
        surfaceTension_ = surfaceTension;
    else
    {
        std::cerr << "Error in setSurfaceTension" << std::endl;
        exit(-1);
    }
}

/*!
 * \return the surface tension of the liquid.
 */
Mdouble LiquidMigrationWilletSpecies::getSurfaceTension() const
{
    return surfaceTension_;
}

/*!
 * \param[in] contactAngle the contact angle between particle and liquid bridge surface.
 */
void LiquidMigrationWilletSpecies::setContactAngle(Mdouble contactAngle)
{
    if (contactAngle >= 0)
    {
        contactAngle_ = contactAngle;
        maxInteractionDistance_ = (1.0 + 0.5 * contactAngle_) * cbrt(liquidBridgeVolumeMax_);
    }
    else
    {
        std::cerr << "Error in setContactAngle" << std::endl;
        exit(-1);
    }
}

/*!
 * \return the contact angle between particle and liquid bridge surface.
 */
Mdouble LiquidMigrationWilletSpecies::getContactAngle() const
{
    return contactAngle_;
}

