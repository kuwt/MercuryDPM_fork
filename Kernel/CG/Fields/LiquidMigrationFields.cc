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

#include "LiquidMigrationFields.h"
#include <Particles/BaseParticle.h>
#include <Particles/LiquidFilmParticle.h>
#include <Interactions/AdhesiveForceInteractions/LiquidMigrationWilletInteraction.h>

namespace CGFields
{

LiquidMigrationFields::LiquidMigrationFields()
{
    setZero();
#ifdef DEBUG_CONSTRUCTOR
    std::cerr << "LiquidMigrationFields::LiquidMigrationFields() finished" << std::endl;
#endif
}

/*!
 * \param[out] os the ostream into which the data is written.
 */
void LiquidMigrationFields::writeNames(std::ostream& os, const unsigned countVariables)
{
    os << countVariables + 1 << ":liquidBridgeVolume ";
    os << countVariables + 2 << ":liquidFilmVolume ";
}

/*!
 * \param[out] os the ostream into which the data is written.
 */
void LiquidMigrationFields::write(std::ostream& os) const
{
    os << liquidBridgeVolume_;
    os << " " << liquidFilmVolume_;
}

/*!
 * \param[out] os the ostream into which the data is written.
 */
void LiquidMigrationFields::output(std::ostream& os) const
{
    os << "liquidBridgeVolume " << liquidBridgeVolume_;
    os << " liquidFilmVolume " << liquidFilmVolume_;
}

void LiquidMigrationFields::setZero()
{
    liquidBridgeVolume_ = 0.0;
    liquidFilmVolume_ = 0.0;
}

/*!
 * \return a CGField containing the square of the values in the current object
 */
LiquidMigrationFields LiquidMigrationFields::getSquared() const
{
    LiquidMigrationFields P;
    P.liquidBridgeVolume_ = mathsFunc::square(liquidBridgeVolume_);
    P.liquidFilmVolume_ = mathsFunc::square(liquidFilmVolume_);
    return P;
}

/*!
 * \param[in] P the CGField that has to be copied
 * \return the CGField into which the values are copied
 */
LiquidMigrationFields& LiquidMigrationFields::operator=(const LiquidMigrationFields& P)
= default;

/*!
* \param[in] P the CGField that has to be added
* \return the CGField to which the values are added
*/
LiquidMigrationFields& LiquidMigrationFields::operator+=(const LiquidMigrationFields& P)
{
    liquidBridgeVolume_ += P.liquidBridgeVolume_;
    liquidFilmVolume_ += P.liquidFilmVolume_;
    return *this;
}

/*!
 * \param[in] P the CGField that has to be subtracted
 * \return the CGField from which the values are subtracted
 */
LiquidMigrationFields& LiquidMigrationFields::operator-=(const LiquidMigrationFields& P)
{
    liquidBridgeVolume_ -= P.liquidBridgeVolume_;
    liquidFilmVolume_ -= P.liquidFilmVolume_;
    return *this;
}

/*!
 * \param[in] a the scalar that we multiply with
 * \return the CGField  to which the multiplied values are written
 */
LiquidMigrationFields LiquidMigrationFields::operator*(const Mdouble a) const
{
    LiquidMigrationFields p;
    p.liquidBridgeVolume_ = liquidBridgeVolume_ * a;
    p.liquidFilmVolume_ = liquidFilmVolume_ * a;
    return p;
}

/*!
 * \param[in] a the scalar that we divide by
 * \return the CGField to which the divided values are written
 */
LiquidMigrationFields& LiquidMigrationFields::operator/=(const Mdouble a)
{
    liquidBridgeVolume_ /= a;
    liquidFilmVolume_ /= a;
    return *this;
}

/*!
 * \param[in] phi the value of the cg function at the current CGPoint
 * \param[in] p the particle which is used in the cg function
 */
void LiquidMigrationFields::addParticleStatistics(Mdouble phi, const LiquidMigrationFields& currentInteraction)
{
    liquidFilmVolume_ += currentInteraction.getLiquidFilmVolume() * phi;
}

void LiquidMigrationFields::addParticleDifferentialStatistics(Vec3D& dphi,
                                                              const LiquidMigrationFields& currentInteraction)
{
}

/*!
 * \param[in] psi the value of the line integral from C to P at the current CGPoint
 * \param[in] c the contact which is used in the line integral
 */
void LiquidMigrationFields::addInteractionStatistics(Mdouble psi, const LiquidMigrationFields& currentInteraction)
{
    liquidBridgeVolume_ += currentInteraction.getLiquidBridgeVolume() * psi;
}

/*!
 * \param[in] phi the value of the cg function for the contact point of c and
 * the current CGPoint
 * \param[in] c the interaction which is used in the cg function
 */
void LiquidMigrationFields::addContactPointStatistics(Mdouble phi UNUSED,
                                                      const LiquidMigrationFields& currentInteraction UNUSED)
{
}

bool LiquidMigrationFields::doInteractionStatistics()
{
    return true;
}

void LiquidMigrationFields::setFields(const BaseInteraction& c, IntegralType type)
{
    auto l = dynamic_cast<const LiquidMigrationWilletInteraction*>(&c);
    logger.assert(l != nullptr,
                  "LiquidMigrationFields::addParticleStatistics: "
                  "interaction type should be LiquidMigrationWilletInteraction");
    liquidBridgeVolume_ = l->getLiquidBridgeVolume();
}
//    if (type==IntegralType::CONTACT_TO_P)
//    {
//        auto p = dynamic_cast<const LiquidFilmParticle*>(c.getI());
//        if (p!=nullptr) {
//            liquidFilmVolume_ = p->getLiquidVolume();
//        } else {
//            liquidBridgeVolume_ = 0;
//        }
//    } else if (type==IntegralType::I_TO_CONTACT)
//    {
//        auto p = dynamic_cast<const LiquidFilmParticle*>(c.getP());
//        if (p!=nullptr) {
//            liquidFilmVolume_ = p->getLiquidVolume();
//        } else {
//            liquidBridgeVolume_ = 0;
//        }
//    }


void LiquidMigrationFields::setFields(const BaseParticle& p)
{
    auto l = dynamic_cast<const LiquidFilmParticle*>(&p);
    logger.assert(l != nullptr,
                  "LiquidMigrationFields::addParticleStatistics: particle type should be LiquidFilmParticle");
    liquidFilmVolume_ = l->getLiquidVolume();
}

void LiquidMigrationFields::setCylindricalFields(const BaseInteraction& c, IntegralType type)
{
    setFields(c, type);
}

void LiquidMigrationFields::setCylindricalFields(const BaseParticle& p)
{
    setFields(p);
}
    
    
}
