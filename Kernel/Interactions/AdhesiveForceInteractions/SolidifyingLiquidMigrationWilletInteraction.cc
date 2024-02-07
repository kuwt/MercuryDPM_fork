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


#include "SolidifyingLiquidMigrationWilletInteraction.h"
#include "Species/AdhesiveForceSpecies/SolidifyingLiquidMigrationWilletSpecies.h"
#include "Particles/BaseParticle.h"
#include "InteractionHandler.h"
#include <iomanip>
#include <fstream>

/*!
 * \param[in] P
 * \param[in] I
 * \param[in] timeStamp
 */
SolidifyingLiquidMigrationWilletInteraction::SolidifyingLiquidMigrationWilletInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp)
        : LiquidMigrationWilletInteraction(P, I, timeStamp)
{
    bonded_ = false;
    solidVolume_ = 0;
    bondForce_ = 0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"SolidifyingLiquidMigrationWilletInteraction::SolidifyingLiquidMigrationWilletInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] p
 */
SolidifyingLiquidMigrationWilletInteraction::SolidifyingLiquidMigrationWilletInteraction(const SolidifyingLiquidMigrationWilletInteraction& p)
        : LiquidMigrationWilletInteraction(p)
{
    ///\todo tw check if the parameters are valid when inserting the species into the handler
    bonded_ = p.bonded_;
    solidVolume_ = p.solidVolume_;
    bondForce_ = p.bondForce_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"SolidifyingLiquidMigrationWilletInteraction::SolidifyingLiquidMigrationWilletInteraction(const SolidifyingLiquidMigrationWilletInteraction &p finished"<<std::endl;
#endif
}

SolidifyingLiquidMigrationWilletInteraction::SolidifyingLiquidMigrationWilletInteraction()
{
#ifdef MERCURY_USE_MPI
    logger(FATAL,"ChargedSolidifyingLiquidMigrationWilletInteractions are currently not implemented in parallel MercuryDPM");
#endif
}

/*!
 *
 */
SolidifyingLiquidMigrationWilletInteraction::~SolidifyingLiquidMigrationWilletInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"SolidifyingLiquidMigrationWilletInteraction::~SolidifyingLiquidMigrationWilletInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in,out] os
 */
void SolidifyingLiquidMigrationWilletInteraction::write(std::ostream& os) const
{
    LiquidMigrationWilletInteraction::write(os);
    os << " bonded " << bonded_;
    os << " solidVolume " << solidVolume_;
    os << " bondForce " << bondForce_;
}

/*!
 * \param[in,out] is
 */
void SolidifyingLiquidMigrationWilletInteraction::read(std::istream& is)
{
    LiquidMigrationWilletInteraction::read(is);
    std::string dummy;
    is >> dummy >> bonded_;
    is >> dummy >> solidVolume_;
    is >> dummy >> bondForce_;
}

/*!
 * \details Uses the most basic adhesion contact model.
 */
void SolidifyingLiquidMigrationWilletInteraction::computeAdhesionForce()
{
    if (bonded_) {
        addForce(getNormal() * (-getBondForce() - getSpecies()->getBondDissipation() * getNormalRelativeVelocity()));
    } else {
        if (getLiquidBridgeVolume()>0 && solidVolume_ >= getLiquidBridgeVolume()) {
            //logger(WARN,"Solidifying contact, V=%", solidVolume_);
            bondInPlace();
        }
    }
}

/*!
 * \return const pointer of SolidifyingLiquidMigrationWilletSpecies*
 */
const SolidifyingLiquidMigrationWilletSpecies* SolidifyingLiquidMigrationWilletInteraction::getSpecies() const
{
    return static_cast<const SolidifyingLiquidMigrationWilletSpecies*>(getBaseSpecies()->getAdhesiveForce());
;
}

/*!
 * \return std::string
 */
std::string SolidifyingLiquidMigrationWilletInteraction::getBaseName() const
{
    return "Solidifying" + LiquidMigrationWilletInteraction::getBaseName();
}

/*!
 * \details Elastic (=potential) energy is defined as the energy gained by separating two interactables.
 * As it costs energy to separate adhesive interactables, the elastic energy is negative.
 * \return the elastic energy stored in the adhesive spring. 
 */
Mdouble SolidifyingLiquidMigrationWilletInteraction::getElasticEnergy() const
{
    if (!bonded_)
        return 0.0;
    else
        return -getSpecies()->getBondForceMax() * getOverlap();
}

bool SolidifyingLiquidMigrationWilletInteraction::getBonded() const
{
    return bonded_;
}

void SolidifyingLiquidMigrationWilletInteraction::setBonded(bool bonded)
{
    bonded_ = bonded;
}


void SolidifyingLiquidMigrationWilletInteraction::bond()
{
    bondForce_=getSpecies()->getBondForceMax();
    bonded_ = true;
}

void SolidifyingLiquidMigrationWilletInteraction::bondInPlace()
{
    bondForce_= getForce().getLength();
    bonded_ = true;
}

void SolidifyingLiquidMigrationWilletInteraction::unbond()
{
    bonded_ = false;
}

unsigned SolidifyingLiquidMigrationWilletInteraction::getNumberOfFieldsVTK() const
{
    return 3;
}

std::string SolidifyingLiquidMigrationWilletInteraction::getTypeVTK(unsigned i) const
{
    return "Float32";
}

std::string SolidifyingLiquidMigrationWilletInteraction::getNameVTK(unsigned i) const
{
    if (i==0) {
        return "liquidBridgeRadius";
    } else if (i==1) {
        return "solidBridgeRadius";
    } else {
        return "solidified";
    }

}

std::vector<Mdouble> SolidifyingLiquidMigrationWilletInteraction::getFieldVTK(unsigned i) const
{
    if (i==0) {
        return std::vector<Mdouble>(1, cbrt(getLiquidBridgeVolume()));
    } else if (i==1) {
        return std::vector<Mdouble>(1, cbrt(solidVolume_));
    } else {
        return std::vector<Mdouble>(1, getBonded());
    }
}

void SolidifyingLiquidMigrationWilletInteraction::setLiquidBridgeVolume(Mdouble liquidBridgeVolume) {
    LiquidMigrationWilletInteraction::setLiquidBridgeVolume(liquidBridgeVolume);
    solidVolume_ = liquidBridgeVolume * getSpecies()->getSolidFraction();
}

void SolidifyingLiquidMigrationWilletInteraction::addLiquidBridgeVolume(Mdouble liquidBridgeVolume) {
    LiquidMigrationWilletInteraction::addLiquidBridgeVolume(liquidBridgeVolume);
    solidVolume_ += liquidBridgeVolume * getSpecies()->getSolidFraction();
}