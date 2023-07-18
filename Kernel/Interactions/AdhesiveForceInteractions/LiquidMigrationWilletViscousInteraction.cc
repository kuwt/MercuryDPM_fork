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


#include "LiquidMigrationWilletViscousInteraction.h"
#include "Species/AdhesiveForceSpecies/LiquidMigrationWilletViscousSpecies.h"
#include "Particles/LiquidFilmParticle.h"
#include "InteractionHandler.h"
#include "DPMBase.h"
#include <iomanip>
#include <fstream>
#include <Interactions/FrictionForceInteractions/SlidingFrictionInteraction.h>

/*!
* \param[in] P
* \param[in] I
* \param[in] timeStamp
*/
LiquidMigrationWilletViscousInteraction::LiquidMigrationWilletViscousInteraction(BaseInteractable* P, BaseInteractable* I,
                                                                                 unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp), LiquidMigrationWilletInteraction(P, I, timeStamp)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "LiquidMigrationWilletViscousInteraction::LiquidMigrationWilletViscousInteraction() finished" << std::endl;
#endif
}

//used for mpi
LiquidMigrationWilletViscousInteraction::LiquidMigrationWilletViscousInteraction()
        : BaseInteraction(), LiquidMigrationWilletInteraction()
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "LiquidMigrationWilletViscousInteraction::LiquidMigrationWilletViscousInteraction() finished" << std::endl;
#endif
}

/*!
 * \param[in] p
 */
LiquidMigrationWilletViscousInteraction::LiquidMigrationWilletViscousInteraction(const LiquidMigrationWilletViscousInteraction& p)
        : BaseInteraction(p), LiquidMigrationWilletInteraction(p)
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "LiquidMigrationWilletViscousInteraction::LiquidMigrationWilletViscousInteraction(const LiquidMigrationWilletViscousInteraction &p finished" << std::endl;
#endif
}

/*!
 *
 */
LiquidMigrationWilletViscousInteraction::~LiquidMigrationWilletViscousInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout << "LiquidMigrationWilletViscousInteraction::~LiquidMigrationWilletViscousInteraction() finished" << std::endl;
#endif
}


/*!
 * This model is a combination of Willet capillary model and a viscous lubrication model, which
 * includes a normal capillary force (same as Willet model), a normal lubrication force and a tangential lubrication force.
 * The description of the viscous lubrication forces can be found here: https://doi.org/10.1016/j.powtec.2019.12.021
 * The capillary force is valid when particles are within rupture distance.
 * The lubrication force is valid when particles are between the limiting distance and the rupture distance.
 * Note: the tangential lubrication force is added to sliding friction force through SlidingFrictionInteraction::addTangentialForce().
 */
void LiquidMigrationWilletViscousInteraction::computeAdhesionForce() {
    LiquidMigrationWilletInteraction::computeAdhesionForce();
    const LiquidMigrationWilletViscousSpecies *species = getSpecies();
    const Mdouble effectiveRadius = getEffectiveRadius();//note here is different from in Willet capillary force. 1/R*=1/RI+1/RP
    if(getOverlap() < 0 && LiquidMigrationWilletInteraction::getWasInContact()){
        Mdouble fdotnl = 0.0;//normal lubrication force in value
        Vec3D fdottl = Vec3D(0.0, 0.0, 0.0);//tangential lubrication force in vector

        if (-getOverlap() > getLimitingDistance()) {
            fdotnl = -6.0 * constants::pi * species->getViscosity() * pow(effectiveRadius, 2.0)
                     * getNormalRelativeVelocity() / (-getOverlap());

            Vec3D tangentialRelativeVelocity = getRelativeVelocity() - getNormal() * getNormalRelativeVelocity();

            fdottl = -6.0 * constants::pi * species->getViscosity() * effectiveRadius * tangentialRelativeVelocity
                     * (8.0 / 15.0 * std::log(effectiveRadius / (-getOverlap())) + 0.9588);

        }
        addForce(getNormal() * fdotnl + fdottl);

        auto slidingFrictionInteraction = dynamic_cast<SlidingFrictionInteraction *>(this);
        slidingFrictionInteraction->addTangentialForce(fdottl);
    }

}

/*!
 * Defines and accesses the minimum distance for the viscous lubrication force to be valid
 */
Mdouble LiquidMigrationWilletViscousInteraction::getLimitingDistance()
{
    return getEffectiveRadius()/100.0;
}

/*!
* \return std::string
*/
std::string LiquidMigrationWilletViscousInteraction::getBaseName() const
{
    return "LiquidMigrationWilletViscous";
}

/*!
 * \return const LiquidMigrationWilletViscousSpecies*
 */
const LiquidMigrationWilletViscousSpecies* LiquidMigrationWilletViscousInteraction::getSpecies() const
{
    return static_cast<const LiquidMigrationWilletViscousSpecies*>(getBaseSpecies()->getAdhesiveForce());
}


