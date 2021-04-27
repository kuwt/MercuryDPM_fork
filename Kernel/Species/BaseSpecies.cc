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

#include "BaseSpecies.h"
#include "SpeciesHandler.h"
#include "DPMBase.h"
#include "Interactions/BaseInteraction.h"
#include<cmath>
#include <Interactions/NormalForceInteractions/LinearPlasticViscoelasticInteraction.h>
#include "Species/NormalForceSpecies/LinearPlasticViscoelasticNormalSpecies.h"

class BaseParticle;

class BaseInteractable;

BaseSpecies::BaseSpecies()
        : BaseObject(), handler_(nullptr), interactionDistance_(0)
{
    logger(DEBUG, "BaseSpecies::BaseSpecies() finished");
}

BaseSpecies::BaseSpecies(BaseNormalForce* normalForce, BaseFrictionForce* frictionForce, BaseAdhesiveForce* adhesiveForce)
        : BaseObject(), normalForce_(normalForce), frictionForce_(frictionForce), adhesiveForce_(adhesiveForce),
        handler_(nullptr), interactionDistance_(0)
{
    normalForce_->setBaseSpecies(this);
    frictionForce_->setBaseSpecies(this);
    adhesiveForce_->setBaseSpecies(this);
    logger(DEBUG, "BaseSpecies::BaseSpecies() finished");
}


/*!
 * \param[in] p the species that is copied
 */
BaseSpecies::BaseSpecies(const BaseSpecies& p)
        : BaseObject(p), handler_(p.handler_), interactionDistance_(p.interactionDistance_)
{
    logger(DEBUG, "BaseSpecies::BaseSpecies(const BaseSpecies &p) finished");
}

BaseSpecies::~BaseSpecies()
{
    ///\todo the BaseSpecies destructor should delete all particles and wall belonging to that species;
    /// however, that will break all codes replacing a species (e.g. Sudeshna and Hao)
    /// so we need a proper way to replace a species
    //for now, I added a particleHandler.clear(); wallHandler.clear() before speciesHandler.clear() in DPMBase::read()

//    //if species gets removed, all particles/ walls of this type need to be removed as well
//    if (getHandler()) {
//        DPMBase* dpm = getHandler()->getDPMBase();
//        for (const auto o : dpm->wallHandler) {
//            if (o->getSpecies()==this) {
//                dpm->wallHandler.removeObject(o->getIndex());
//            }
//        }
//        for (const auto o : dpm->particleHandler) {
//            if (o->getSpecies()==this) {
//                dpm->particleHandler.removeObject(o->getIndex());
//            }
//        }
//    }
    logger(DEBUG, "BaseSpecies::~BaseSpecies() finished");
}

/*!
 * \param[in] the pointer to the handler to which this species belongs.
 */
void BaseSpecies::setHandler(SpeciesHandler* const handler)
{
    handler_ = handler;
}

/*!
 * \return the pointer to the handler to which this species belongs.
 */
SpeciesHandler* BaseSpecies::getHandler() const
{
    return handler_;
}

/*!
 * \detail Returns the harmonic mean of two variables.
 * This function is used to define default mixed species.
 * \param[in] a,b The two variables you want to average
 * \return The harmonic mean of a and b, \f$\frac{2}{1/a+1/b}\f$
 */
Mdouble BaseSpecies::average(Mdouble a, Mdouble b)
{
    //the second algorithm seems to have a better accuracy, at least for the case average(2e5,2e5)
    //return (a + b) != 0.0 ? (2. * (a * b) / (a + b)) : 0;
    return (a==0||b==0) ? 0.0 : (2. / (1.0 / a + 1.0 / b));
}

/*! 
 * \detail Returns the harmonic mean of two variables, returning inf if either is inf.
 */
Mdouble BaseSpecies::averageInf(Mdouble a, Mdouble b)
{
    if (a == constants::inf || b == constants::inf)
        return constants::inf;
    else
        return average(a, b);
}

/*!
 * \param[out] os output stream (typically the restart file)
 */
void BaseSpecies::write(std::ostream& os) const
{
    //BaseObject::write(os);
    if (normalForce_->getConstantRestitution()) os << " constantRestitution " << normalForce_->getConstantRestitution();
}

/*!
 * \param[in] is input stream (typically the restart file)
 */
void BaseSpecies::read(std::istream& is)
{
    //BaseObject::read(is);
    bool constantRestitution;
    if (helpers::readOptionalVariable(is, "constantRestitution", constantRestitution)) {
        normalForce_->setConstantRestitution(constantRestitution);
    }
}

/**
 * Sets BaseSpecies::interactionDistance_.
 * This function should not be called by the user, only by functions in classes derived from BaseSpecies (in particular, the adhesive-force species).
 * This function gets called every time a variable is set on which the interaction distance depends.
 * See for example LiquidBridgeWilletSpecies::setLiquidBridgeVolume.
 * @param interactionDistance
 */
void BaseSpecies::setInteractionDistance(Mdouble interactionDistance) {
    interactionDistance_ = interactionDistance;
    
    SpeciesHandler* handler = getHandler();
    if (handler == nullptr) return;
    
    for (auto mixedSpecies : handler->getMixedObjects()) {
        if (mixedSpecies == this) {
            // get the two particlespecies id's
            unsigned mixedId = mixedSpecies->getIndex();
            unsigned maxId= 1;
            unsigned maxMixedId = (maxId * (maxId + 1)) / 2;
            while (maxMixedId<mixedId) {
                ++maxId;
                maxMixedId = (maxId * (maxId + 1)) / 2;
            }
            unsigned minId = (mixedId + maxId) - maxMixedId;
            handler->getObject(minId)->setMaxInteractionDistance(interactionDistance);
            handler->getObject(maxId)->setMaxInteractionDistance(interactionDistance);
            //logger(INFO,"setInteractionDistance(%) mixed % handler %",interactionDistance, getIndex(), getHandler());
            return;
        }
    }
    handler->getObject(getId())->setMaxInteractionDistance(interactionDistance);
    
    //logger(INFO,"setInteractionDistance(%) species % handler %",interactionDistance, getIndex(), getHandler());
}
