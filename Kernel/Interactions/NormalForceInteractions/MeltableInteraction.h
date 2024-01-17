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

#ifndef MERCURY_MELTABLEINTERACTION_H
#define MERCURY_MELTABLEINTERACTION_H

#include "Interactions/BaseInteraction.h"
#include "ParticleHandler.h"
#include "InteractionHandler.h"

class BaseParticle;
class MeltableParticle;
class BaseInteractable;
class MeltableNormalSpecies;

class MeltableInteraction : public virtual BaseInteraction {
public:
    MeltableInteraction(BaseInteractable *P, BaseInteractable *I, unsigned timeStamp)
            : BaseInteraction(P, I, timeStamp) {}

    MeltableInteraction(const MeltableInteraction &p) = default;

    MeltableInteraction() = default;

    void read(std::istream &is) override;

    void write(std::ostream &os) const override;

    void computeNormalForce();

    void actionsAfterTimeStep() override;

    Mdouble getElasticEnergy() const override;

    std::string getBaseName() const {
        return "Meltable";
    }

    /*!
     * Returns the type of normal species required for this kind of interaction
     */
    const MeltableNormalSpecies* getMeltableNormalSpecies() const;

    Mdouble getBondingOverlap() const
    {
        return bondingOverlap_;
    }

    void setBondingOverlap(Mdouble bondingOverlap)
    {
        bondingOverlap_ = bondingOverlap;
    }

    void addBondingOverlap(Mdouble bondingOverlap)
    {
        bondingOverlap_ += bondingOverlap;
    }

    Mdouble getNeckRadius() const {
        return std::sqrt(neckRadius2_);
    };

    void addNeckRadius2(Mdouble diffNeckRadius2) {
        neckRadius2_ = std::max(neckRadius2_+diffNeckRadius2,0.0);
    }

    Mdouble getOverlapGrowthRate() const  {
        return Vec3D::dot(getI()->getVelocity()-getP()->getVelocity(),getNormal());
    }

private:
    Mdouble bondingOverlap_ = 0.0;

    Mdouble neckRadius2_ = constants::NaN;
};

#endif
