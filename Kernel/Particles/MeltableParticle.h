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

#ifndef MeltableParticle_H
#define MeltableParticle_H

#include "ThermalParticle.h"
#include "Species/NormalForceSpecies/MeltableNormalSpecies.h"

/*!
 * \class MeltableParticle
 * \brief
 */
class MeltableParticle : public ThermalParticle {
public:
    /*!
     * \brief Basic Particle constructor, creates a particle at (0,0,0) with radius, mass and inertia equal to 1
     */
    MeltableParticle() = default;

    /*!
     * \brief Particle copy constructor, which accepts as input a reference to a Particle. It creates a copy of this Particle and all it's information. Usually it is better to use the copy() function for polymorphism.
     */
    MeltableParticle(const MeltableParticle &p) = default;

    /*!
     * \brief Particle copy method. It calls to copy constructor of this Particle, useful for polymorphism
     */
    MeltableParticle *copy() const override;

    //void write(std::ostream &os) const override;

    std::string getName() const override;

    Mdouble getParticleProjectedArea();

    const MeltableNormalSpecies* getMeltableSpecies() const {
        return meltableSpecies_;
    }

    void actionsBeforeTimeStep() override;

    void actionsAfterTimeStep() override;

    unsigned getNumberOfFieldsVTK() const override { return 3; }

    std::string getTypeVTK(unsigned i) const override { return "Float32"; }

    std::string getNameVTK(unsigned i) const override;

    std::vector<Mdouble> getFieldVTK(unsigned i) const override;

    void setSpecies(const ParticleSpecies* species) override {
        BaseParticle::setSpecies(species);
        meltableSpecies_ = dynamic_cast<const MeltableNormalSpecies*>(species);
    }

    void addHeat(double heat) {
        heat_ += heat;
    }

    double getMeltRate(double solidRadius) const {
        return heat_/(4.0*constants::pi*solidRadius*solidRadius*getSpecies()->getDensity()*getMeltableSpecies()->getEffectiveLatentHeat());
    }

    double getSolidRadius() const {
        return meltableSpecies_->getRelativeSolidRadius(temperature_)*getRadius();
    }

    double getMoltenLayerThickness() const {
        return getRadius() - getSolidRadius();
    }

private:
    double heat_ = 0;
    const MeltableNormalSpecies* meltableSpecies_ = nullptr;
};

#endif
