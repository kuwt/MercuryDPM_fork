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


#include <Species/NormalForceSpecies/MeltableNormalSpecies.h>
#include "Particles/MeltableParticle.h"

/*!
 * \details Copy method. Uses copy constructor to create a copy on the heap. 
 *          Useful for polymorphism.
 * \return pointer to the particle's copy
 */
MeltableParticle* MeltableParticle::copy() const
{
    return new MeltableParticle(*this);
}

/*!
 * \details Returns the name of the object; in this case 'MeltableParticle'.
 * \return The object name.
 */
std::string MeltableParticle::getName() const
{
    return "MeltableParticle";
}

Mdouble MeltableParticle::getParticleProjectedArea ()
{
    return constants::pi * getRadius()*getRadius();
}

/**
 * Compute heat
 */
void MeltableParticle::actionsBeforeTimeStep()
{
    const auto s = meltableSpecies_;
    heat_ = s->getHeatInput(this);
    // add thermal convection
    heat_ += -s->getThermalConvectionCoefficient() * getParticleProjectedArea()
            * (getTemperature() - s->getAmbientTemperature());
    // add thermal radiation
    Mdouble radiationCoefficient = constants::stefanBoltzmanConstant
            * s->getMaterialEmissivity()
            * (getTemperature() + s->getAmbientTemperature())
            * (mathsFunc::square(getTemperature())
            + mathsFunc::square(s->getAmbientTemperature()));
    heat_ += -radiationCoefficient * getParticleProjectedArea() * (getTemperature() - s->getAmbientTemperature());
}

/**
 * Compute heat
 */
void MeltableParticle::actionsAfterTimeStep()
{
    // apply heat
    double heatCapacity = meltableSpecies_->getEffectiveHeatCapacity(getTemperature());
    double timeStep = getHandler()->getDPMBase()->getTimeStep();
    double temperatureStep = heat_*getInvMass()/heatCapacity*timeStep;
    addTemperature(temperatureStep);

//    // add bonding radius
//    for (auto i : getInteractions()) {
//        if (i->getP()==this) {
//            auto m = dynamic_cast<MeltableInteraction*>(i);
//            m->addBondingOverlap()
//        }
//    }

    // thermal expansion
    if (meltableSpecies_->getThermalExpansionCoefficient()!=0)
    {
        setRadius(getRadius()+meltableSpecies_->getThermalExpansionCoefficient()*getRadius()*temperatureStep);
        ///\todo mass not kept!
    }
}

std::string MeltableParticle::getNameVTK(unsigned i) const
{
    if (i==0)
        return "MoltenLayerThickness";
    else if (i==1)
        return "SolidRadius";
    else /*i=2*/
        return "Temperature";
}

std::vector<Mdouble> MeltableParticle::getFieldVTK(unsigned i) const
{
    double solidRadius = getSolidRadius();
    double moltenLayerThickness = getRadius() - solidRadius;
    if (i==0) {
        return std::vector<Mdouble>(1, moltenLayerThickness);
    } else if (i==1) {
        return std::vector<Mdouble>(1, solidRadius);
    } else {
        return std::vector<Mdouble>(1, temperature_);
    }
}
