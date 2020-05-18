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

#include <Species/LinearViscoelasticSpecies.h>
#include "InteractionHandler.h"
#include "Particles/ThermalParticle.h"
#include "Interactions/BaseInteraction.h"
#include "Species/ParticleSpecies.h"
#include "ParticleHandler.h"
#include "DPMBase.h"

/*!
 * \details default constructor, creates an Particle at (0,0,0) with radius, 
 * mass and heatinertia equal to 1
 */
ThermalParticle::ThermalParticle()
{
    temperature_ = 0;
    //temperatureDependentDensity_
}

/*!
 * \details Constructor that copies most of the properties of the given particle.
 *          Please note that not everything is copied, for example the position 
 *          in the HGrid is not determined yet by the end of this constructor. 
 *          It also does not copy the interactions and the pointer to the handler
 *          that handles this particle. Use with care.
 * \param[in,out] p  Reference to the ThermalParticle this one should become a copy of.
 */
ThermalParticle::ThermalParticle(const ThermalParticle& p)
        : BaseParticle(p)
{
    temperature_ = p.temperature_;
    timeDependentTemperature_ = p.timeDependentTemperature_;
}

/*!
 * \details Destructor. It asks the ParticleHandler to check if this was the 
 *          smallest or largest particle and adjust itself accordingly.
 */
ThermalParticle::~ThermalParticle()
= default;

/*!
 * \details Copy method. Uses copy constructor to create a copy on the heap. 
 *          Useful for polymorphism.
 * \return pointer to the particle's copy
 */
ThermalParticle* ThermalParticle::copy() const
{
    return new ThermalParticle(*this);
}

/*!
 * \details ThermalParticle print method, which accepts an os std::ostream as 
 *          input. It prints human readable ThermalParticle information to the 
 *          std::ostream.
 * \param[in,out] os    stream to which the info is written
 */
void ThermalParticle::write(std::ostream& os) const
{
    BaseParticle::write(os);
    os << " temperature " << temperature_;
}

/*!
 * \details Returns the name of the object; in this case 'ThermalParticle'.
 * \return The object name.
 */
std::string ThermalParticle::getName() const
{
    return "ThermalParticle";
}

//todo Does mass and interaction radius change when a liquid film is added?

/*!
 * \details Particle read function. Has an std::istream as argument, from which 
 *          it extracts the radius_, invMass_ and invInertia_, respectively. 
 *          From these the mass_ and inertia_ are deduced. An additional set of 
 *          properties is read through the call to the parent's method
 *          BaseParticle::read().
 * \param[in,out] is    input stream with particle properties.
 */
void ThermalParticle::read(std::istream& is)
{
    BaseParticle::read(is);
    std::string dummy;
    is >> dummy >> temperature_;
}

Mdouble ThermalParticle::getTemperature() const
{
    return temperature_;
}

void ThermalParticle::setTemperature(Mdouble temperature)
{
    temperature_ = temperature;
}

void ThermalParticle::addTemperature(Mdouble temperature)
{
    temperature_ += temperature;
}

void ThermalParticle::actionsAfterTimeStep()
{
    if (timeDependentTemperature_)
    {
        temperature_ = timeDependentTemperature_(getHandler()->getDPMBase()->getTime());
    }
    if (getSpecies()->getTemperatureDependentDensity())
    {
        const Mdouble density = getSpecies()->getTemperatureDependentDensity()(temperature_);
        radius_ = getRadius() * cbrt(getMass() / (getVolume() * density));
    }
}

const std::function<double(double)>& ThermalParticle::getTimeDependentTemperature() const
{
    return timeDependentTemperature_;
}

void ThermalParticle::setTimeDependentTemperature(const std::function<double(double)>& timeDependentTemperature)
{
    timeDependentTemperature_ = timeDependentTemperature;
    temperature_ = timeDependentTemperature(0);
    logger(INFO, "Setting initial temperature to %", temperature_);
}
