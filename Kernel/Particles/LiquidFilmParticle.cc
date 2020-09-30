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
#include "Particles/LiquidFilmParticle.h"
#include "Interactions/BaseInteraction.h"
#include "Species/ParticleSpecies.h"
#include "ParticleHandler.h"
#include "DPMBase.h"

/*!
 * \details default constructor, creates an Particle at (0,0,0) with radius, 
 * mass and inertia equal to 1
 */
LiquidFilmParticle::LiquidFilmParticle()
{
    liquidVolume_ = 0;
}

/*!
 * \details Constructor that copies most of the properties of the given particle.
 *          Please note that not everything is copied, for example the position 
 *          in the HGrid is not determined yet by the end of this constructor. 
 *          It also does not copy the interactions and the pointer to the handler
 *          that handles this particle. Use with care.
 * \param[in,out] p  Reference to the LiquidFilmParticle this one should become a copy of.
 */
LiquidFilmParticle::LiquidFilmParticle(const LiquidFilmParticle& p)
        : BaseParticle(p)
{
    liquidVolume_ = p.liquidVolume_;
}

/*!
 * \details Destructor. It asks the ParticleHandler to check if this was the 
 *          smallest or largest particle and adjust itself accordingly.
 */
LiquidFilmParticle::~LiquidFilmParticle()
= default;

/*!
 * \details Copy method. Uses copy constructor to create a copy on the heap. 
 *          Useful for polymorphism.
 * \return pointer to the particle's copy
 */
LiquidFilmParticle* LiquidFilmParticle::copy() const
{
    return new LiquidFilmParticle(*this);
}

/*!
 * \details LiquidFilmParticle print method, which accepts an os std::ostream as 
 *          input. It prints human readable LiquidFilmParticle information to the 
 *          std::ostream.
 * \param[in,out] os    stream to which the info is written
 */
void LiquidFilmParticle::write(std::ostream& os) const
{
    BaseParticle::write(os);
    os << " liquidVolume " << liquidVolume_;
}

/*!
 * \details Returns the name of the object; in this case 'LiquidFilmParticle'.
 * \return The object name.
 */
std::string LiquidFilmParticle::getName() const
{
    return "LiquidFilmParticle";
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
void LiquidFilmParticle::read(std::istream& is)
{
    BaseParticle::read(is);
    std::string dummy;
    is >> dummy >> liquidVolume_;
    // a fix to allow reading of restart files pre-nonspherical
    if (dummy == "invInertia")
    {
        is >> dummy >> liquidVolume_;
    }
}

Mdouble LiquidFilmParticle::getLiquidVolume() const
{
    return liquidVolume_;
}

void LiquidFilmParticle::setLiquidVolume(Mdouble liquidVolume)
{
    liquidVolume_ = liquidVolume;
}

void LiquidFilmParticle::addLiquidVolume(Mdouble liquidVolume)
{
    liquidVolume_ += liquidVolume;
}


unsigned LiquidFilmParticle::getNumberOfFieldsVTK() const
{
    return 3;
}

std::string LiquidFilmParticle::getTypeVTK(unsigned i) const
{
    return "Float32";
}

std::string LiquidFilmParticle::getNameVTK(unsigned i) const
{
    if (i==1)
        return "liquidFilmVolume";
    else if (i==2)
        return "liquidBridgeVolume";
    else /*i=0*/
        return "fullLiquidVolume";
}

std::vector<Mdouble> LiquidFilmParticle::getFieldVTK(unsigned i) const
{
    if (i==1) {
        return std::vector<Mdouble>(1, liquidVolume_);
    } else /*i=2 or 0*/ {
        Mdouble fullLiquidVolume = (i==2)?0:liquidVolume_;
        for (auto k : getInteractions()) {
            auto j = dynamic_cast<LiquidMigrationWilletInteraction*>(k);
            if (j && j->getLiquidBridgeVolume()) {
                fullLiquidVolume += 0.5*j->getLiquidBridgeVolume();
//            } else {
//                logger(WARN,"All contacts of % need to be LiquidMigrationWilletInteraction",i);
            }
        }
        return std::vector<Mdouble>(1, fullLiquidVolume);
    }
}
