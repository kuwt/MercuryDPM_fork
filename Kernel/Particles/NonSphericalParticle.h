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

#ifndef NonSphericalParticle_H
#define NonSphericalParticle_H


#include "BaseParticle.h"

/*!
 * \class NonSphericalParticle
 * \brief Base class for all non-spherical particle types.
 */
class NonSphericalParticle : public BaseParticle
{
public:
    /*!
     * Default constructor.
     */
    NonSphericalParticle() = default;
    
    /*!
     * Default copy constructor.
     */
    NonSphericalParticle(const NonSphericalParticle& p) = default;
    
    /*!
     * \brief Base class copy constructor.
     * Creates a NonSphericalParticle particle from a BaseParticle.
     */
    NonSphericalParticle(const BaseParticle& p) : BaseParticle(p) {}
    
    /**
     * Default destructor.
     */
    ~NonSphericalParticle() override = default;

    /*!
     * Particle copy method. Pure virtual as this needs to set in the derived class.
     */
    NonSphericalParticle* copy() const override = 0;

    /*!
     * Returns the name of the object. Pure virtual as this needs to set in the derived class.
     */
    std::string getName() const override = 0;

    /*!
     * This property is the main characteristic of non-spherical particles. Set here for all derived classes.
     */
    bool isSphericalParticle() const override {
        return false;
    }

    /*!
    * The following redefines functions of BaseParticles as virtual to make them available in child Clump class
    */

    virtual Mdouble getKineticEnergy() const override{
        return BaseParticle::getKineticEnergy();
    }

    virtual Mdouble getRotationalEnergy() const override{
        return BaseParticle::getRotationalEnergy();
    }

};

#endif
