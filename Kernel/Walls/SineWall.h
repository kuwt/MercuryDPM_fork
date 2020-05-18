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

#ifndef SINEWALL_H
#define SINEWALL_H

#include "BaseWall.h"
#include "Math/Vector.h"

/*! A half-infinite surface which has a sinusoidal direction in the x direction,
 * and no variation in the y direction.
 *
 * The equation of the surface is
 *   z = sw_amp * sin( sw_wavn * x + sw_phshift )  for x < l
 * Because this is nonlinear, we must solve for the interaction location using
 * Newton's method.
 * */
class SineWall : public BaseWall
{
public:
    /*!
     * \brief Default constructor, sets a chute with default parameters.
     */
    SineWall();
    
    /*!
     * \brief Copy constructor
     */
    SineWall(const SineWall& other);
    
    /*!
     * \brief Constructor in which all parameters are set.
     */
    SineWall(Mdouble length, Mdouble sw_wavn, Mdouble sw_phshift, Mdouble sw_amp);
    
    ~SineWall() override;
    
    void set(Mdouble length, Mdouble sw_wavn, Mdouble sw_phshift, Mdouble sw_amp);
    
    SineWall* copy() const override;
    
    bool getDistanceAndNormal(const BaseParticle& P, Mdouble& distance, Vec3D& normal_return) const override;
    
    BaseInteraction*
    getInteractionWith(BaseParticle* p, unsigned timeStamp, InteractionHandler* interactionHandler) override;
    
    void read(std::istream& is) override;
    
    void write(std::ostream& os) const override;
    
    std::string getName() const override;

private:
    Mdouble l_;
    Mdouble sw_amp_; // amplitude of sinusoidal variations
    Mdouble sw_wavn_; // (angular) wavenumber of sinusoidal variations
    Mdouble sw_phshift_; // phase shift of sinusoidal variations
};

#endif
