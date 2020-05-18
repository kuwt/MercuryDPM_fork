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

#ifndef LiquidFilmParticle_H
#define LiquidFilmParticle_H

#include "BaseParticle.h"

/*!
 * \class LiquidFilmParticle
 * \brief
 */
class LiquidFilmParticle final : public BaseParticle
{
public:
    /*!
     * \brief Basic Particle constructor, creates an Particle at (0,0,0) with radius, mass and inertia equal to 1
     */
    LiquidFilmParticle();
    
    /*!
     * \brief Particle copy constructor, which accepts as input a reference to a Particle. It creates a copy of this Particle and all it's information. Usually it is better to use the copy() function for polymorfism.
     */
    LiquidFilmParticle(const LiquidFilmParticle& p);
    
    /*!
     * \brief Particle destructor, needs to be implemented and checked if it removes tangential spring information
     */
    ~LiquidFilmParticle() override;
    
    /*!
     * \brief Particle copy method. It calls to copy constructor of this Particle, useful for polymorfism
     */
    LiquidFilmParticle* copy() const override;
    
    
    void write(std::ostream& os) const override;
    
    std::string getName() const override;
    
    void read(std::istream& is) override;
    
    Mdouble getLiquidVolume() const;
    
    void setLiquidVolume(Mdouble liquidVolume);
    
    void addLiquidVolume(Mdouble liquidVolume);
    
    unsigned getNumberOfFieldsVTK() const override;
    
    std::string getTypeVTK(unsigned i) const override;
    
    std::string getNameVTK(unsigned i) const override;
    
    std::vector<Mdouble> getFieldVTK(unsigned i) const override;

    bool isSphericalParticle() const override {return true;}

private:
    
    Mdouble liquidVolume_;
};

#endif
