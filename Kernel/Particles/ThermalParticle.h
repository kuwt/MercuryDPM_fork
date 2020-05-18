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

#ifndef ThermalParticle_H
#define ThermalParticle_H

#include "BaseParticle.h"

/*!
 * \class ThermalParticle
 * \brief
 */
class ThermalParticle final : public BaseParticle
{
public:
    /*!
     * \brief Basic Particle constructor, creates a particle at (0,0,0) with radius, mass and inertia equal to 1
     */
    ThermalParticle();
    
    /*!
     * \brief Particle copy constructor, which accepts as input a reference to a Particle. It creates a copy of this Particle and all it's information. Usually it is better to use the copy() function for polymorfism.
     */
    ThermalParticle(const ThermalParticle& p);
    
    /*!
     * \brief Particle destructor, needs to be implemented and checked if it removes tangential spring information
     */
    ~ThermalParticle() override;
    
    /*!
     * \brief Particle copy method. It calls to copy constructor of this Particle, useful for polymorfism
     */
    ThermalParticle* copy() const override;
    
    
    void write(std::ostream& os) const override;
    
    std::string getName() const override;
    
    void read(std::istream& is) override;
    
    Mdouble getTemperature() const;
    
    void setTemperature(Mdouble temperature);
    
    void addTemperature(Mdouble temperature);
    
    void setTemperatureDependentDensity(const std::function<double(double)>& temperatureDependentDensity);
    
    const std::function<double(double)>& getTemperatureDependentDensity() const;
    
    const std::function<double(double)>& getTimeDependentTemperature() const;
    
    void setTimeDependentTemperature(const std::function<double(double)>& timeDependentTemperature);
    
    void actionsAfterTimeStep() override;

    bool isSphericalParticle() const override {return true;}

private:
    
    /*!
     * Change this function to let the temperature be time-dependent.
     */
    std::function<double(double temperature)> timeDependentTemperature_;
    
    Mdouble temperature_;
};

#endif
