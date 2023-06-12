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

#ifndef ThermalParticle_H
#define ThermalParticle_H
#include "DPMBase.h"

/*!
 * \class Thermal
 * \brief
 */
template<class Particle>
class Thermal : public Particle
{
public:
    /*!
     * \brief Basic Particle constructor, creates a particle at (0,0,0) with radius, mass and inertia equal to 1.
     * \details default constructor
     */
    Thermal()
    {
        temperature_ = 0;
        //temperatureDependentDensity_
    }
    
    /*!
     * \brief Particle copy constructor, which accepts as input a reference to a Particle. It creates a copy of this Particle and all it's information. Usually it is better to use the copy() function for polymorphism.
     * \details Constructor that copies most of the properties of the given particle.
     *          Please note that not everything is copied, for example the position
     *          in the HGrid is not determined yet by the end of this constructor.
     *          It also does not copy the interactions and the pointer to the handler
     *          that handles this particle. Use with care.
     * \param[in,out] p  Reference to the Thermal this one should become a copy of.
     */
    Thermal(const Thermal& p) : Particle(p)
    {
        temperature_ = p.temperature_;
        timeDependentTemperature_ = p.timeDependentTemperature_;
    }
    
    /*!
     * \brief Particle destructor, needs to be implemented and checked if it removes tangential spring information
     * \details Destructor. It asks the ParticleHandler to check if this was the
     *          smallest or largest particle and adjust itself accordingly.
     */
    ~Thermal() override
    = default;
    
    /*!
     * \brief Particle copy method. It calls to copy constructor of this Particle, useful for polymorfism
     * \details Copy method. Uses copy constructor to create a copy on the heap.
     *          Useful for polymorphism.
     * \return pointer to the particle's copy
     */
    Thermal* copy() const override
    {
        return new Thermal(*this);
    }

    /*!
     * \details Thermal print method, which accepts an os std::ostream as
     *          input. It prints human readable Thermal information to the
     *          std::ostream.
     * \param[in,out] os    stream to which the info is written
     */
    void write(std::ostream& os) const override
    {
        Particle::write(os);
        os << " temperature " << temperature_;
    }
    
    std::string getName() const override
    {
        return "Thermal" + Particle::getName();
    }
    
    void read(std::istream& is) override;
    
    Mdouble getTemperature() const
    {
        return temperature_;
    }
    
    void setTemperature(Mdouble temperature)
    {
        logger.assert_debug(temperature>0, "Temperature has to be positive");
        temperature_ = temperature;
    }
    
    void addTemperature(Mdouble temperature)
    {
        temperature_ += temperature;
        logger.assert_debug(temperature_>0, "Temperature has to be positive");
    }
    
    void setTemperatureDependentDensity(const std::function<double(double)>& temperatureDependentDensity);
    
    const std::function<double(double)>& getTemperatureDependentDensity() const;
    
    const std::function<double(double)>& getTimeDependentTemperature() const
    {
        return timeDependentTemperature_;
    }
    
    void setTimeDependentTemperature(const std::function<double(double)>& timeDependentTemperature);
    
    void actionsAfterTimeStep() override;

    bool isSphericalParticle() const override {return true;}

private:
    
    /*!
     * Change this function to let the temperature be time-dependent.
     */
    std::function<double(double temperature)> timeDependentTemperature_;

protected:
    Mdouble temperature_;
};

/*!
 * \details Particle read function. Has an std::istream as argument, from which
 *          it extracts the radius_, invMass_ and invInertia_, respectively.
 *          From these the mass_ and inertia_ are deduced. An additional set of
 *          properties is read through the call to the parent's method
 *          BaseParticle::read().
 * \param[in,out] is    input stream with particle properties.
 */
template<class Particle>
void Thermal<Particle>::read(std::istream& is)
{
    Particle::read(is);
    std::string dummy;
    is >> dummy >> temperature_;
}

template<class Particle>
void Thermal<Particle>::actionsAfterTimeStep()
{
    if (timeDependentTemperature_)
    {
        temperature_ = timeDependentTemperature_(this->getHandler()->getDPMBase()->getTime());
    }
    if (this->getSpecies()->getTemperatureDependentDensity())
    {
        const Mdouble density = this->getSpecies()->getTemperatureDependentDensity()(temperature_);
        this->radius_ = this->getRadius() * cbrt(this->getMass() / (this->getVolume() * density));
    }
}

template<class Particle>
void Thermal<Particle>::setTimeDependentTemperature(const std::function<double(double)>& timeDependentTemperature)
{
    timeDependentTemperature_ = timeDependentTemperature;
    temperature_ = timeDependentTemperature(0);
    logger(INFO, "Setting initial temperature to %", temperature_);
}

typedef Thermal<SphericalParticle> ThermalParticle;
#endif
