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
     * \brief Basic Particle constructor, creates an Particle at (0,0,0) with radius, mass and inertia equal to 1.
     * \details default constructor
     */
    LiquidFilmParticle()
    {
        liquidVolume_ = 0;
    }
    
    /*!
     * \brief Particle copy constructor, which accepts as input a reference to a Particle. It creates a copy of this Particle and all it's information. Usually it is better to use the copy() function for polymorfism.
     * \details Constructor that copies most of the properties of the given particle.
     *          Please note that not everything is copied, for example the position
     *          in the HGrid is not determined yet by the end of this constructor.
     *          It also does not copy the interactions and the pointer to the handler
     *          that handles this particle. Use with care.
     * \param[in,out] p  Reference to the LiquidFilmParticle this one should become a copy of.
     */
    LiquidFilmParticle(const LiquidFilmParticle& p) : BaseParticle(p)
    {
        liquidVolume_ = p.liquidVolume_;
    }
    
    /*!
     * \brief Particle destructor, needs to be implemented and checked if it removes tangential spring information
     * \details Destructor. It asks the ParticleHandler to check if this was the
     *          smallest or largest particle and adjust itself accordingly.
     */
    ~LiquidFilmParticle() override
    = default;
    
    /*!
     * \brief Particle copy method. It calls to copy constructor of this Particle, useful for polymorfism
     * \details Copy method. Uses copy constructor to create a copy on the heap.
     *          Useful for polymorphism.
     * \return pointer to the particle's copy
     */
    LiquidFilmParticle* copy() const override
    {
        return new LiquidFilmParticle(*this);
    }

    /*!
     * \details LiquidFilmParticle print method, which accepts an os std::ostream as
     *          input. It prints human readable LiquidFilmParticle information to the
     *          std::ostream.
     * \param[in,out] os    stream to which the info is written
     */
    void write(std::ostream& os) const override
    {
        BaseParticle::write(os);
        os << " liquidVolume " << liquidVolume_;
    }

    /*!
     * \details Returns the name of the object; in this case 'LiquidFilmParticle'.
     * \return The object name.
     */
    std::string getName() const override
    {
        return "LiquidFilmParticle";
    }
    
    void read(std::istream& is) override;
    
    Mdouble getLiquidVolume() const
    {
        return liquidVolume_;
    }
    
    void setLiquidVolume(Mdouble liquidVolume)
    {
        liquidVolume_ = liquidVolume;
    }
    
    void addLiquidVolume(Mdouble liquidVolume)
    {
        liquidVolume_ += liquidVolume;
    }
    
    unsigned getNumberOfFieldsVTK() const override
    {
        return 3;
    }
    
    std::string getTypeVTK(unsigned i) const override
    {
        return "Float32";
    }
    
    std::string getNameVTK(unsigned i) const override;
    
    std::vector<Mdouble> getFieldVTK(unsigned i) const override;

    bool isSphericalParticle() const override {return true;}

private:
    
    Mdouble liquidVolume_;
};

#endif
