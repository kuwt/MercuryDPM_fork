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

#ifndef LiquidFilmParticle_H
#define LiquidFilmParticle_H
#include "DPMBase.h"

/*!
 * \class LiquidFilm
 * \brief
 */
template<class Particle>
class LiquidFilm : public Particle
{
public:
    /*!
     * \brief Basic Particle constructor, creates an Particle at (0,0,0) with radius, mass and inertia equal to 1.
     * \details default constructor
     */
    LiquidFilm()
    {
        liquidVolume_ = 0;
        totalEvaporatedLiquidVolume_ = 0;
    }
    
    /*!
     * \brief Particle copy constructor, which accepts as input a reference to a Particle. It creates a copy of this Particle and all it's information. Usually it is better to use the copy() function for polymorfism.
     * \details Constructor that copies most of the properties of the given particle.
     *          Please note that not everything is copied, for example the position
     *          in the HGrid is not determined yet by the end of this constructor.
     *          It also does not copy the interactions and the pointer to the handler
     *          that handles this particle. Use with care.
     * \param[in,out] p  Reference to the LiquidFilm this one should become a copy of.
     */
    LiquidFilm(const LiquidFilm& p) : Particle(p)
    {
        liquidVolume_ = p.liquidVolume_;
        totalEvaporatedLiquidVolume_ = p.totalEvaporatedLiquidVolume_;
    }
    
    /*!
     * \brief Particle destructor, needs to be implemented and checked if it removes tangential spring information
     * \details Destructor. It asks the ParticleHandler to check if this was the
     *          smallest or largest particle and adjust itself accordingly.
     */
    ~LiquidFilm() override
    = default;
    
    /*!
     * \brief Particle copy method. It calls to copy constructor of this Particle, useful for polymorfism
     * \details Copy method. Uses copy constructor to create a copy on the heap.
     *          Useful for polymorphism.
     * \return pointer to the particle's copy
     */
    LiquidFilm* copy() const override
    {
        return new LiquidFilm(*this);
    }

    /*!
     * \details LiquidFilm print method, which accepts an os std::ostream as
     *          input. It prints human readable LiquidFilm information to the
     *          std::ostream.
     * \param[in,out] os    stream to which the info is written
     */
    void write(std::ostream& os) const override
    {
        Particle::write(os);
        os << " liquidVolume " << liquidVolume_ << " totalEvaporatedLiquidVolume " << totalEvaporatedLiquidVolume_;
    }

    /*!
     * \details Returns the name of the object; in this case 'LiquidFilm'.
     * \return The object name.
     */
    std::string getName() const override
    {
        return "LiquidFilm" + Particle::getName();
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

    Mdouble getFullLiquidVolume() const
    {
        return liquidVolume_ + getLiquidBridgeVolume();
    }

    Mdouble getLiquidBridgeVolume() const;

    Mdouble getTotalEvaporatedLiquidVolume() const
    {
        return totalEvaporatedLiquidVolume_;
    }

    void setTotalEvaporatedLiquidVolume(Mdouble liquidVolume)
    {
        totalEvaporatedLiquidVolume_ = liquidVolume;
    }

    void addTotalEvaporatedLiquidVolume(Mdouble liquidVolume)
    {
        totalEvaporatedLiquidVolume_ += liquidVolume;
    }
    
    unsigned getNumberOfFieldsVTK() const override
    {
        return 4;
    }
    
    std::string getTypeVTK(unsigned i) const override
    {
        return "Float32";
    }
    
    std::string getNameVTK(unsigned i) const override;
    
    std::vector<Mdouble> getFieldVTK(unsigned i) const override;

    bool isSphericalParticle() const override {return true;}

protected:
    
    Mdouble liquidVolume_, totalEvaporatedLiquidVolume_;
};


//todo Does mass and interaction radius change when a liquid film is added?
/*!
 * \details Particle read function. Has an std::istream as argument, from which
 *          it extracts the radius_, invMass_ and invInertia_, respectively.
 *          From these the mass_ and inertia_ are deduced. An additional set of
 *          properties is read through the call to the parent's method
 *          BaseParticle::read().
 * \param[in,out] is    input stream with particle properties.
 */
template<class Particle>
void LiquidFilm<Particle>::read(std::istream& is)
{
    Particle::read(is);
    std::string dummy;
    is >> dummy >> liquidVolume_;
    // a fix to allow reading of restart files pre-nonspherical
    if (dummy == "invInertia")
    {
        is >> dummy >> liquidVolume_;
    }
    is >> dummy >> totalEvaporatedLiquidVolume_;
}

template<class Particle>
std::string LiquidFilm<Particle>::getNameVTK(unsigned i) const
{
    if (i==1)
        return "liquidFilmVolume";
    else if (i==2)
        return "liquidBridgeVolume";
    else if (i==3)
        return "totalEvaporatedLiquidVolume";
    else /*i=0*/
        return "fullLiquidVolume";
}

template<class Particle>
std::vector<Mdouble> LiquidFilm<Particle>::getFieldVTK(unsigned i) const
{
    if (i==1)
        return { liquidVolume_ };
    else if (i==2)
        return { getLiquidBridgeVolume() };
    else if (i==3)
        return { totalEvaporatedLiquidVolume_ };
    else /*i=0*/
        return { getFullLiquidVolume() };
}

template<class Particle>
Mdouble LiquidFilm<Particle>::getLiquidBridgeVolume() const
{
    // Sums the volume of all liquid bridges. Half when interacting with a particle, full when interacting with a wall.
    Mdouble volume = 0.0;
    for (auto i : this->getInteractions())
    {
        auto j = dynamic_cast<LiquidMigrationWilletInteraction*>(i);
        if (j)
        {
            if (dynamic_cast<BaseParticle*>(j->getI()))
                volume += 0.5 * j->getLiquidBridgeVolume();
            else
                volume += j->getLiquidBridgeVolume();
        }
        auto k = dynamic_cast<LiquidMigrationLSInteraction*>(i);
        if (k)
        {
            if (dynamic_cast<BaseParticle*>(k->getI()))
                volume += 0.5 * k->getLiquidBridgeVolume();
            else
                volume += k->getLiquidBridgeVolume();
        }
    }
    return volume;
}

typedef LiquidFilm<SphericalParticle> LiquidFilmParticle;

#endif