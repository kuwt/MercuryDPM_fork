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

#ifndef SINTERINTERACTION_H
#define SINTERINTERACTION_H

#include "Interactions/BaseInteraction.h"

class SinterNormalSpecies;

class BaseInteractable;

/*!
 * \class SinterInteraction
 * \brief Computes normal forces in case of a linear plastic visco-elastic interaction.
 */
class SinterInteraction : public virtual BaseInteraction
{
public:
    /*!
     * \brief An alias for the corresponding species
     */
    typedef SinterNormalSpecies SpeciesType;
    
    /*!
     * \brief Constructor.
     */
    SinterInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp);
    
    /*!
     * \brief Copy constructor.
     */
    SinterInteraction(const SinterInteraction& p);
    
    /// \todo: MX add to MPI
    SinterInteraction();
    
    /*!
     * \brief Destructor.
     */
    ~SinterInteraction() override;
    /*!
     * \brief Creates a copy of an object of this class. (Deep copy)
     */
    //BaseInteraction* copy() const;
    /*!
     * \brief Computes the normal forces due to linear plastic visco elastic interaction.
     */
    void computeNormalForce();
    
    /*!
     * \brief Interaction read function, which accepts an std::istream as input.
     */
    void read(std::istream& is) override;
    
    /*!
     * \brief Interaction write function, which accepts an std::ostream as input.
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief Returns the name of the interaction.
     */
    virtual std::string getBaseName() const;
    
    /*!
     * \brief Computes and returns the amount of elastic energy stored in the spring.
     */
    Mdouble getElasticEnergy() const override;
    
    /*!
     * \brief
     */
    const SinterNormalSpecies* getSpecies() const;
    
    /*!
     * \brief
     */
    Mdouble getPlasticOverlap() const;
    
    /*!
     * \brief
     */
    void setPlasticOverlap(Mdouble plasticOverlap);
    
    /*!
     * \brief
     */
    Mdouble getUnloadingStiffness() const;
    
    unsigned getNumberOfFieldsVTK() const override;
    
    std::string getTypeVTK(unsigned i) const override;
    
    std::string getNameVTK(unsigned i) const override;
    
    std::vector<Mdouble> getFieldVTK(unsigned i) const override;

private:
    
    //set in integrate, used in compute force
    Mdouble plasticOverlap_;
};

#endif
