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

#ifndef HertzianBSHPViscoelasticInteraction_H
#define HertzianBSHPViscoelasticInteraction_H

#include "Interactions/NormalForceInteractions/HertzianViscoelasticInteraction.h"

class BaseInteractable;

class HertzianBSHPViscoelasticNormalSpecies;

/*!
 * \class HertzianBSHPViscoelasticInteraction
 * \brief Computes normal forces for a Herztian visco-elastic interaction.
 */
class HertzianBSHPViscoelasticInteraction : public virtual HertzianViscoelasticInteraction
{
public:
    /*!
     * \brief An alias for HertzianViscoelasticNormalSpecies.
     */
    typedef HertzianBSHPViscoelasticNormalSpecies SpeciesType;
    
    /*!
     * \brief Constructor.
     */
    HertzianBSHPViscoelasticInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp);
    
    //used for mpi
    HertzianBSHPViscoelasticInteraction();
    
    /*!
     * \brief Copy constructor.
     */
    HertzianBSHPViscoelasticInteraction(const HertzianBSHPViscoelasticInteraction& p);
    
    /*!
     * \brief Destructor.
     */
    ~HertzianBSHPViscoelasticInteraction() override;
    //
    //    HertzianBSHPViscoelasticInteraction* copy() const;
    //
    /*!
     * \brief Computes the amount of normal force due to an Hertzian visco-elastic interaction.
     */
    void computeNormalForce();
    
    explicit HertzianBSHPViscoelasticInteraction(const HertzianViscoelasticInteraction& p) : HertzianViscoelasticInteraction(p)
    {
        
    }
    
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
    std::string getBaseName() const;
    
    /*!
     * \brief Returns a const pointer of type HerztianBSHPViscoelasticNormalSpecies (static-cast).
     */
    const HertzianBSHPViscoelasticNormalSpecies* getSpecies() const;
};

#endif
