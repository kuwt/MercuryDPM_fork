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

#ifndef SPHINTERACTION_H
#define SPHINTERACTION_H

#include "Interactions/BaseInteraction.h"

class BaseInteractable;

class SPHNormalSpecies;

/*!
 * \class SPHInteraction
 * \brief Enables the use of SPH particles within MercuryDPM
 */
class SPHInteraction : public virtual BaseInteraction
{
public:
    /*!
     * \brief An alias for SPHNormalSpecies
     */
    typedef SPHNormalSpecies SpeciesType;
    
    /*!
     * \brief Constructor.
     */
    SPHInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp);
    
    //used for mpi
    SPHInteraction();
    
    /*!
     * \brief Copy constructor.
     */
    SPHInteraction(const SPHInteraction& p);
    
    /*!
     * \brief Destructor.
     */
    ~SPHInteraction() override;
    /*!
     * \brief Creates a copy of an object of this class. (Deep copy)
     */
    //LinearViscoelasticInteraction* copy() const;
    /*!
     * \brief Computes the normal force caused due to normal collision between particles or particle-wall.
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
     * \brief Returns the type/name of interaction (linear visco-elastic interaction).
     */
    std::string getBaseName() const;
    
    /*!
     * \brief Returns the amount of energy stored in the spring due to head on collision.
     */
    Mdouble getElasticEnergy() const override;
    
    /*!
     * \brief Returns a const pointer of type LinearViscoelasticNormalSpecies*
     */
    const SPHNormalSpecies* getSpecies() const;
    
    void copyHistoryInteraction();

    Mdouble getElasticEnergyAtEquilibrium(Mdouble adhesiveForce) const override;
    
};

#endif
