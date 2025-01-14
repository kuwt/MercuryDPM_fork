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

#ifndef LINEARPLASTICVISCOELASTICINTERACTIONCUTOFF_H
#define LINEARPLASTICVISCOELASTICINTERACTIONCUTOFF_H

#include "Interactions/BaseInteraction.h"
class LinearPlasticViscoelasticNormalSpeciesExtended;
class BaseInteractable;
/*!
* \class LinearPlasticViscoelasticInteractionCutoff
* \brief Computes normal forces in case of a linear plastic visco-elastic interaction.
*/
class LinearPlasticViscoelasticInteractionCutoff : public virtual BaseInteraction
{
public:
   /*!
   * \brief An alias for the species (needed for e.g. ThermalSpecies)
   */
   typedef LinearPlasticViscoelasticNormalSpeciesExtended SpeciesType;
   /*!
   * \brief Constructor.
   */
   LinearPlasticViscoelasticInteractionCutoff(BaseInteractable* P, BaseInteractable* I, Mdouble timeStamp);
   /*!
   * \brief Copy constructor.
   */
   LinearPlasticViscoelasticInteractionCutoff(const LinearPlasticViscoelasticInteractionCutoff &p);
   /*!
   * \brief Destructor.
   */
   virtual ~LinearPlasticViscoelasticInteractionCutoff();
   /*!
   * \brief Creates a copy of an object of this class. (Deep copy)
   */
   //BaseInteraction* copy() const;
   /*!
   * \brief Computes the normal forces due to linear plastic visco elastic interaction.
   */
   void computeLinearPlasticViscoelasticForce();
   /*!
   * \brief Calls computeLinearPlasticViscoElasticForce().
   */
   void computeNormalForce();
   /*!
   * \brief Interaction read function, which accepts an std::istream as input.
   */
   virtual void read(std::istream& is);
   /*!
   * \brief Interaction write function, which accepts an std::ostream as input.
   */
   virtual void write(std::ostream& os) const;
   /*!
   * \brief Returns the name of the interaction.
   */
   virtual std::string getBaseName() const;
   /*!
   * \brief Computes and returns the amount of elastic energy stored in the spring.
   */
   Mdouble getElasticEnergy() const;
   /*!
   * \brief
   */
   const LinearPlasticViscoelasticNormalSpeciesExtended* getSpecies() const;
   /*!
   * \brief
   */
   Mdouble getMaxOverlap() const;
   /*!
   * \brief
   */
   void setMaxOverlap(const Mdouble maxOverlap);
   /*!
   * \brief
   */
   Mdouble getUnloadingStiffness() const;

private:

   //set in integrate, used in compute force
   Mdouble maxOverlap_;

   // bool bound;
};
#endif
