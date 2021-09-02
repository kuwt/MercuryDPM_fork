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

#ifndef MINDLINSPECIES_H
#define MINDLINSPECIES_H

#include "Species/FrictionForceSpecies/BaseFrictionForce.h"
#include "Math/ExtendedMath.h"
#include "Interactions/FrictionForceInteractions/MindlinInteraction.h"

/*!
 * \brief MindlinSpecies contains the parameters used to describe sliding friction.
 * \details See MindlinInteraction::computeForce for a description of the force law.
 */
class MindlinSpecies : public BaseFrictionForce
{
public:
    ///\brief The correct Interaction type for this FrictionForceSpecies
    typedef MindlinInteraction InteractionType;
    
    ///\brief The default constructor.
    MindlinSpecies();
    
    ///\brief The default copy constructor.
    MindlinSpecies(const MindlinSpecies& s);
    
    ///\brief The default destructor.
    ~MindlinSpecies();
    
    /// \brief Reads the species properties from an input stream.
    void read(std::istream& is);
    
    /// \brief Writes the species properties to an output stream.
    void write(std::ostream& os) const;
    
    /// \brief Used in Species::getName to obtain a unique name for each Species.
    virtual std::string getBaseName() const;

//setters and getters
    
    ///Allows the tangential viscosity to be changed
    void setSlidingDissipation(Mdouble new_dispt);
    
    ///Allows the tangential viscosity to be accessed
    Mdouble getSlidingDissipation() const;
    
    ///Allows the (dynamic) Coulomb friction coefficient to be changed; also sets mu_s by default
    //mu has to be set to allow tangential forces (sets dispt=disp as default)
    void setSlidingFrictionCoefficient(Mdouble new_mu);
    
    ///Allows the (dynamic) Coulomb friction coefficient to be accessed
    Mdouble getSlidingFrictionCoefficient() const;
    
    ///Allows the static Coulomb friction coefficient to be changed
    void setSlidingFrictionCoefficientStatic(Mdouble new_mu);
    
    ///Allows the static Coulomb friction coefficient to be accessed
    Mdouble getSlidingFrictionCoefficientStatic() const;
    
    ///\brief Allows the poisson ratio to be changed
    MERCURY_DEPRECATED void setPoissonRatio(Mdouble poissonRatio);

    ///\brief allows the shear modulus to be changed
    void setEffectiveShearModulus(Mdouble shearModulus);
    
    ///\brief Allows the shear modulus to be accessed
    Mdouble getEffectiveShearModulus() const;

    /*!
     * \brief Returns true if torques have to be calculated.
     */
    bool getUseAngularDOFs() const override;
    
    ///\brief creates default values for mixed species
    void mix(MindlinSpecies* S, MindlinSpecies* T);

private:
    
    /*! 
     * \brief tangential dissipation coefficient. 
     * \details Typically set to 2/7 of the normal force dissipation, as both
     * the tangential and the normal spring have the same restitution coefficients 
     * (if the tangential and normal stiffnesses also have a ratio of 2/7). 
     * should always satisfy \f$4*dispt*dt<mass, dispt \approx disp\f$, otherwise 
     * the restitution is zero and the particles stick.
     */
    Mdouble slidingDissipation_;
    
    /// (dynamic) Coulomb friction coefficient
    Mdouble slidingFrictionCoefficient_;
    
    /// static Coulomb friction coefficient (by default set equal to mu)
    Mdouble slidingFrictionCoefficientStatic_;
  
    /// tangential spring constant
    Mdouble shearModulus_;
    
    
};

#endif
