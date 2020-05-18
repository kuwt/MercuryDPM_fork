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

#ifndef HertzianSinterNormalSpecies_H
#define HertzianSinterNormalSpecies_H

#include "Species/NormalForceSpecies/BaseNormalForce.h"
#include "Math/ExtendedMath.h"
#include "Interactions/NormalForceInteractions/HertzianSinterInteraction.h"

/*!
 * \brief HertzianSinterNormalSpecies contains the parameters used to describe a plastic-cohesive normal force (Stefan Ludings plastic-cohesive force model).
 * \details See SinterNormalInteraction::computeForce for a description of the force law.
 */
class HertzianSinterNormalSpecies : public BaseNormalForce
{
public:
    ///\brief The correct Interaction type for this FrictionForceSpecies
    typedef HertzianSinterInteraction InteractionType;
    
    ///\brief The default constructor.
    HertzianSinterNormalSpecies();
    
    ///\brief The default copy constructor.
    HertzianSinterNormalSpecies(const HertzianSinterNormalSpecies& p);
    
    ///\brief The default destructor.
    ~HertzianSinterNormalSpecies();
    
    /// \brief Reads the species properties from an input stream.
    void read(std::istream& is);
    
    /// \brief Writes the species properties to an output stream.
    void write(std::ostream& os) const;
    
    /// \brief Used in Species::getName to obtain a unique name for each Species.
    std::string getBaseName() const;

// Species-specific functions
    
    ///\brief creates default values for mixed species
    void mix(HertzianSinterNormalSpecies* S, HertzianSinterNormalSpecies* T);


//setters and getters
    
    /*!
     * \brief Sets all parameters of the linear plastic-viscoelastic normal force at once.
     */
    void setPlasticParameters(Mdouble loadingModulus, Mdouble unloadingModulusMax, Mdouble cohesionModulus,
                              Mdouble penetrationDepthMax);
    
    /*!
     * \brief Returns the loading stiffness of the linear plastic-viscoelastic normal force.
     */
    Mdouble getLoadingModulus() const;
    
    /*!
     * \brief Returns the maximum unloading stiffness of the linear plastic-viscoelastic normal force.
     */
    Mdouble getUnloadingModulusMax() const;
    
    /*!
     * \brief Returns the cohesive stiffness of the linear plastic-viscoelastic normal force.
     */
    Mdouble getCohesionModulus() const;
    
    /*!
     * \brief Returns the maximum penetration depth of the linear plastic-viscoelastic normal force.
     */
    Mdouble getPenetrationDepthMax() const;
    
    /*!
     * \brief Sets the loading stiffness of the linear plastic-viscoelastic normal force.
     */
    void setLoadingModulus(Mdouble loadingModulus);
    
    /*!
     * \brief Sets the maximum unloading stiffness of the linear plastic-viscoelastic normal force.
     */
    void setUnloadingModulusMax(Mdouble unloadingModulusMax);
    
    /*!
     * \brief Sets the cohesive stiffness of the linear plastic-viscoelastic normal force.
     */
    void setCohesionModulus(Mdouble cohesionModulus);
    
    /*!
     * \brief Sets the maximum penetration depth of the linear plastic-viscoelastic normal force.
     */
    void setPenetrationDepthMax(Mdouble penetrationDepthMax);
    
    /*!
     * \brief Sets the linear dissipation coefficient of the linear plastic-viscoelastic normal force.
     */
    void setDissipation(Mdouble dissipation);
    
    /*!
     * \brief Allows the normal dissipation to be accessed.
     */
    Mdouble getDissipation() const;
    
    /*!
     * \brief Sets the linear dissipation coefficient of the linear plastic-viscoelastic normal force.
     */
    void setSinterRate(Mdouble sinterRate);
    
    /*!
     * \brief Allows the normal dissipation to be accessed.
     */
    Mdouble getSinterRate() const;
    
    /*!
     * \brief Returns the optimal time step to resolve a collision of two particles of a given mass.
     */
    Mdouble computeTimeStep(Mdouble mass);


private:
    ///(normal) spring constant (k_1)
    Mdouble loadingModulus_;
    
    ///the maximum elastic constant (k_2^max) for plastic deformations
    Mdouble unloadingModulusMax_;
    
    ///the adhesive spring constant (k^c) for plastic deformations
    Mdouble cohesionModulus_;
    
    ///the depth (relative to the normalized radius) at which k_2^max is used (phi_f)
    Mdouble penetrationDepthMax_;
    
    ///linear dissipation coefficient
    Mdouble dissipation_;
    
    Mdouble sinterRate_;
};

#endif
