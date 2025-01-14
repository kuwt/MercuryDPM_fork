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

#ifndef LinearPlasticViscoelasticNormalSpecies_H
#define LinearPlasticViscoelasticNormalSpecies_H

#include "Species/NormalForceSpecies/BaseNormalForce.h"
#include "Math/ExtendedMath.h"
#include "Interactions/NormalForceInteractions/LinearPlasticViscoelasticInteraction.h"
#include "Species/FrictionForceSpecies/SlidingFrictionSpecies.h"
#include "Math/Helpers.h"


/*!
 * \brief LinearPlasticViscoelasticNormalSpecies contains the parameters used to describe a plastic-cohesive normal force (Stefan Ludings plastic-cohesive force model).
 * \details See LinearPlasticViscoelasticNormalInteraction::computeForce for a description of the force law.
 */
class LinearPlasticViscoelasticNormalSpecies : public BaseNormalForce
{
public:
    ///\brief The correct Interaction type for this FrictionForceSpecies
    typedef LinearPlasticViscoelasticInteraction InteractionType;
    
    ///\brief The default constructor.
    LinearPlasticViscoelasticNormalSpecies();
    
    ///\brief The default copy constructor.
    LinearPlasticViscoelasticNormalSpecies(const LinearPlasticViscoelasticNormalSpecies& p);
    
    ///\brief The default destructor.
    ~LinearPlasticViscoelasticNormalSpecies();
    
    /// \brief Reads the species properties from an input stream.
    void read(std::istream& is);
    
    /// \brief Writes the species properties to an output stream.
    void write(std::ostream& os) const;
    
    /// \brief Used in Species::getName to obtain a unique name for each Species.
    std::string getBaseName() const;

// Species-specific functions
    
    ///\brief creates default values for mixed species
    void mix(LinearPlasticViscoelasticNormalSpecies* S, LinearPlasticViscoelasticNormalSpecies* T);
    
    ///Set k, disp such that is matches a given tc and eps for a collision of two different masses.
    ///Recall the resitution constant is a function of k, disp and the mass of each particle in the collision
    /// See also setCollisionTimeAndRestitutionCoefficient(Mdouble tc, Mdouble eps, Mdouble mass)
    void setCollisionTimeAndRestitutionCoefficient(Mdouble tc, Mdouble eps, Mdouble mass);
    
    ///Sets k, disp such that it matches a given tc and eps for a collision of two copies of P
    void setStiffnessAndRestitutionCoefficient(Mdouble k_, Mdouble eps, Mdouble mass);
    
    ///Sets disp to obtain a restitution coefficient eps for a collision of two particles of mass m
    void setRestitutionCoefficient(double eps, Mdouble mass);
    
    /*!
 * \brief Calculates collision time for two copies of a particle of given disp, k, mass
 */
    Mdouble getCollisionTime(Mdouble mass) const;
    
    ///Calculates restitution coefficient for two copies of given disp, k, mass
    Mdouble getRestitutionCoefficient(Mdouble mass) const;

//setters and getters
    
    /*!
     * \brief Sets all parameters of the linear plastic-viscoelastic normal force at once.
     */
    void setPlasticParameters(Mdouble loadingStiffness, Mdouble unloadingStiffnessMax, Mdouble cohesionStiffness,
                              Mdouble penetrationDepthMax);
    
    /*!
     * \brief Returns the loading stiffness of the linear plastic-viscoelastic normal force.
     */
    Mdouble getLoadingStiffness() const;
    
    /*!
     * \brief Returns the maximum unloading stiffness of the linear plastic-viscoelastic normal force.
     */
    Mdouble getUnloadingStiffnessMax() const;
    
    /*!
     * \brief Returns the cohesive stiffness of the linear plastic-viscoelastic normal force.
     */
    Mdouble getCohesionStiffness() const;
    
    /*!
     * \brief Returns the maximum penetration depth of the linear plastic-viscoelastic normal force.
     */
    Mdouble getPenetrationDepthMax() const;
    
    /*!
     * \brief Sets the loading stiffness of the linear plastic-viscoelastic normal force.
     */
    void setLoadingStiffness(Mdouble loadingStiffness);
    
    /*!
     * \brief Sets the maximum unloading stiffness of the linear plastic-viscoelastic normal force.
     */
    void setUnloadingStiffnessMax(Mdouble unloadingStiffnessMax);
    
    /*!
     * \brief Sets the cohesive stiffness of the linear plastic-viscoelastic normal force.
     */
    void setCohesionStiffness(Mdouble cohesionStiffness);
    
    /*!
     * \brief Sets the maximum penetration depth of the linear plastic-viscoelastic normal force.
     */
    void setPenetrationDepthMax(Mdouble penetrationDepthMax);
    
    /*!
     * \brief Sets the linear dissipation coefficient of the linear plastic-viscoelastic normal force.
     */
    void setDissipation(Mdouble dissipation);
    
    /*!
     * \brief Allows the spring and dissipation constants to be changed simultaneously.
     */
    MERCURYDPM_DEPRECATED void setLoadingStiffnessAndDissipation(helpers::KAndDisp new_);
    
    /*!
     * \brief Returns the optimal time step to resolve a collision of two particles of a given mass.
     */
    Mdouble computeTimeStep(Mdouble mass);
    
    /*!
     * \brief Allows the normal dissipation to be accessed.
     */
    Mdouble getDissipation() const;

    /*!
     * \brief
     * 1) Computes the maximum plastic overlap
     *    delta_p* = phi*r
     * 2) Computes the overlap at which the maximum adhesive force is generated:
     *    delta_c* = delta_p* / (1+k_c/k_2*)
     * 3) Computes the maximum adhesive force
     *    f_c* = k_c * delta_c*
     * 4) Computes the maximum bond number
     *    Bo* = f_c* / (m*g)
     */
    Mdouble computeBondNumberMax(Mdouble harmonicMeanRadius, Mdouble gravitationalAcceleration) const;

    bool getDoConstantUnloadingStiffness() const {return doConstantUnloadingStiffness_;}

    void setDoConstantUnloadingStiffness(bool doConstantUnloadingStiffness) {doConstantUnloadingStiffness_ = doConstantUnloadingStiffness;}

private:
    ///(normal) spring constant (k_1)
    Mdouble loadingStiffness_;
    
    ///the maximum elastic constant (k_2^max) for plastic deformations
    Mdouble unloadingStiffnessMax_;
    
    ///the adhesive spring constant (k^c) for plastic deformations
    Mdouble cohesionStiffness_;
    
    ///the depth (relative to the normalized radius) at which k_2^max is used (phi_f)
    Mdouble penetrationDepthMax_;
    
    ///linear dissipation coefficient
    Mdouble dissipation_;

    //whether unloading stiffness is variable (Luding) or constant (WaltonBraun)
    bool doConstantUnloadingStiffness_ = false;
};

#endif
