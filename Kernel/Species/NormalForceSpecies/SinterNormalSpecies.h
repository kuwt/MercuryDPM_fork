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

#ifndef SinterNormalSpecies_H
#define SinterNormalSpecies_H

#include "Species/NormalForceSpecies/BaseNormalForce.h"
#include "Math/ExtendedMath.h"
#include "Interactions/NormalForceInteractions/SinterInteraction.h"
#include "Math/Helpers.h"

enum class SINTERTYPE : unsigned char
{
    PARHAMI_MCKEEPING = 0,
    CONSTANT_RATE = 1,
    TEMPERATURE_DEPENDENT_FRENKEL = 2,
    REGIME_SINTERING = 3
};

/*!
 * \brief SinterNormalSpecies contains the parameters used to describe a plastic-cohesive normal force (Stefan Ludings plastic-cohesive force model).
 * \details See SinterNormalInteraction::computeForce for a description of the force law.
 */
class SinterNormalSpecies : public BaseNormalForce
{
public:
    ///\brief The correct Interaction type for this FrictionForceSpecies
    typedef SinterInteraction InteractionType;
    
    ///\brief The default constructor.
    SinterNormalSpecies();
    
    ///\brief The default copy constructor.
    SinterNormalSpecies(const SinterNormalSpecies& p);
    
    ///\brief The default destructor.
    ~SinterNormalSpecies();
    
    /// \brief Reads the species properties from an input stream.
    void read(std::istream& is);
    
    /// \brief Writes the species properties to an output stream.
    void write(std::ostream& os) const;
    
    /// \brief Used in Species::getName to obtain a unique name for each Species.
    std::string getBaseName() const;

// Species-specific functions
    
    ///\brief creates default values for mixed species
    void mix(SinterNormalSpecies* S, SinterNormalSpecies* T);
    
    ///Set k, disp such that is matches a given tc and eps for a collision of two different masses.
    ///Recall the resitution constant is a function of k, disp and the mass of each particle in the collision
    /// See also setCollisionTimeAndRestitutionCoefficient(Mdouble tc, Mdouble eps, Mdouble mass)
    void setCollisionTimeAndRestitutionCoefficient(Mdouble tc, Mdouble eps, Mdouble mass);
    
    ///Sets k, disp such that it matches a given tc and eps for a collision of two copies of P
    void setStiffnessAndRestitutionCoefficient(Mdouble k_, Mdouble eps, Mdouble mass);
    
    /*!
 * \brief Calculates collision time for two copies of a particle of given disp, k, mass
 */
    Mdouble getCollisionTime(Mdouble mass);

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
     * \brief Allows the normal dissipation to be accessed.
     */
    Mdouble getDissipation() const;
    
    /*!
     * \brief Sets sinterAdhesion_.
     */
    void setSinterAdhesion(Mdouble sinterAdhesion);
    
    /*!
     * \brief Accesses sinterAdhesion_.
     */
    Mdouble getSinterAdhesion() const;
    
    /*!
     * \brief Sets inverseSinterViscosity_.
     */
    void setInverseSinterViscosity(Mdouble inverseSinterViscosity);
    
    /*!
     * \brief Accesses inverseSinterViscosity_.
     */
    Mdouble getInverseSinterViscosity() const;
    
    /*!
     * \brief Allows the spring and dissipation constants to be changed simultaneously.
     */
    MERCURY_DEPRECATED void setLoadingStiffnessAndDissipation(helpers::KAndDisp new_);
    
    /*!
     * \brief Returns the optimal time step to resolve a collision of two particles of a given mass.
     */
    Mdouble computeTimeStep(Mdouble mass);
    
    void setSinterForceAndTime(Mdouble adhesionForce, Mdouble sinterTime, Mdouble radius);
    
    /*!
     * \brief Sets the sinterAdhesion_ and inverseSinterViscosity_ based on the Parhami-McKeeping parameters.
     */
    void setParhamiMcKeeping
            (Mdouble alpha, Mdouble beta, Mdouble atomicVolume /*Omega*/, Mdouble surfaceEnergy /*gamma_s*/,
             Mdouble thicknessDiffusion /*deltaB*D0B*/, Mdouble activationEnergy /*QB*/, Mdouble temperature /*T*/);
    
    /*!
     * \brief Sets sinterRate_
     */
    void setSinterRate(Mdouble sinterRate);
    
    /*!
     * \brief Sets sinterRate_
     */
    void setSinterType(SINTERTYPE sinterType);
    
    /*!
     * \brief Accesses sinterRate_.
     */
    Mdouble getSinterRate() const;
    
    SINTERTYPE getSinterType() const;
    
    std::function<double(double temperature)> getTemperatureDependentSinterRate() const;
    
    double getTemperatureDependentSinterRate(double temperature) const;
    
    void setTemperatureDependentSinterRate(std::function<double(double temperature)> temperatureDependentSinterRate);

    void setComplianceZero(Mdouble complianceZero);

    Mdouble getComplianceZero() const;

    void setSurfTension(Mdouble complianceZero);

    Mdouble getSurfTension() const;

    void setConstantC1(Mdouble constantC1_);

    Mdouble getConstantC1() const;

    void setSeparationDis(Mdouble separationDis);

    Mdouble getSeparationDis() const;


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
    
    /*!
     * Determines sinter adhesion force fa=sinterAdhesion_*radius in sinter rate: d(delta0)/dt = (fa+fep)/nu
     */
    Mdouble sinterAdhesion_;
    
    /*!
     * Determines sinter viscosity nu = contactRadius^4/inverseSinterViscosity_ in sinter rate: d(delta0)/dt = (fa+fep)/nu
     */
    Mdouble inverseSinterViscosity_;
    
    /*!
     * Determines sinter rate: d(delta0)/dt = (fa+fep)/nu
     */
    Mdouble sinterRate_;

    /*!
    * Compliance 0, C_0, corresponds to the inverse of stiffness at the instantaneous time.
    */
    Mdouble complianceZero_;

    /*!
    * Material constant
    */

    Mdouble surfTension_;

    /*!
    * Material constant
    */

    Mdouble constantC1_;

    /*!
    * Separation distance comes from the fracture theory.
    */
    Mdouble separationDis_;
    
    /*!
     * \brief sinterType options determin how the rate of sintering, d(equilibriumOverlap)/dt is computed
     * PARHAMI_MCKEEPING: sinter rate given by (sinterAdhesion+normalForce)/sinterViscosity
     * CONSTANT_RATE: sinter rate given by sinterRate
    */
    SINTERTYPE sinterType_;
    
    /*!
     *
     */
    std::function<double(double temperature)> temperatureDependentSinterRate_ = [this](
            double temperature) { return sinterRate_; };
};

#endif
