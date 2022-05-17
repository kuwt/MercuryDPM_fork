//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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

#ifndef SinterLinNormalSpecies_H
#define SinterLinNormalSpecies_H

#include "Species/BaseSpecies.h"
#include "Species/NormalForceSpecies/BaseNormalForce.h"
#include "Math/ExtendedMath.h"
#include "Interactions/NormalForceInteractions/SinterLinInteraction.h"

/*!
 * There are two sintering approaches in this class. FRENKEL = 0, which represents the sintering driven by surface tension,
 * VISCOELASTIC_CONTACT = 1, sintering driven by intersuface forces and surface tension.
*/

enum class SINTER_APPROACH : unsigned char
{
    FRENKEL = 0,
    VISCOELASTIC_CONTACT = 1,
};

/*!
 * \brief SinterLinNormalSpecies contains the parameters used to describe a plastic-cohesive normal force (Stefan Ludings plastic-cohesive force model) based on three different sintering mechanisms.
 * \details Ref: Visco-elastic sintering kinetics in virgin and aged polymer powders. Powder Technology, 2020.
 */

class SinterLinNormalSpecies : public BaseNormalForce
{
public:
    ///\brief The correct Interaction type for this FrictionForceSpecies
    typedef SinterLinInteraction InteractionType;
    
    ///\brief The default constructor.
    SinterLinNormalSpecies();
    
    ///\brief The default copy constructor.
    SinterLinNormalSpecies(const SinterLinNormalSpecies& p);
    
    ///\brief The default destructor.
    ~SinterLinNormalSpecies();
    
    /// \brief Reads the species properties from an input stream.
    void read(std::istream& is);
    
    /// \brief Writes the species properties to an output stream.
    void write(std::ostream& os) const;
    
    /// \brief Used in Species::getName to obtain a unique name for each Species.
    std::string getBaseName() const;

// Species-specific functions
    
    ///\brief creates default values for mixed species
    void mix(SinterLinNormalSpecies* S, SinterLinNormalSpecies* T);
    
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

    MERCURY_DEPRECATED void setLoadingStiffnessAndDissipation(helpers::KAndDisp new_);
    
    /*!
     * \brief Returns the optimal time step to resolve a collision of two particles of a given mass.
     */
    Mdouble computeTimeStep(Mdouble mass);

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
    * \brief Sets sinterRate_
    */
    void setSinterRate(Mdouble sinterRate);

    /*!
    * \brief Accesses sinterRate_.
    */
    Mdouble getSinterRate() const;

    /*!
    * \brief sets the instantaneous compliance (compliance zero)
    */

    void setComplianceZero(Mdouble complianceZero);

    /*!
    * \brief Accesses the instantaneous compliance (compliance zero)
    */

    Mdouble getComplianceZero() const;

    /*!
    * \brief sets the surface tension.
    */

    void setSurfTension(Mdouble complianceZero);

    /*!
    * \brief accesses the surface tension.
    */

    Mdouble getSurfTension() const;

    /*!
    * \brief sets the fluidity (inverse of viscosity).
    */

    void setFluidity(Mdouble complianceOne);

    /*!
    * \brief accesses the fluidity (inverse of viscosity).
    */

    Mdouble getFluidity() const;

    /*!
    * \brief sets the critical separation distance. This mimics inter-surface forces.
    */

    void setSeparationDis(Mdouble separationDis);

    /*!
    * \brief accesses the critical separation distance. This mimics inter-surface forces.
    */

    Mdouble getSeparationDis() const;

    /*!
    * \brief sets the sintering approach selected.
    */

    void setSinterType(SINTER_APPROACH sinterType);

    // Sinter declarations:
    SINTER_APPROACH getSinterType() const;


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
    * Determines sinter rate: d(delta0)/dt = (fa+fep)/nu
    */
    Mdouble sinterRate_;

    /*!
    * Material property. C0 = (1-\nu^2)/E, E: young's Modulus, nu poisson ratio.
    */

    Mdouble complianceZero_;

    /*!
    * Material property.
    */
    Mdouble surfTension_;

    /*!
    * Material property, which represents the inverse of viscosity
    */
    Mdouble fluidity;

    /*!
    * Critical separation distance, to describe inter-surface forces.
    */
    Mdouble separationDis_;

    ///Type of sintering
    SINTER_APPROACH sinterType_;

};

#endif
