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
#ifndef StandardFields_H
#define StandardFields_H

#include <Math/Matrix.h>
#include <Math/MatrixSymmetric.h>
#include <CG/Functions/IntegralType.h>
#include <array>

class BaseParticle;

class BaseInteraction;

namespace CGFields
{

/*!
 * \brief Contains the computed field values, like density, momentum and stress.
 * \details CGPoints inherits from this class; CGPoints::evaluate adds to the 
 * values of these variables.
 * \todo These are currently the only fields that are computed. 
 * However, this class is destined to be extended to
 * contain additional information such as fabric, energy, local angular 
 * momentum. Also, a simpler version is planned, where only particle statistics 
 * are evaluated (density and momentum).
 */
class StandardFields
{
public:
    
    /*!
     * \brief Default constructor, sets all field values to zero.
     */
    StandardFields();
    
    /*!
     * \brief Default copy constructor, copies the values of all fields.
     */
    StandardFields(const StandardFields& P) = default;
    
    /*!
     * \brief Destructor, it simply destructs the StandardFields and all the objects
     * it contains.
     */
    ~StandardFields() = default;
    
    static void writeNames(std::ostream& os, unsigned countVariables);
    
    /*!
     * \brief Writes class content into an output stream, typically a stat file.
     */
    void write(std::ostream& os) const;
    
    /*!
     * \brief Writes human-readable class content into an output stream, typically a stat file.
     */
    void output(std::ostream& os) const;
    
    /*!
     * \brief Sets all fields to zero.
     */
    void setZero();
    
    /*!
     * \brief Returns the square of all field values (to calculate standard deviation).
     */
    StandardFields getSquared() const;
    
    /*!
     * \brief Copies all field values.
     */
    StandardFields& operator=(const StandardFields& P);
    
    /*!
     * \brief Adds the field values on the RHS to the LHS of the equation.
     */
    StandardFields& operator+=(const StandardFields& P);
    
    /*!
     * \brief Subtracts the field values on the RHS from the LHS of the equation.
     */
    StandardFields& operator-=(const StandardFields& P);
    
    /*!
     * \brief Divides the field values on the LHS by the RHS of the equation.
     */
    StandardFields& operator/=(Mdouble a);
    
    /*!
     * \brief Multiplies the field values on the left of the '*' by the
     * scalar value on the right of the '*' and returns the answer.
     */
    StandardFields operator*(Mdouble a) const;
    
    /*!
     * \brief This function should be called from within a loop over all
     * particles to compute all the fields that are defined as a sum over all
     * particles (e.g. density, momentum).
     */
    void addParticleStatistics(Mdouble phi, const StandardFields& currentInteraction);
    
    void addParticleDifferentialStatistics(Vec3D& dphi, const StandardFields& currentInteraction);
    
    /*!
     * \brief This function should be called from within a loop over all
     * Interactions to compute all the fields that are defined as a sum over all
     * Interactions (e.g. stress).
     */
    void addInteractionStatistics(Mdouble psi, const StandardFields& currentInteraction);
    
    /*!
     * \brief This function should be called from within a loop over all
     * Interactions to compute all the fields that are defined as a sum over all
     * Interactions with external objects (e.g. IFD).
     */
    void addContactPointStatistics(Mdouble phi, const StandardFields& currentInteraction);
    
    void setFields(const BaseInteraction& c, IntegralType type);
    
    void setCylindricalFields(const BaseInteraction& c, IntegralType type);
    
    void setFields(const BaseParticle& p);
    
    void setCylindricalFields(const BaseParticle& p);
    
    
    /*!
     * \brief Returns true if the class contains fields that are defined as a
     * sum over all Interactions (e.g. stress), else returns false.
     */
    static bool doInteractionStatistics();
    
    Mdouble getVolumeFraction() const
    {
        return volumeFraction_;
    }
    
    Mdouble getDensity() const
    {
        return density_;
    }
    
    Vec3D getMomentum() const
    {
        return momentum_;
    }
    
    MatrixSymmetric3D getMomentumFlux() const
    {
        return momentumFlux_;
    }
    
    Matrix3D getContactStress() const
    {
        return contactStress_;
    }
    
    Vec3D getInteractionForceDensity() const
    {
        return interactionForceDensity_;
    }
    
    Mdouble getParticleSizeDensity(size_t i) const
    {
        return particleSizeDensity_[i];
    }
    
    std::array<Mdouble, 6> getParticleSizeDensity() const
    {
        return particleSizeDensity_;
    }
    
    std::array<Mdouble, 6> getParticleSizeMomenta() const;
    
    std::array<Mdouble, 6> getCentralParticleSizeMomenta() const;
    
    std::array<Mdouble, 6> getStandardisedParticleSizeMomenta() const;
    
    void outputStandardisedParticleSizeMomenta(std::ostream& os) const;
    
    static bool evaluateFixedParticles()
    {
        return false;
    }
    
    /*!
     * A bool that determines if the derivative of the CG function has
     * to be computed
     */
    static bool isDifferentialField()
    {
        return false;
    }

private:
    
    /*!
     * Particle volume fraction, computed as the sum over all particles i
     * \f[\nu(\vec r,t)=\sum_i V_i \phi(\vec r,\vec r_i),\f]
     * with particle volume V_i and cg function \f$\phi(\vec r,\vec r_i)\f$,
     * see CGFunctions::Gauss::evaluateCGFunction.
     */
    Mdouble volumeFraction_;
    
    /*!
     * (Mass) density, computed as the sum over all particles i
     * \f[\rho(\vec r,t)=\sum_i m_i \phi(\vec r,\vec r_i)\f]
     * with particle mass m_i and cg function \f$\phi(\vec r,\vec r_i),\f$,
     * see CGFunctions::Gauss::evaluateCGFunction.
     */
    Mdouble density_;
    
    /*!
     * Momentum, computed as the sum over all particles i
     * \f[\vec j(\vec r,t)=\sum_i m_i \vec v_i\phi(\vec r,\vec r_i),\f]
     * with particle momentum \f$m_i\vec v_i\f$ and cg function \f$\phi(\vec r,\vec r_i)\f$,
     * see CGFunctions::Gauss::evaluateCGFunction.
     *
     * Velocity can be calculated in a post-processing step from momentum and
     * density as
     * \f[\vec V(\vec r,t) = \frac{\vec j(\vec r,t)}{\rho(\vec r,t)}.\f]
     */
    Vec3D momentum_;
    
    /*!
     * Momentum, computed as the sum over all particles i
     * \f[\mathbf k(\vec r,t)=\sum_i m_i \vec v_i \otimes \vec v_i\phi(\vec r,\vec r_i),\f]
     * with particle momentum flux  \f$m_i\vec v_i \otimes \vec v_i\f$ and cg function
     * \f$\phi(\vec r,\vec r_i)\f$, see CGFunctions::Gauss::evaluateCGFunction.
     *
     * Kinetic stress can be calculated in a post-processing step from momentum
     * flux and density as
     * \f[\mathbf \sigma^{k} = \mathbf k - \rho \vec V \otimes  \vec V.\f]
     */
    MatrixSymmetric3D momentumFlux_;
    
    /*!
     * Contact stress, computed as the sum over all contacts between particle i
     * and particle/wall/fixed particle j
     * \f[\mathbf \sigma^{c}(\vec r,t)=\sum_{ij} \vec f_{ij} \otimes \vec l_{ij} \psi(\vec r,\vec r_i,\vec r_j),\f]
     * with contact force \f$\vec f_{ij}\f$,
     * branch vector \f$\vec l_{ij}= \vec c_{ij}-\vec r_i\f$,
     * particle position \f$\vec r_i\f$, contact point \f$\vec c_{ij}\f$,
     * and cg line integral \f$\psi(\vec r,\vec r_i,\vec r_j)\f$,
     * see CGFunctions::Gauss::evaluateCGIntegral.
     */
    Matrix3D contactStress_;
    
    /*!
     * Interaction force density, computed as the sum over all contacts between
     * particle i and (external) wall/fixed particle j
     * \f[\vec{IFD}(\vec r,t)=\sum_{ij} f_{ij} \phi(\vec r,\vec c_{ij}),\f]
     * with contact force \f$\vec f_{ij}\f$
     * and cg function \f$\phi(\vec r,\vec r_i)\f$,
     * see CGFunctions::Gauss::evaluateCGFunction.
     */
    Vec3D interactionForceDensity_;
    
    /*!
     * Density of particle size, and powers of particle size,
     * \f[ \vec{PS}_k(\vec r,t)= \sum_i r_i^k \phi(\vec r),\f]
     * with radius \f$r_i\f$.
     *
     * Used to compute the first five moments of the particle size distribution,
     * \f[ <r^k>(\vec r,t)= \vec{PS}_k / \vec{PS}_0, \f]
     * with \f$\vec{PS}_0\f$ denoting the number density.
     */
    std::array<Mdouble, 6> particleSizeDensity_;
    
};
    
}
#endif
