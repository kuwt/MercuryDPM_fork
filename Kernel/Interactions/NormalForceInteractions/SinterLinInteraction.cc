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


#include "SinterLinInteraction.h"
#include "Particles/BaseParticle.h"
#include "DPMBase.h"
#include <iomanip>
#include <Species/NormalForceSpecies/SinterLinNormalSpecies.h>
#include <cmath>    // std::max
using constants::pi;

/*!
 * \param[in] P
 * \param[in] I
 * \param[in] timeStamp
 */
SinterLinInteraction::SinterLinInteraction(BaseInteractable* P, BaseInteractable* I,
                                                                           unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp)
{
    plasticOverlap_ = 0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearPlasticViscoelasticInteraction::LinearPlasticViscoelasticInteraction() finished"<<std::endl;
#endif
}

//used for mpi
SinterLinInteraction::SinterLinInteraction()
        : BaseInteraction()
{
    plasticOverlap_ = 0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearPlasticViscoelasticInteraction::LinearPlasticViscoelasticInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] p
 */
SinterLinInteraction::SinterLinInteraction(
        const SinterLinInteraction& p)
        : BaseInteraction(p)
{
    plasticOverlap_ = p.plasticOverlap_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"LinearPlasticViscoelasticInteraction::LinearPlasticViscoelasticInteraction(const LinearPlasticViscoelasticInteraction &p finished"<<std::endl;
#endif
}

/*!
 *
 */
SinterLinInteraction::~SinterLinInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"SinterLinInteraction::~SinterLinInteraction() finished"<<std::endl;
#endif
}

/*!
 * \details Calls the write function of BaseInteraction().
 * \param[in,out] os
 */
void SinterLinInteraction::write(std::ostream& os) const
{
    BaseInteraction::write(os);
    os << " plasticOverlap " << plasticOverlap_;
}

/*!
 * \details Calls the read function of BaseInteraction().
 * \param[in,out] is
 */
void SinterLinInteraction::read(std::istream& is)
{
    BaseInteraction::read(is);
    std::string dummy;
    is >> dummy >> plasticOverlap_;
    //helpers::readOptionalVariable<Mdouble>(is, "maxOverlap", maxOverlap_);
}

/*!
 * \return std::string
 */
std::string SinterLinInteraction::getBaseName() const
{
    return "Sinter";
}

/*!
 *
 */
void SinterLinInteraction::computeNormalForce()
{

    // Compute the relative velocity vector of particle P w.r.t. I
    setRelativeVelocity(
            getP()->getVelocityAtContact(getContactPoint()) - getI()->getVelocityAtContact(getContactPoint()));
    // Compute the projection of vrel onto the normal (can be negative)
    setNormalRelativeVelocity(Vec3D::dot(getRelativeVelocity(), getNormal()));

    if (getOverlap() > 0) //if contact forces
    {
        const SinterLinNormalSpecies* species = getSpecies();

        // calculate the effective diameter, equal to the radius for two equal-sized particles
        const Mdouble effectiveDiameter = 2.0 * getEffectiveRadius();

        //[1] Compute  delta fluid of equilibrium. In this model, it is the same than delta*
        const Mdouble d_fluid_0 = (species->getUnloadingStiffnessMax()
                                   / (species->getUnloadingStiffnessMax() - species->getLoadingStiffness()))
                                  * species->getPenetrationDepthMax() * effectiveDiameter;

        //[2] Compute the rate d(k2/k1)/d(delta0). To obtain this parameter, the linear relationship
        //between unloading stiffness and maximum overlap is used.
        const Mdouble dk = (species->getUnloadingStiffnessMax()/species->getLoadingStiffness() - 1.0)*(1.0/d_fluid_0);

        //Increase max overlap if necessary
        //Here, two relationships are used. The linear relationship between unloadingstiffness max and
        // the equilibrium overlap based on the unloading stiffness.
        //k2 = k1*(1+dk*d), k1*d = k2*(d-d0)
        //k1*d = k1*(1 + dk*d)*(d-d0)
        //d = d - d0 +dk*d^2 -dk*d*d0
        //dk*d^2 = dk *d * d0 + d0
        //dk*d^2 = (dk*d + 1)*d0
        //dkd^2/(dk*d + 1) = d0
        //d0 = d/(1 + 1/(dk*d))

        //[3] Compute the equilibirum overlap is:
        const Mdouble d0 = getOverlap()/(1.0 + 1.0/(dk*getOverlap()));

        const Mdouble minPlasticOverlap = std::min(d0,d_fluid_0);

        //[4] Determine the plastic overlap
        plasticOverlap_ = std::max(minPlasticOverlap,plasticOverlap_);

        //[5] Compute the unloading Stiffness \hat{k2}.
        const Mdouble unloadingStiffness = species->getLoadingStiffness() * (1.0 + dk * plasticOverlap_);

        //[6] Compute the elastic force
        Mdouble normalForce = unloadingStiffness * (getOverlap() - plasticOverlap_);

        //[7] Add cohesive force (distinct from sinteirng)
        //Decrease plastic overlap if necessary
        Mdouble nonSinterAdhesiveForce = -species->getCohesionStiffness() * getOverlap();

        if (normalForce < nonSinterAdhesiveForce)
        {
            plasticOverlap_ = (1.0 + species->getCohesionStiffness() / unloadingStiffness) * getOverlap();
            normalForce = nonSinterAdhesiveForce;
        }

        //[[8] Add dissipative force
        normalForce -= species->getDissipation() * getNormalRelativeVelocity();

        //[9] Sintering effect as adhesive force:
        Mdouble adhesiveForce = species->getSinterAdhesion() * effectiveDiameter;
//        Mdouble adhesiveForce = species->getSinterAdhesion();

        //[10] now set the interaction force equal to this normal force (friction and adhesive forces will be added later)
        setForce(getNormal() * ((normalForce - adhesiveForce)));
        setTorque(Vec3D(0.0, 0.0, 0.0));
        //used for tangential force calculations; don't add adhesive force components
        setAbsoluteNormalForce(std::abs(normalForce));

        //[11] Checking the overlap to change sintering model
        //Before this point, the material is unrelaxed, and the contact radious is given by the JKR theory.
        //Ref: The role of viscoelastic adhesive contact in the sintering of polymeric particles
        //Author: Y.Y Lin et al.
        const Mdouble baseNum = (9.0/2.0)*pi*species->getComplianceZero()*species->getSurfTension()/(getEffectiveRadius());
        const Mdouble a0_R = std::pow(2.0*baseNum,1.0/3.0);

        const Mdouble realOverlap = std::sqrt(getOverlap()/effectiveDiameter);

        DPMBase* dpmBase = getHandler()->getDPMBase();
        Mdouble rateOverlap;
        Mdouble rateOverlap2;

        if (species->getSinterType() == SINTER_APPROACH::FRENKEL)
        {
            rateOverlap = 2.0*normalForce * species->getSinterRate() / species->getSinterAdhesion();
        }
        else if (species->getSinterType() == SINTER_APPROACH::VISCOELASTIC_CONTACT)
        {
            if ((getContactRadius())/getEffectiveRadius() < a0_R)
            {
                rateOverlap = 0.0; //normalForce / species->getSinterAdhesion();
                //ToDo: Measure the evolution of the creep compliance: Art.Contact creep compliance of viscoelastic material via nanoindentation
            } else {
                const Mdouble C1 = std::pow((63 * pow(pi, 3.0) / 16.0), 2.0 / 7.0);
                const Mdouble C2 = pow((species->getSeparationDis()/ getEffectiveRadius()), 2.0 / 7.0);
                const Mdouble C3 = species->getFluidity()* species->getSurfTension()/getEffectiveRadius();

                Mdouble time_s = std::pow(8.0*getOverlap()/(getEffectiveRadius()*C1*C2*std::pow(C3,2.0/7.0)),7.0/2.0);

                rateOverlap = getEffectiveRadius()*C1*C2*((2.0*std::pow(C3,2.0/7.0))/(7.0*std::pow(time_s,5.0/7.0))) * (normalForce / species->getSinterAdhesion());

                const Mdouble aVisc_R = std::pow(63.0 * std::pow(pi, 3.0), 1.0 / 5.0) *
                                        std::pow((1.0/8.0)*(species->getSeparationDis())/ (0.5*getEffectiveRadius()), 2.0 / 5.0);


                if (getContactRadius()/getEffectiveRadius() > aVisc_R)
                {
                    Mdouble t  = (getOverlap())*(1.0/(4.0* species->getFluidity()* species->getSurfTension()));

                    Mdouble theta = atan(std::pow((8.0* species->getFluidity()* species->getSurfTension()*t)/getEffectiveRadius(),0.5));

                    Mdouble dtheta = (8.0* species->getFluidity()* species->getSurfTension()/getEffectiveRadius());

                    dtheta *= (std::pow(2.0,-5.0/3.0))*cos(theta)*sin(theta);
                    Mdouble K1 = tan(theta)/2.0 - (sin(theta)/6.0)*((2.0*(2.0-cos(theta))+(1.0+cos(theta)))/((1.0+cos(theta)))*(2.0-cos(theta)));
                    dtheta /= std::pow(K1,2.0);

                    dtheta /= (std::pow(2.0-cos(theta),5.0/3.0))*(std::pow(1.0+cos(theta),4.0/3.0));
                    //++++++
                    rateOverlap2 = 2.0*getEffectiveRadius()*std::pow(4.0/(std::pow(1.0+cos(theta),2.0)*(2.0-cos(theta))),2.0/3.0);
                    rateOverlap2 *= sin(theta)*cos(theta)* dtheta;

                    Mdouble term1 = 4.0*std::pow(2.0,1.0/3.0)*std::pow(sin(theta),3.0)*dtheta*(cos(theta)-1.0);
                    Mdouble term2 = std::pow(cos(theta)+1.0,7.0/3.0)*std::pow(-cos(theta)+2.0,5.0/3.0);
                    rateOverlap2 -= term1/term2;
                    rateOverlap2 *= 2.0*getEffectiveRadius()*(normalForce/species->getSinterAdhesion());

                    rateOverlap = rateOverlap2;
                }
            }
        }
        else{
            rateOverlap = 0.0;
        }
        //[12] Increase plastic overlap due to sintering
        plasticOverlap_ = std::max(0.0, std::min(d_fluid_0, plasticOverlap_ + rateOverlap * dpmBase->getTimeStep()));
    }
    else
    {
        setAbsoluteNormalForce(0.0);
        setForce(Vec3D(0.0, 0.0, 0.0));
        setTorque(Vec3D(0.0, 0.0, 0.0));
    }
}
/*!
 * \return Mdouble
 */
Mdouble SinterLinInteraction::getElasticEnergy() const
{
    Mdouble energy = getOverlap() > 0 ? 0.5 * (getSpecies()->getLoadingStiffness() * mathsFunc::square(getOverlap())) : 0.0;
    if (getSpecies()->getConstantRestitution()) energy *= 2.0*getEffectiveMass();
    return energy;
    ///\todo TW this is not correct; we should count the return energy
}

/*!
 * \return const SinterLinNormalSpecies*
 */
const SinterLinNormalSpecies* SinterLinInteraction::getSpecies() const
{
    return dynamic_cast<const SinterLinNormalSpecies*>(getBaseSpecies());
}

/*!
 * \return Mdouble plasticOverlap_
 */
Mdouble SinterLinInteraction::getMaxOverlap() const
{
    return maxOverlap_;
}

/*!
 * \return Mdouble plasticOverlap_
 */
Mdouble SinterLinInteraction::getPlasticOverlap() const
{
    return plasticOverlap_;
}



/*!
 * \param[in] maxOverlap
 */
void SinterLinInteraction::setMaxOverlap(const Mdouble maxOverlap)
{
    maxOverlap_ = maxOverlap;
}

/*!
 * \param[in] maxOverlap
 */
void SinterLinInteraction::setPlasticOverlap(const Mdouble plasticOverlap)
{
    plasticOverlap_ = plasticOverlap;
}


/*!
 * \return Mdouble
 */
Mdouble SinterLinInteraction::getUnloadingStiffness() const
{

    const SinterLinNormalSpecies* species = getSpecies();
    const Mdouble effectiveDiameter = 2.0 * getEffectiveRadius();

    Mdouble d_max_fluid = (species->getUnloadingStiffnessMax() / (species->getUnloadingStiffnessMax()
                                                                  - species->getLoadingStiffness()))* species->getPenetrationDepthMax() *
                          effectiveDiameter;

    if (getOverlap() > d_max_fluid)
        return species->getUnloadingStiffnessMax();
    else
        return species->getLoadingStiffness() +
               (species->getUnloadingStiffnessMax() - species->getLoadingStiffness()) * getMaxOverlap() / d_max_fluid;
}