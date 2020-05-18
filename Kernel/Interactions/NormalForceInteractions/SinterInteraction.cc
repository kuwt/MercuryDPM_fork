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


#include "SinterInteraction.h"
#include "InteractionHandler.h"
#include "Particles/BaseParticle.h"
#include "DPMBase.h"
#include <Species/NormalForceSpecies/SinterNormalSpecies.h>
#include <Particles/ThermalParticle.h>

/*!
 * \param[in] P
 * \param[in] I
 * \param[in] timeStamp
 */
SinterInteraction::SinterInteraction(BaseInteractable* P, BaseInteractable* I, unsigned timeStamp)
        : BaseInteraction(P, I, timeStamp)
{
    plasticOverlap_ = 0;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"SinterInteraction::SinterInteraction() finished"<<std::endl;
#endif
}

/*!
 * \param[in] p 
 */
SinterInteraction::SinterInteraction(const SinterInteraction& p)
        : BaseInteraction(p)
{
    plasticOverlap_ = p.plasticOverlap_;
#ifdef DEBUG_CONSTRUCTOR
    std::cout<<"SinterInteraction::SinterInteraction(const SinterInteraction &p finished"<<std::endl;
#endif
}

SinterInteraction::SinterInteraction()
= default;

/*!
 *
 */
SinterInteraction::~SinterInteraction()
{
#ifdef DEBUG_DESTRUCTOR
    std::cout<<"SinterInteraction::~SinterInteraction() finished"<<std::endl;
#endif
}

/*!
 * \details Calls the write function of BaseInteraction().
 * \param[in,out] os
 */
void SinterInteraction::write(std::ostream& os) const
{
    BaseInteraction::write(os);
    os << " plasticOverlap " << plasticOverlap_;
}

/*!
 * \details Calls the read function of BaseInteraction().
 * \param[in,out] is
 */
void SinterInteraction::read(std::istream& is)
{
    BaseInteraction::read(is);
    std::string dummy;
    is >> dummy >> plasticOverlap_;
}

/*!
 * \return std::string
 */
std::string SinterInteraction::getBaseName() const
{
    return "Sinter";
}

/*!
 *
 */
void SinterInteraction::computeNormalForce()
{
    // Compute the relative velocity vector of particle P w.r.t. I
    setRelativeVelocity(
            getP()->getVelocityAtContact(getContactPoint()) - getI()->getVelocityAtContact(getContactPoint()));
    // Compute the projection of vrel onto the normal (can be negative)
    setNormalRelativeVelocity(Vec3D::dot(getRelativeVelocity(), getNormal()));
    
    if (getOverlap() > 0) //if contact forces
    {
        const SinterNormalSpecies* species = getSpecies();
        // calculate the effective diameter, equal to the radius for two equal-sized particles
        const Mdouble effectiveDiameter = 2.0 * getEffectiveRadius();
        //calculate the overlap above which the max. unloading stiffness 
        //becomes active (the 'fluid branch')
        const Mdouble deltaStar = species->getPenetrationDepthMax() * effectiveDiameter;
        //compute the rate d(k2/k1)/d(delta0)
        const Mdouble dk = (species->getUnloadingStiffnessMax() / species->getLoadingStiffness() - 1.0) / deltaStar;
        
        //increase max overlap if necessary
        //k2*(d-d0)=k1*d, k2=k1*(1+dk*d)
        //(1+dk*d)*(d-d0)=d
        //dk*d^2=(1+dk*d)*d0
        //d0=d/(1+1/(dk*d))
        const Mdouble minPlasticOverlap = std::min(getOverlap() / (1.0 + 1.0 / (dk * getOverlap())), deltaStar);
        ///\todo basing the stiffness increase on the plastic overlap does not work:
        ///k2*(d-d0)=k1*d, k2=k1*(1+dk*d0)
        ///(1+dk*d0)*(d-d0)=d
        ///d+dk*d*d0-d0-dk*d0*d0=d
        ///-dk*d0^2+(dk*d-1)*d0+d=d
        ///d0=(dk*d-1)/dk
        //const Mdouble minPlasticOverlap = std::min(getOverlap()-1.0/dk,deltaStar);
        //logger(INFO,"minPlasticOverlap % overlap % dk % M %",minPlasticOverlap,getOverlap(),dk, dk*getOverlap());
        plasticOverlap_ = std::max(minPlasticOverlap, plasticOverlap_);
        
        //calculate the unloading stiffness (only linear in maxOverlap, not equilibriumOverlap))
        const Mdouble unloadingStiffness = species->getLoadingStiffness() * (1.0 + dk * plasticOverlap_);
        
        //compute elastic force
        Mdouble normalForce = unloadingStiffness * (getOverlap() - plasticOverlap_);
        
        //add cohesive forces (distinct from sintering)
        Mdouble nonSinterAdhesiveForce = -species->getCohesionStiffness() * getOverlap();
        if (normalForce < nonSinterAdhesiveForce)
        {
            plasticOverlap_ = (1.0 + species->getCohesionStiffness() / unloadingStiffness) * getOverlap();
            normalForce = nonSinterAdhesiveForce;
        }
        
        //add dissipative force (distinct from sintering)
        normalForce -= species->getDissipation() * getNormalRelativeVelocity();
        
        //Here comes the sintering effect:
        Mdouble adhesiveForce = species->getSinterAdhesion() * effectiveDiameter;
        
        //now set the interaction force equal to this normal force (friction and adhesive forces will be added later)
        setForce(getNormal() * ((normalForce - adhesiveForce)));
        setTorque(Vec3D(0.0, 0.0, 0.0));
        //used for tangential force calculations; don't add adhesive force components
        setAbsoluteNormalForce(std::abs(normalForce));
        
        //increase plastic overlap due to sintering
        DPMBase* dpmBase = getHandler()->getDPMBase();
        Mdouble rateOverlap;
        // sinter adhesion force fa=sinterAdhesion_*radius in sinter rate:
        if (species->getSinterType() == SINTERTYPE::PARHAMI_MCKEEPING)
        {
            rateOverlap = normalForce * species->getSinterRate() /
                          (0.375 * species->getSinterAdhesion() * mathsFunc::square(
                                  getOverlap() / effectiveDiameter)); ///\todo adhesive force only, or add normalForce?
            if (species->getSinterRate() == 0) rateOverlap = 0;
        }
        else if (species->getSinterType() == SINTERTYPE::CONSTANT_RATE)
        {
            rateOverlap = 2.0 * normalForce * species->getSinterRate() / species->getSinterAdhesion();
            if (species->getSinterRate() == 0) rateOverlap = 0;
        }
        else if (species->getSinterType() == SINTERTYPE::TEMPERATURE_DEPENDENT_FRENKEL)
        {
            ThermalParticle* tp = dynamic_cast<ThermalParticle*>(getP());
            ThermalParticle* ti = dynamic_cast<ThermalParticle*>(getI());
            logger.assert(tp && ti,
                          "warning contact partners have to be ThermalParticle's if this sinter species is used");
            double temperature =
                    2.0 * tp->getTemperature() * ti->getTemperature() / (tp->getTemperature() + ti->getTemperature());
            rateOverlap = 2.0 * normalForce * species->getTemperatureDependentSinterRate(temperature) /
                          species->getSinterAdhesion();
        }
        else
        {
            rateOverlap = 0;
            //missing: add the sintering model 'modified Frenkel' of the Pokula paper
        }
        plasticOverlap_ = std::max(0.0, std::min(deltaStar, plasticOverlap_ + rateOverlap * dpmBase->getTimeStep()));
        
        /*//change particle radius by dr
        Mdouble dr;
        BaseParticle* PParticle = dynamic_cast<BaseParticle*>(getP());
        BaseParticle* IParticle = dynamic_cast<BaseParticle*>(getI());
        if (dpmBase->getSystemDimensions()==2) {
            //2D: increase the radius of each particle such that the particle area
            //increases by the same amount that the contact area decreases
            //Particle circumference C = 2 pi r increased by dr => dA = 2 pi r dr
            //Contact line L = 2*sqrt(2*r*o) indented by do/2 => dA = sqrt(2*r*o) do
            //Thus, dr = sqrt(0.5*o/r)/pi do.
            dr = sqrt(0.5*plasticOverlap_/effectiveDiameter)/3.14 *doverlap;
        } else {
            //3D: increase the radius of each sphere such that the particle volume
            //increases by the same amount that the contact volume decreases
            //Particle surface area S = 4 pi r^2 increased by dr => dA = 4 pi r^2 dr
            //Contact area L = pi 2*r*o indented by do/2 => dA = pi r o do
            //Thus, dr = 0.25*o/r do
            dr = 0.25*plasticOverlap_/effectiveDiameter *doverlap;
        }
        if (PParticle==nullptr) { //if P is a wall
            IParticle->setRadius(IParticle->getRadius()+dr);//should be twice that amount
        } else if (IParticle==nullptr) { //if I is a wall
            PParticle->setRadius(PParticle->getRadius()+dr);
        } else { //if both P and I are particles
            PParticle->setRadius(PParticle->getRadius()+dr);
            IParticle->setRadius(IParticle->getRadius()+dr);
        }*/
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
Mdouble SinterInteraction::getElasticEnergy() const
{
    if (getOverlap() > 0)
        return 0.5 * (getSpecies()->getLoadingStiffness() * mathsFunc::square(getOverlap()));
    else
        return 0.0;
    ///\todo TW this is not correct; we should count the return energy
}

/*!
 * \return const SinterNormalSpecies*
 */
const SinterNormalSpecies* SinterInteraction::getSpecies() const
{
    return static_cast<const SinterNormalSpecies*>(getBaseSpecies()->getNormalForce());
;
}

/*!
 * \return Mdouble plasticOverlap_
 */
Mdouble SinterInteraction::getPlasticOverlap() const
{
    return plasticOverlap_;
}

/*!
 * \param[in] maxOverlap
 */
void SinterInteraction::setPlasticOverlap(const Mdouble plasticOverlap)
{
    plasticOverlap_ = plasticOverlap;
}

/*!
 * \return Mdouble
 */
Mdouble SinterInteraction::getUnloadingStiffness() const
{
    const SinterNormalSpecies* species = getSpecies();
    Mdouble effectiveDiameter = 2.0 * getEffectiveRadius();
    Mdouble deltaMaxFluid = species->getPenetrationDepthMax() * effectiveDiameter /
                            (1.0 - species->getLoadingStiffness() / species->getUnloadingStiffnessMax());
    if (getOverlap() > deltaMaxFluid)
        return species->getUnloadingStiffnessMax();
    else
        return species->getLoadingStiffness() + (species->getUnloadingStiffnessMax() - species->getLoadingStiffness()) *
                                                getPlasticOverlap() / deltaMaxFluid;
}

unsigned SinterInteraction::getNumberOfFieldsVTK() const
{
    return 2;
}

std::string SinterInteraction::getTypeVTK(unsigned i) const
{
    return "Float32";
}

std::string SinterInteraction::getNameVTK(unsigned i) const
{
    if (i == 0)
        return "plasticOverlap";
    else
        return "neckRadius";
}

std::vector<Mdouble> SinterInteraction::getFieldVTK(unsigned i) const
{
    if (i == 0)
        return std::vector<Mdouble>(1, plasticOverlap_);
    else
        return std::vector<Mdouble>(1, sqrt(2.0 * getEffectiveRadius() * plasticOverlap_));
}


//set xlabel 't'; set ylabel 'f_n = f_ep-f_a'; p 'SingleParticleCT.fstat' u 1:9 title 'f_g=f_a', 'SingleParticleCT_noGravity.fstat' u 1:9 title 'f_g=0'/

//set xlabel 't'; set ylabel 'delta/r'; p 'SingleParticleCT.fstat' u 1:($7/1e-6) title 'f_g=f_a', 'SingleParticleCT_noGravity.fstat' u 1:($7/1e-6) title 'f_g=0'
