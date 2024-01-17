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

#include "MeltableInteraction.h"
#include "Species/NormalForceSpecies/MeltableNormalSpecies.h"
#include "Particles/MeltableParticle.h"

void MeltableInteraction::write(std::ostream& os UNUSED) const
{
    BaseInteraction::write(os);
    os << " bondingOverlap " << bondingOverlap_;
    os << " neckRadius2 " << neckRadius2_;
}

void MeltableInteraction::read(std::istream& is UNUSED)
{
    BaseInteraction::read(is);
    std::string dummy;
    is >> dummy >> bondingOverlap_;
    is >> dummy >> neckRadius2_;
}

void MeltableInteraction::computeNormalForce() {

    setRelativeVelocity(getP()->getVelocityAtContact(getContactPoint()) - getI()->getVelocityAtContact(getContactPoint()));
    setNormalRelativeVelocity(Vec3D::dot(getRelativeVelocity(),getNormal()));

    using mathsFunc::square;
    using constants::pi;

    auto pParticle = dynamic_cast<MeltableParticle*>(getP());
    auto iParticle = dynamic_cast<MeltableParticle*>(getI());
    logger.assert_debug(pParticle || iParticle,"Particles need to be of type MeltableParticle");

    double normalForce = 0;
    // In case two particles in contacts
    if (getOverlap() >= 0)
    {
        const MeltableNormalSpecies* species = getMeltableNormalSpecies();
        // Thermal conduction
        // pi*a^2/|rij|
        //Mdouble neckRadius = std::min(sqrt(2.0*harmonicMeanRadius * getOverlap()),
        //                              constants::cbrt_2*harmonicMeanRadius);
        if (std::isnan(neckRadius2_)) {
            neckRadius2_ = 2.0*getEffectiveRadius()*getOverlap();
        }
        Mdouble contactAreaOverDist = pi * neckRadius2_ / getDistance();
        double iTemperature = iParticle
                ? iParticle->getTemperature()
                : (species->getWallTemperature()<0 ? pParticle->getTemperature() : species->getWallTemperature());
        double thermalConductionP_ = species->getThermalConductivityCoefficient() * contactAreaOverDist * (iTemperature - pParticle->getTemperature());
        pParticle->addHeat(thermalConductionP_);
        if (iParticle) iParticle->addHeat(-thermalConductionP_);

        // solid radius
        Mdouble solidRadiusP = pParticle->getSolidRadius();
        Mdouble solidRadiusI = iParticle?iParticle->getSolidRadius():solidRadiusP;
        Mdouble meanSolidRadius = 0.5*(solidRadiusP+solidRadiusI);
        Mdouble moltenLayerThicknessP = pParticle->getRadius() - solidRadiusP;
        Mdouble moltenLayerThicknessI = iParticle?iParticle->getRadius() - solidRadiusI:0.0;
        Mdouble solidOverlap = getOverlap() - moltenLayerThicknessP - moltenLayerThicknessI;
        Mdouble solidContactRadius = sqrt(meanSolidRadius * solidOverlap - 0.25*solidOverlap*solidOverlap);

        // temperature dependent elastic modulus:
        //Mdouble tempDependentElasticModulus = (effectiveElasticModulus / 2.0) * (1 + tanh(-(species->getMeltingTemperature() - aveTemp) / (species->getDeltaT())));

        // if melt -> vis + surfaceTension forces
        if (moltenLayerThicknessI> 0 | moltenLayerThicknessP > 0)
        {
            const Mdouble harmonicMeanRadius = 2.0*getEffectiveRadius();
            const Mdouble neckRadius = getNeckRadius();
            //neckRadius = std::sqrt(harmonicMeanRadius*getOverlap());
            const Mdouble meanTemperature = (pParticle->getTemperature() + iTemperature) / 2.0;
            const Mdouble jagotaTerm = harmonicMeanRadius / getDistance(); //1.0;
            const Mdouble viscosity = species->getRefViscosity()
                    * mathsFunc::exp(species->getActivationEnergy() / (constants::R * meanTemperature));
            const Mdouble viscousForce = 3.0 * pi * viscosity * neckRadius * jagotaTerm * getNormalRelativeVelocity();
            const Mdouble surfaceTensionForce = (3.0/2.0) * 0.83 * pi * species->getSurfaceTension() * neckRadius;
            normalForce -= viscousForce;
            normalForce -= surfaceTensionForce;
            //std::cout << viscousForce << " " << surfaceTensionForce << std::endl;
        }

        // 1. solid-solid
        if (solidOverlap > 0) {
            const Mdouble effectiveElasticModulus = species->getEffectiveElasticModulus();
            const Mdouble bondingContactRadius = sqrt(meanSolidRadius * getBondingOverlap());//-0.25*getBondingOverlap()*getBondingOverlap());
            if (solidOverlap > getBondingOverlap()) {
                //Calculate the damping coefficient as: d =  2*eta*sqrt(m_eff * stiffness)
                Mdouble stiffness = 2.0 * effectiveElasticModulus * solidContactRadius;
                Mdouble dissipationCoefficient =
                        2.0 * species->getDissipation() * sqrt(2.0 * getEffectiveMass() * stiffness);
                Mdouble elasticForce = 4.0 / 3.0 * effectiveElasticModulus
                        * (solidContactRadius * solidOverlap - bondingContactRadius * getBondingOverlap());
                Mdouble dampingForce = dissipationCoefficient * getNormalRelativeVelocity();
                normalForce += elasticForce - dampingForce;
                //std::cout << "s " << solidOverlap << " " << elasticForce << " " << dampingForce << std::endl;
            }
            else {
                //Calculate the damping coefficient as: d =  2*eta*sqrt(m_eff * stiffness)
                Mdouble stiffness = 2.0 * effectiveElasticModulus * bondingContactRadius;
                Mdouble dissipationCoefficient = 2.0 * species->getDissipation()
                        * sqrt(2.0 * getEffectiveMass() * stiffness);
                Mdouble elasticForce = 2.0 * effectiveElasticModulus
                        * bondingContactRadius * (solidOverlap - getBondingOverlap());
                Mdouble dampingForce = dissipationCoefficient * getNormalRelativeVelocity();
                normalForce += elasticForce - dampingForce;
                //std::cout << "s " << elasticForce << " " << dampingForce << std::endl;
            }
        }
    }
    setAbsoluteNormalForce(std::abs(normalForce));
    setForce(getNormal() * normalForce);
    setTorque(Vec3D(0.0, 0.0, 0.0));
}

void MeltableInteraction::actionsAfterTimeStep()
{
    // bonding radius
    auto pParticle = dynamic_cast<const MeltableParticle*>(getP());
    auto iParticle = dynamic_cast<const MeltableParticle*>(getI());

    // calculate the bonding overlap
    const Mdouble solidRadiusP = pParticle->getSolidRadius();
    const Mdouble solidRadiusI = iParticle?iParticle->getSolidRadius():solidRadiusP;
    const Mdouble moltenLayerThicknessP = pParticle->getRadius() - solidRadiusP;
    const Mdouble moltenLayerThicknessI = iParticle?iParticle->getRadius() - solidRadiusI:0.0;
    const Mdouble solidOverlap = getOverlap() - moltenLayerThicknessP - moltenLayerThicknessI;
    const Mdouble timeStep = getHandler()->getDPMBase()->getTimeStep();
    if (solidOverlap >= 0) {
        Mdouble meltRate = pParticle->getMeltRate(solidRadiusP)
                + (iParticle?iParticle->getMeltRate(solidRadiusI):0.0);
        if (meltRate < 0) {
            addBondingOverlap(-meltRate*timeStep);
        }
    } else {
        setBondingOverlap(0.0);
    }

    // calculate neck radius
    using mathsFunc::square;
    const Mdouble radius = 2.0*getEffectiveRadius();
    const Mdouble dOverlap = timeStep*getOverlapGrowthRate();
    if (moltenLayerThicknessP==0 and moltenLayerThicknessI==0) {
        if (getBondingOverlap()>0) {
            addNeckRadius2(radius*dOverlap);
        }
    } else {
        addNeckRadius2(2.0*radius*dOverlap);
    }
}

Mdouble MeltableInteraction::getElasticEnergy() const
{
    ///\todo
    // CALCULATE FOR THE MODEL:
    //
/*    if (getOverlap() >= 0)
    {
        //return 8. / 15. * getSpecies()->getElasticModulus() * std::sqrt(getEffectiveRadius() * getOverlap()) * mathsFunc::square(getOverlap());
        //
        if ((getSolidOverlap() > getBondingOverlap() & getSolidOverlap() > 0))
        {
            return ((4.0 / 3.0) * (2.0 / 5.0) * getSpecies()->getElasticModulus() * getSolidContactRadius() * mathsFunc::square(getSolidOverlap()))
                   - ((4.0 / 3.0) * getSpecies()->getElasticModulus() * getBondingContactRadius() * getBondingOverlap() * getSolidOverlap());
        }
        else if (getSolidOverlap() <= getBondingOverlap() & getSolidOverlap() > 0)
        {
            return ((4.0 / 3.0) * (3.0 / 4.0) * getSpecies()->getElasticModulus() * getBondingContactRadius() * mathsFunc::square(getSolidOverlap()))
                   - ((4.0 / 3.0) * (3.0 / 2.0) * getSpecies()->getElasticModulus() * getBondingContactRadius() * getBondingOverlap() * getSolidOverlap());
        }
        else
            return 0.0;
    }
    else
    {
        return 0.0;
    }*/
return 0.0;
}

const MeltableNormalSpecies* MeltableInteraction::getMeltableNormalSpecies() const
{
    return dynamic_cast<const MeltableNormalSpecies*>(getBaseSpecies()); //downcast
}