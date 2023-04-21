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

#include "Species/NormalForceSpecies/HertzianViscoelasticNormalSpecies.h"
#include "Species/NormalForceSpecies/HeatFluidCoupledSpecies.h"
#include "InteractionHandler.h"
#include "Particles/HeatFluidCoupledParticle.h"
#include "Interactions/BaseInteraction.h"
#include "Species/ParticleSpecies.h"
#include "ParticleHandler.h"
#include "DPMBase.h"


//todo Does mass and interaction radius change when a liquid film is added?

/*!
 * \details Particle read function. Has an std::istream as argument, from which
 *          it extracts the radius_, invMass_ and invInertia_, respectively.
 *          From these the mass_ and inertia_ are deduced. An additional set of
 *          properties is read through the call to the parent's method
 *          ThermalParticle::read().
 * \param[in,out] is    input stream with particle properties.
 */
void HeatFluidCoupledParticle::read(std::istream& is)
{
    ThermalParticle::read(is);
    std::string dummy;
    is >> dummy >> liquidVolume_;
}

std::string HeatFluidCoupledParticle::getNameVTK(unsigned i) const
{
    if (i==1)
        return "liquidFilmVolume";
    else if (i==2)
        return "liquidBridgeVolume";
    else if (i==0)
        return "fullLiquidVolume";
    else /* i=3 */
        return "temperature";
}

std::vector<Mdouble> HeatFluidCoupledParticle::getFieldVTK(unsigned i) const
{
    if (i==1) {
        return std::vector<Mdouble>(1, liquidVolume_);
    } else if (i==2 || i==0) {
        Mdouble fullLiquidVolume = (i==2)?0:liquidVolume_;
        for (auto k : getInteractions()) {
            auto j = dynamic_cast<LiquidMigrationWilletInteraction*>(k);
            if (j && j->getLiquidBridgeVolume()) {
                fullLiquidVolume += 0.5*j->getLiquidBridgeVolume();
//            } else {
//                logger(WARN,"All contacts of % need to be LiquidMigrationWilletInteraction",i);
            }
        }
        return std::vector<Mdouble>(1, fullLiquidVolume);
    } else {
        return std::vector<Mdouble>(1, temperature_);
    }
}

void HeatFluidCoupledParticle::actionsAfterTimeStep()
{
    //Runge–Kutta method
    double dt=getHandler()->getDPMBase()->getTimeStep();
    double k11 = f1(liquidVolume_,temperature_);
    double k21 = f2(liquidVolume_,temperature_);
    double k12 = f1(liquidVolume_+dt*0.5*k11,temperature_+dt*0.5*k21 );
    double k22 = f2(liquidVolume_+dt*0.5*k11,temperature_+dt*0.5*k21);
    double k13 = f1(liquidVolume_+dt*0.5*k12,temperature_+dt*0.5*k22);
    double k23 = f2(liquidVolume_+dt*0.5*k12,temperature_+dt*0.5*k22);
    double k14 = f1(liquidVolume_+dt*k13,temperature_+dt*k23);
    double k24 = f2(liquidVolume_+dt*k13,temperature_+dt*k23);
    liquidVolume_ = liquidVolume_+dt*(k11+2*k12+2*k13+k14)/6.0;
    temperature_ = temperature_+dt*(k21+2*k22+2*k23+k24)/6.0;
}


double HeatFluidCoupledParticle::f1(double liquidVolume_,double temperature_){
    //dynamic_cast should change (The method to access species of HeatFluidCoupledSpecies should change)
    auto species = dynamic_cast<const HeatFluidCoupledSpecies<HertzianViscoelasticNormalSpecies>*>(getSpecies());
    double delE_vb=-constants::R*species->getAmbientTemperature()*log(species->getAmbientHumidity());                     //Activation energy
    double delE_v=delE_vb*(liquidVolume_*species->getLiquidDensity()/getMass()-species->getAmbientEquilibriumMoistureContent());           //delE_v=delE_vb*f(Y-ambientEquilibriumMoistureContent_);  f(Y-ambientEquilibriumMoistureContent_)=Y-ambientEquilibriumMoistureContent_;    f is the function of difference in the moisture content
    double rho_vsat_temperature_=4.844e-9*pow((temperature_-273),4)-1.4807e-7*pow((temperature_-273),3)+2.6572e-5*pow((temperature_-273),2)-4.8613e-5*(temperature_-273)+8.342e-3;       //The saturated vapour concentration at the surface temperature temperature_ (kg/m^3)
    double rho_vs=exp(-delE_v/(constants::R*temperature_))*rho_vsat_temperature_;     //The vapour concentrations at the particle medium interface (kg/m^3)
    double S_fluid=-species->getMassTransferCoefficient()*getSurfaceArea()*(rho_vs-species->getAmbientVapourConcentration());             //The exchanged moisture with surrounding drying medium

    return (S_fluid)/species->getLiquidDensity();
}

double HeatFluidCoupledParticle::f2(double liquidVolume_,double temperature_){
    //dynamic_cast should change (The method to access species of HeatFluidCoupledSpecies should change)
    auto species = dynamic_cast<const HeatFluidCoupledSpecies<HertzianViscoelasticNormalSpecies>*>(getSpecies());
    double delE_vb=-constants::R*species->getAmbientTemperature()*log(species->getAmbientHumidity());                       //Activation energy
    double delE_v=delE_vb*(liquidVolume_*species->getLiquidDensity()/getMass()-species->getAmbientEquilibriumMoistureContent());             //delE_v=delE_vb*f(Y-ambientEquilibriumMoistureContent_);  f(Y-ambientEquilibriumMoistureContent_)=Y-ambientEquilibriumMoistureContent_;    f is the function of difference in the moisture content
    double rho_vsat_temperature_=4.844e-9*pow((temperature_-273),4)-1.4807e-7*pow((temperature_-273),3)+2.6572e-5*pow((temperature_-273),2)-4.8613e-5*(temperature_-273)+8.342e-3;       //The saturated vapour concentration at the surface temperature temperature_ (kg/m^3)
    double rho_vs=exp(-delE_v/(constants::R*temperature_))*rho_vsat_temperature_;      // The vapour concentrations at the particle medium interface (kg/m^3)
    double S_fluid=-species->getMassTransferCoefficient()*getSurfaceArea()*(rho_vs-species->getAmbientVapourConcentration());              //The exchanged moisture with surrounding drying medium

    return (species->getLatentHeatVaporization()*(1+species->getEvaporationCoefficientA()*exp((species->getEvaporationCoefficientB()*species->getLiquidDensity()/getMass())*liquidVolume_))*(S_fluid))/(getMass()*species->getHeatCapacity());
}

