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

#include "MeltableNormalSpecies.h"
#include <Logger.h>

void MeltableNormalSpecies::write(std::ostream &os) const
{
    os << " elasticModulus " << elasticModulus_
       << " poissonRatio " << poissonRatio_
       << " dissipation " << dissipation_
       << " latentHeat " << latentHeat_
       << " thermalConductivityCoefficient " << thermalConductivityCoefficient_
       << " thermalConvectionCoefficient " << thermalConvectionCoefficient_
       << " materialEmissivity " << materialEmissivity_
       << " heatCapacity " << solidHeatCapacity_
       << " meltingTemperature " << meltingTemperature_
       << " ambientTemperature " << ambientTemperature_
       << " heatInput " << heatInput_
       << " heatInputFunction " << (heatInputFunction_!= nullptr)
       << " materialAbsorptivity " << materialAbsorptivity_
       << " deltaT " << deltaT_
       << " liquidHeatCapacity " << liquidHeatCapacity_
       << " thermalExpansionCoefficient " << thermalExpansionCoefficient_
       << " activationEnergy " << activationEnergy_
       << " surfaceTension " << surfaceTension_
       << " refViscosity " << refViscosity_;
}

void MeltableNormalSpecies::read(std::istream &is)
{
    std::string dummy;
    bool hasHeatInputFunction;
    is >> dummy >> elasticModulus_
       >> dummy >> poissonRatio_
       >> dummy >> dissipation_
       >> dummy >> latentHeat_
       >> dummy >> thermalConductivityCoefficient_
       >> dummy >> thermalConvectionCoefficient_
       >> dummy >> materialEmissivity_
       >> dummy >> solidHeatCapacity_
       >> dummy >> meltingTemperature_
       >> dummy >> ambientTemperature_
       >> dummy >> heatInput_
       >> dummy >> hasHeatInputFunction
       >> dummy >> materialAbsorptivity_
       >> dummy >> deltaT_
       >> dummy >> liquidHeatCapacity_
       >> dummy >> thermalExpansionCoefficient_
       >> dummy >> activationEnergy_
       >> dummy >> surfaceTension_
       >> dummy >> refViscosity_;
    if (hasHeatInputFunction) {
        logger(WARN,"Could not restart heat input function; heat input set to %", heatInput_);
    }
}

std::string MeltableNormalSpecies::getBaseName() const
{
    return "Meltable";
}

void MeltableNormalSpecies::setElasticModulus(Mdouble elasticModulus)
{
    if (elasticModulus >= 0)
        elasticModulus_ = elasticModulus;
    else logger(ERROR,"Error in setElasticModulus(%)",elasticModulus);
}

Mdouble MeltableNormalSpecies::getElasticModulus() const
{
    return elasticModulus_;
}

void MeltableNormalSpecies::setPoissonRatio(Mdouble poissonRatio)
{
    if (poissonRatio >= -1 and poissonRatio <= 0.5)
        poissonRatio_ = poissonRatio;
    else logger(ERROR,"Error in setPoissonRatio(%)",poissonRatio);
}

Mdouble MeltableNormalSpecies::getPoissonRatio() const
{
    return poissonRatio_;
}

void MeltableNormalSpecies::setDissipation(Mdouble dissipation)
{
    if (dissipation >= 0)
        dissipation_ = dissipation;
    else logger(ERROR,"Error in setDissipation(%)",dissipation);
}

Mdouble MeltableNormalSpecies::getDissipation() const
{
    return dissipation_;
}

void MeltableNormalSpecies::setLatentHeat(Mdouble latentHeat)
{
    if (latentHeat >=0)
        latentHeat_ = latentHeat;
    else logger(ERROR,"Error in setLatentHeat(%)",latentHeat);
}

Mdouble MeltableNormalSpecies::getLatentHeat() const
{
    return latentHeat_;
}

void MeltableNormalSpecies::setThermalConductivityCoefficient(Mdouble thermalConductivityCoefficient)
{
    if (thermalConductivityCoefficient >= 0)
        thermalConductivityCoefficient_ = thermalConductivityCoefficient;
    else logger(ERROR,"Error in setThermalConductivityCoefficient(%)",thermalConductivityCoefficient);
}

Mdouble MeltableNormalSpecies::getThermalConductivityCoefficient() const
{
    return thermalConductivityCoefficient_;
}

void MeltableNormalSpecies::setThermalConvectionCoefficient(Mdouble thermalConvectionCoefficient)
{
    if (thermalConvectionCoefficient >= 0)
        thermalConvectionCoefficient_ = thermalConvectionCoefficient;
    else logger(ERROR,"Error in setThermalConvectionCoefficient(%)",thermalConvectionCoefficient);
}

Mdouble MeltableNormalSpecies::getThermalConvectionCoefficient() const
{
    return thermalConvectionCoefficient_;
}

void MeltableNormalSpecies::setMaterialEmissivity(Mdouble materialEmissivity)
{
    if (materialEmissivity >= 0)
        materialEmissivity_ = materialEmissivity;
    else logger(ERROR,"Error in setMaterialEmissivity(%)",materialEmissivity);
}

Mdouble MeltableNormalSpecies::getMaterialEmissivity() const
{
    return materialEmissivity_;
}

void MeltableNormalSpecies::setSolidHeatCapacity(Mdouble solidHeatCapacity)
{
    if (solidHeatCapacity >= 0)
        solidHeatCapacity_ = solidHeatCapacity;
    else logger(ERROR,"Error in setSolidHeatCapacity(%)",solidHeatCapacity);
}

Mdouble MeltableNormalSpecies::getSolidHeatCapacity() const
{
    return solidHeatCapacity_;
}

void MeltableNormalSpecies::setAmbientTemperature(Mdouble ambientTemperature)
{
    if (ambientTemperature >= 0)
        ambientTemperature_ = ambientTemperature;
    else logger(ERROR,"Error in setAmbientTemperature(%)",ambientTemperature);
}

Mdouble MeltableNormalSpecies::getAmbientTemperature() const
{
    return ambientTemperature_;
}

void MeltableNormalSpecies::setMeltingTemperature(Mdouble meltingTemperature)
{
    if (meltingTemperature >= 0)
        meltingTemperature_ = meltingTemperature;
    else logger(ERROR,"Error in setMeltingTemperature(%)",meltingTemperature);
}

Mdouble MeltableNormalSpecies::getMeltingTemperature() const
{
    return meltingTemperature_;
}

void MeltableNormalSpecies::setHeatInput(Mdouble heatInput)
{
    heatInput_ = heatInput;
}

Mdouble MeltableNormalSpecies::getHeatInput(const BaseParticle* p) const
{
    if (heatInputFunction_) {
        return heatInputFunction_(p);
    } else {
        return heatInput_;
    }
}

//double MeltableNormalSpecies::computeMaxTimeStep(Mdouble minRadius, Mdouble minDensity, Mdouble minElasticModulus)
//{
//    return 0.5*minRadius*sqrt(minDensity/minElasticModulus);
//}

void MeltableNormalSpecies::setMaterialAbsorptivity(Mdouble materialAbsorptivity)
{
    if (materialAbsorptivity >= 0)// && materialAbsorptivity <= 1.0)
        materialAbsorptivity_ = materialAbsorptivity;
    else logger(ERROR,"Error in setMaterialAbsorptivity(%)",materialAbsorptivity);
}

Mdouble MeltableNormalSpecies::getMaterialAbsorptivity() const
{
    return materialAbsorptivity_;
}

void MeltableNormalSpecies::setDeltaT(Mdouble deltaT)
{
    if (deltaT > 0)
        deltaT_ = deltaT;
    else logger(ERROR,"Error in setDeltaT(%)",deltaT);
}

Mdouble MeltableNormalSpecies::getDeltaT() const
{
    return deltaT_;
}

void MeltableNormalSpecies::setLiquidHeatCapacity(Mdouble liquidHeatCapacity)
{
    if (liquidHeatCapacity >= 0)
        liquidHeatCapacity_ = liquidHeatCapacity;
    else logger(ERROR,"Error in setLiquidHeatCapacity(%)",liquidHeatCapacity);
}

Mdouble MeltableNormalSpecies::getLiquidHeatCapacity() const
{
    return liquidHeatCapacity_;
}

void MeltableNormalSpecies::setThermalExpansionCoefficient(Mdouble thermalExpansionCoefficient)
{
    thermalExpansionCoefficient_ = thermalExpansionCoefficient;
}

Mdouble MeltableNormalSpecies::getThermalExpansionCoefficient() const
{
    return thermalExpansionCoefficient_;
}

void MeltableNormalSpecies::setActivationEnergy(Mdouble activationEnergy)
{
    if (activationEnergy >= 0)
        activationEnergy_ = activationEnergy;
    else logger(ERROR,"Error in setSurfaceTension(%)",activationEnergy_);
}

Mdouble MeltableNormalSpecies::getActivationEnergy() const
{
    return activationEnergy_;
}

void MeltableNormalSpecies::setSurfaceTension(Mdouble surfaceTension)
{
    if (surfaceTension >= 0)
        surfaceTension_ = surfaceTension;
    else logger(ERROR,"Error in setSurfaceTension(%)",surfaceTension);
}

Mdouble MeltableNormalSpecies::getSurfaceTension() const
{
    return surfaceTension_;
}

void MeltableNormalSpecies::setRefViscosity(Mdouble refViscosity)
{
    if (refViscosity >= 0)
        refViscosity_ = refViscosity;
    else logger(ERROR,"Error in setRefViscosity(%)",refViscosity);
}

Mdouble MeltableNormalSpecies::getRefViscosity() const
{
    return refViscosity_;
}

void MeltableNormalSpecies::mix(MeltableNormalSpecies* const S, MeltableNormalSpecies* const T)
{
    elasticModulus_ = BaseSpecies::average(S->getElasticModulus(), T->getElasticModulus());
    poissonRatio_ = BaseSpecies::average(S->getPoissonRatio(), T->getPoissonRatio());
    ///\todo
}

Mdouble MeltableNormalSpecies::getEffectiveHeatCapacity(double temperature) const {
    if (temperature<meltingTemperature_-0.5*deltaT_) {
        return solidHeatCapacity_;
    } else if (temperature<meltingTemperature_+0.5*deltaT_) {
        return 0.5*(solidHeatCapacity_+liquidHeatCapacity_)+latentHeat_/deltaT_;
    } else {
        return liquidHeatCapacity_;
    }
}

Mdouble MeltableNormalSpecies::getEffectiveLatentHeat() const {
    return latentHeat_ + 0.5*(solidHeatCapacity_+liquidHeatCapacity_)*deltaT_;
}

Mdouble MeltableNormalSpecies::getRelativeSolidRadius(double temperature) const {
    return std::max(minRelativeSolidRadius_,std::min(1.0,std::cbrt(1-(temperature-meltingTemperature_+0.5*deltaT_)/deltaT_)));
}

void MeltableNormalSpecies::analyseTimeScales(double radius, double density, double temperature) const {
    using mathsFunc::cubic;
    using mathsFunc::square;
    using constants::pi;
    using constants::R;
    using constants::sqrt_2;
    using constants::stefanBoltzmanConstant;
    double mass = density*(4./3.*pi*cubic(radius));
    double te = sqrt(mass/(4./3.*elasticModulus_*radius));
    double td = mass/(2.*dissipation_*sqrt(2.*mass*elasticModulus_*radius));
    double ts = sqrt(mass/surfaceTension_);
    double viscosity = refViscosity_*exp(activationEnergy_/R/temperature);
    double tv = mass/(viscosity*sqrt_2*radius);
    double tcd = mass*solidHeatCapacity_/(thermalConductivityCoefficient_*0.5*pi*radius);
    double tcv = mass*solidHeatCapacity_/(thermalConvectionCoefficient_*pi*square(radius));
    double tr = mass*solidHeatCapacity_/(4*pi*square(radius)*materialEmissivity_*stefanBoltzmanConstant* cubic(temperature));
    logger(INFO, "Time scales:\n"
                 "- elastic %\n"
                 "- dissipative %\n"
                 "- surface tension %\n"
                 "- viscous %\n"
                 "- conduction %\n"
                 "- convection %\n"
                 "- radiation %\n"
                 "viscosity %", te, td, ts, tv, tcd, tcv, tr, viscosity);
}
