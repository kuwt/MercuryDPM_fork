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

#include "Helpers/FormulaHelpers.h"
#include "Logger.h"

/*!
 * \details The effective mass is an important parameter in a collision.
 * E.g. the collision time and the restitution coefficient are functions of the effective mass.
 * \param[in] mass0 The mass of the first particle.
 * \param[in] mass1 The mass of the second particle.
 * \return the effective mass of the particle pair.
 */
Mdouble helpers::getEffectiveMass(Mdouble mass0, Mdouble mass1)
{
    return mass0 * mass1 / (mass0 + mass1);
}

/*!
 * \brief Returns the Rayleigh time step for a Hertz contact law.
 * \details An accepted time step for Hertz is 10-20% of the Rayleigh time step.
 * See \cite Marigo2015
 */
Mdouble helpers::getRayleighTime(Mdouble radius, Mdouble shearModulus, Mdouble poisson, Mdouble density) 
{
    return constants::pi * radius * sqrt(density / shearModulus)/(0.1631 * poisson + 0.8766);
}

helpers::KAndDispAndKtAndDispt
helpers::computeDisptFromCollisionTimeAndRestitutionCoefficientAndTangentialRestitutionCoefficientAndEffectiveMass(
        Mdouble tc, Mdouble r, Mdouble beta, Mdouble mass)
{
    helpers::KAndDispAndKtAndDispt ans;
    ans.disp = -2.0 * mass * log(r) / tc;
    ans.k = mass * (mathsFunc::square(constants::pi / tc) + mathsFunc::square(ans.disp / (2.0 * mass)));
    ans.kt = 2.0 / 7.0 * ans.k * (mathsFunc::square(constants::pi) + mathsFunc::square(log(beta))) /
             (mathsFunc::square(constants::pi) + mathsFunc::square(log(r)));
    if (beta != 0.0)
        ans.dispt = -2 * log(beta) *
                    sqrt(1.0 / 7.0 * mass * ans.kt / (mathsFunc::square(constants::pi) + mathsFunc::square(log(beta))));
    else
        ans.dispt = 2. * sqrt(1.0 / 7.0 * mass * ans.kt);
    return ans;
}

Mdouble helpers::getMaximumVelocity(Mdouble k, Mdouble disp, Mdouble radius, Mdouble mass)
{
    // note: underestimate based on energy argument,
    // Ekin = 2*1/2*m*(v/2)^2 = 1/2*k*(2*r)^2, gives
    // return radius * sqrt(8.0*k/mass);

    // with dissipation, see S. Luding, Collisions & Contacts between two particles, eq 34
    Mdouble w = sqrt(k / mass - mathsFunc::square(disp / (2.0 * mass)));
    Mdouble w0 = sqrt(k / mass);
    Mdouble DispCorrection = exp(-disp / mass / w) * asin(w / w0);
    //std::cout << "DispCorrection" << DispCorrection << std::endl;
    return radius * sqrt(8.0 * k / mass) / DispCorrection;
}

/*!
 * \details MercuryDPM uses the DPMBase::setSaveCount to determine how often output is written.
 * However, if the user wants to set the total number of saves instead of the saveCount,
 * he can use this function to calculate the correct saveCount, assuming that the
 * final time and the mean time step is known.
 *
 * Example of use:
 * > setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(numberOfSaves, getTimeMax(), getTimeStep()));
 *
 * \param[in] numberOfSaves the total number of output files the user wants at the end of the simulation.
 * \param[in] timeMax       the final time of the simulation
 * \param[in] timeStep      the mean time step used during the simulation
 * \return the saveCount value that should be used to get the desired number of saves.
 */
unsigned int helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(
    unsigned int numberOfSaves,
    Mdouble timeMax,
    Mdouble timeStep)
{
    if (numberOfSaves > 0 && timeMax > 0 && timeStep > 0)
    {
        return static_cast<unsigned int>(ceil(
                (timeMax + timeStep) / timeStep / static_cast<double>(numberOfSaves - 1)));
    }
    else
    {
        logger(ERROR,
               "[Helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep()] numberOfSaves: %, timeMax: %, "
               "timestep: %\n Arguments need to be positive",
               numberOfSaves, timeMax, timeStep);
    }
    return 0;
}

/*!
 * \brief Returns the shear modulus calculated from Young's modulus and Poisson's ratio. 
 */
Mdouble helpers::getShearModulus(Mdouble youngsModulus, Mdouble poisson){
    return youngsModulus / (2.0 * (1 + poisson));
}
