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

/*!
 * \todo gmb Remove this file?
 */

#include "Helpers.h"

/*!
 * \todo gmb Remove these? If not move them to Helpers/FormulaHelpers.cc and use Logger()
 */
//helpers::KAndDisp helpers::computeKAndDispFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(Mdouble tc, Mdouble r, Mdouble mass)
//{
//    helpers::KAndDisp ans;
//    ans.k = computeKFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(tc, r, mass);
//    ans.disp = computeDispFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(tc, r, mass);
//    return ans;
//}
//
//Mdouble helpers::computeCollisionTimeFromKAndDispAndEffectiveMass(Mdouble k, Mdouble disp, Mdouble mass)
//{
//    if (k <= 0)
//    {
//        std::cerr << "Error in helpers::computeCollisionTimeFromKandDispAndEffectiveMass(Mdouble k, Mdouble disp, Mdouble mass) stiffness is not set or has an unexpected value, (stiffness=" << k << ")" << std::endl;
//        exit(-1);
//    }
//    if (disp < 0)
//    {
//        std::cerr << "Error in helpers::computeCollisionTimeFromKandDispAndEffectiveMass(Mdouble k, Mdouble disp, Mdouble mass) dissipation is not set or has an unexpected value, (dissipation=" << disp << ")" << std::endl;
//        exit(-1);
//    }
//    if (mass <= 0)
//    {
//        std::cerr << "Error in helpers::computeCollisionTimeFromKandDispAndEffectiveMass(Mdouble k, Mdouble disp, Mdouble mass) mass is not set or has an unexpected value (mass=" << mass << ")" << std::endl;
//        exit(-1);
//    }
//    if (4.0 * k / mass - mathsFunc::square(disp / mass) < 0)
//    {
//        std::cerr << "Error in helpers::computeCollisionTimeFromKandDispAndEffectiveMass(Mdouble k, Mdouble disp, Mdouble mass) values for stiffness, dissipation and mass lead to an overdamped system (stiffness=" << k << " dissipation=" << disp << " mass=" << mass << ")" << std::endl;
//        exit(-1);
//    }
//    return 2.0 * constants::pi / std::sqrt(4.0 * k / mass - mathsFunc::square(disp / mass));
//}
//
//Mdouble helpers::computeRestitutionCoefficientFromKAndDispAndEffectiveMass(Mdouble k, Mdouble disp, Mdouble mass)
//{
//    if (k <= 0)
//    {
//        std::cerr << "Error in helpers::computeRestitutionCoefficientFromKandDispAndEffectiveMass(Mdouble k, Mdouble disp, Mdouble mass) stiffness is not set or has an unexpected value, (stiffness=" << k << ")" << std::endl;
//        exit(-1);
//    }
//    if (disp < 0)
//    {
//        std::cerr << "Error in helpers::computeRestitutionCoefficientFromKandDispAndEffectiveMass(Mdouble k, Mdouble disp, Mdouble mass) dissipation is not set or has an unexpected value, (dissipation=" << disp << ")" << std::endl;
//        exit(-1);
//    }
//    if (mass <= 0)
//    {
//        std::cerr << "Error in helpers::computeRestitutionCoefficientFromKandDispAndEffectiveMass(Mdouble k, Mdouble disp, Mdouble mass) mass is not set or has an unexpected value (mass=" << mass << ")" << std::endl;
//        exit(-1);
//    }
//    if (4.0 * mass * k - mathsFunc::square(disp) < 0)
//    {
//        std::cerr << "Error in helpers::computeRestitutionCoefficientFromKandDispAndEffectiveMass(Mdouble k, Mdouble disp, Mdouble mass) values for stiffness, dissipation and mass lead to an overdamped system (stiffness=" << k << " dissipation=" << disp << " mass=" << mass << ")" << std::endl;
//        exit(-1);
//    }
//    return std::exp(-disp * constants::pi / std::sqrt(4.0 * mass * k - mathsFunc::square(disp)));
//}
//
//Mdouble helpers::computeDispFromKAndRestitutionCoefficientAndEffectiveMass(Mdouble k, Mdouble r, Mdouble mass)
//{
//    if (k <= 0)
//    {
//        std::cerr << "Error in helpers::computeDispFromKAndRestitutionCoefficientAndEffectiveMass(Mdouble k, Mdouble r, Mdouble mass) stiffness is not set or has an unexpected value, (stiffness=" << k << ")" << std::endl;
//        exit(-1);
//    }
//    if (r < 0)
//    {
//        std::cerr << "Error in helpers::computeDispFromKAndRestitutionCoefficientAndEffectiveMass(Mdouble k, Mdouble r, Mdouble mass) restitution coefficient is not set or has an unexpected value, (restitution coefficient=" << r << ")" << std::endl;
//        exit(-1);
//    }
//    if (mass <= 0)
//    {
//        std::cerr << "Error in helpers::computeDispFromKAndRestitutionCoefficientAndEffectiveMass(Mdouble k, Mdouble r, Mdouble mass) mass is not set or has an unexpected value (mass=" << mass << ")" << std::endl;
//        exit(-1);
//    }
//    return -2.0 * sqrt(mass * k / (constants::sqr_pi + mathsFunc::square(std::log(r)))) * std::log(r);
//}
//
//Mdouble helpers::computeCollisionTimeFromKAndRestitutionCoefficientAndEffectiveMass(Mdouble k, Mdouble r, Mdouble mass)
//{
//    if (k <= 0)
//    {
//        std::cerr << "Error in helpers::computeCollisionTimeFromKAndRestitutionCoefficientAndEffectiveMass(Mdouble k, Mdouble r, Mdouble mass) stiffness is not set or has an unexpected value, (stiffness=" << k << ")" << std::endl;
//        exit(-1);
//    }
//    if (r < 0)
//    {
//        std::cerr << "Error in helpers::computeCollisionTimeFromKAndRestitutionCoefficientAndEffectiveMass(Mdouble k, Mdouble r, Mdouble mass) restitution coefficient is not set or has an unexpected value, (restitution coefficient=" << r << ")" << std::endl;
//        exit(-1);
//    }
//    if (mass <= 0)
//    {
//        std::cerr << "Error in helpers::computeCollisionTimeFromKAndRestitutionCoefficientAndEffectiveMass(Mdouble k, Mdouble r, Mdouble mass) mass is not set or has an unexpected value (mass=" << mass << ")" << std::endl;
//        exit(-1);
//    }
//    return sqrt(mass / k * (constants::sqr_pi + mathsFunc::square(std::log(r))));
//}
//
//Mdouble helpers::computeDispFromKAndCollisionTimeAndEffectiveMass(Mdouble k, Mdouble tc, Mdouble mass)
//{
//    if (k <= 0)
//    {
//        std::cerr << "Error in helpers::computeDispFromKAndCollisionTimeAndEffectiveMass(Mdouble k, Mdouble tc, Mdouble mass) stiffness is not set or has an unexpected value, (stiffness=" << k << ")" << std::endl;
//        exit(-1);
//    }
//    if (tc <= 0)
//    {
//        std::cerr << "Error in helpers::computeDispFromKAndCollisionTimeAndEffectiveMass(Mdouble k, Mdouble tc, Mdouble mass) collision time is not set or has an unexpected value, (collision time=" << tc << ")" << std::endl;
//        exit(-1);
//    }
//    if (mass <= 0)
//    {
//        std::cerr << "Error in helpers::computeDispFromKAndCollisionTimeAndEffectiveMass(Mdouble k, Mdouble tc, Mdouble mass) mass is not set or has an unexpected value (mass=" << mass << ")" << std::endl;
//        exit(-1);
//    }
//    if (mass * k - constants::sqr_pi * mathsFunc::square(mass / tc) < 0)
//    {
//        std::cerr << "Error in helpers::computeDispFromKAndCollisionTimeAndEffectiveMass(Mdouble k, Mdouble disp, Mdouble mass) values for stiffness, collision time and mass lead to an overdamped system (stiffness=" << k << " collision time=" << tc << " mass=" << mass << ")" << std::endl;
//        exit(-1);
//    }
//    return 2.0 * std::sqrt(mass * k - constants::sqr_pi * mathsFunc::square(mass / tc));
//}
//
//Mdouble helpers::computeRestitutionCoefficientFromKAndCollisionTimeAndEffectiveMass(Mdouble k, Mdouble tc, Mdouble mass)
//{
//    if (k <= 0)
//    {
//        std::cerr << "Error in helpers::computeRestitutionCoefficientFromKAndCollisionTimeAndEffectiveMass(Mdouble k, Mdouble tc, Mdouble mass) stiffness is not set or has an unexpected value, (stiffness=" << k << ")" << std::endl;
//        exit(-1);
//    }
//    if (tc <= 0)
//    {
//        std::cerr << "Error in helpers::computeRestitutionCoefficientFromKAndCollisionTimeAndEffectiveMass(Mdouble k, Mdouble tc, Mdouble mass) collision time is not set or has an unexpected value, (collision time=" << tc << ")" << std::endl;
//        exit(-1);
//    }
//    if (mass <= 0)
//    {
//        std::cerr << "Error in helpers::computeRestitutionCoefficientFromKAndCollisionTimeAndEffectiveMass(Mdouble k, Mdouble tc, Mdouble mass) mass is not set or has an unexpected value (mass=" << mass << ")" << std::endl;
//        exit(-1);
//    }
//    if (k / mass * mathsFunc::square(tc) - constants::sqr_pi < 0)
//    {
//        std::cerr << "Error in helpers::computeRestitutionCoefficientFromKAndCollisionTimeAndEffectiveMass(Mdouble k, Mdouble disp, Mdouble mass) values for stiffness, collision time and mass lead to an overdamped system (stiffness=" << k << " collision time=" << tc << " mass=" << mass << ")" << std::endl;
//        exit(-1);
//    }
//    return std::exp(-std::sqrt(k / mass * mathsFunc::square(tc) - constants::sqr_pi));
//}
//
//Mdouble helpers::computeDispFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(Mdouble tc, Mdouble r, Mdouble mass)
//{
//    if (tc <= 0)
//    {
//        std::cerr << "Error in helpers::computeDispFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(Mdouble tc, Mdouble r, Mdouble mass) collision time is not set or has an unexpected value, (collision time=" << tc << ")" << std::endl;
//        exit(-1);
//    }
//    if (r < 0)
//    {
//        std::cerr << "Error in helpers::computeDispFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(Mdouble tc, Mdouble r, Mdouble mass) restitution coefficient not set or has an unexpected value, (restitution coefficient=" << r << ")" << std::endl;
//        exit(-1);
//    }
//    if (mass <= 0)
//    {
//        std::cerr << "Error in helpers::computeDispFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(Mdouble tc, Mdouble r, Mdouble mass) mass is not set or has an unexpected value (mass=" << mass << ")" << std::endl;
//        exit(-1);
//    }
//    return -2.0 * mass * std::log(r) / tc;
//}
//
//Mdouble helpers::computeKFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(Mdouble tc, Mdouble r, Mdouble mass)
//{
//    if (tc <= 0)
//    {
//        std::cerr << "Error in helpers::computeKFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(Mdouble tc, Mdouble r, Mdouble mass) collision time is not set or has an unexpected value, (collision time=" << tc << ")" << std::endl;
//        exit(-1);
//    }
//    if (r < 0)
//    {
//        std::cerr << "Error in helpers::computeKFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(Mdouble tc, Mdouble r, Mdouble mass) restitution coefficient not set or has an unexpected value, (restitution coefficient=" << r << ")" << std::endl;
//        exit(-1);
//    }
//    if (mass <= 0)
//    {
//        std::cerr << "Error in helpers::computeKFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(Mdouble tc, Mdouble r, Mdouble mass) mass is not set or has an unexpected value (mass=" << mass << ")" << std::endl;
//        exit(-1);
//    }
//    return mass * (mathsFunc::square(constants::pi / tc) + mathsFunc::square(-std::log(r) / tc));
//}
//
//Mdouble helpers::computeKFromCollisionTimeAndDispAndEffectiveMass(Mdouble tc, Mdouble disp, Mdouble mass)
//{
//    if (tc <= 0)
//    {
//        std::cerr << "Error in helpers::computeKFromCollisionTimeAndDispAndEffectiveMass(Mdouble tc, Mdouble disp, Mdouble mass) collision time is not set or has an unexpected value, (collision time=" << tc << ")" << std::endl;
//        exit(-1);
//    }
//    if (disp < 0)
//    {
//        std::cerr << "Error in helpers::computeKFromCollisionTimeAndDispAndEffectiveMass(Mdouble tc, Mdouble disp, Mdouble mass) dissipation not set or has an unexpected value, (dissipation=" << disp << ")" << std::endl;
//        exit(-1);
//    }
//    if (mass <= 0)
//    {
//        std::cerr << "Error in helpers::computeKFromCollisionTimeAndDispAndEffectiveMass(Mdouble tc, Mdouble disp, Mdouble mass) mass is not set or has an unexpected value (mass=" << mass << ")" << std::endl;
//        exit(-1);
//    }
//    return 0.25 * mathsFunc::square(disp) / mass + constants::sqr_pi * mass / mathsFunc::square(tc);
//}
//
//Mdouble helpers::computeRestitutionCoefficientFromCollisionTimeAndDispAndEffectiveMass(Mdouble tc, Mdouble disp, Mdouble mass)
//{
//    if (tc <= 0)
//    {
//        std::cerr << "Error in helpers::computeRestitutionCoefficientFromCollisionTimeAndDispAndEffectiveMass(Mdouble tc, Mdouble disp, Mdouble mass) collision time is not set or has an unexpected value, (collision time=" << tc << ")" << std::endl;
//        exit(-1);
//    }
//    if (disp < 0)
//    {
//        std::cerr << "Error in helpers::computeRestitutionCoefficientFromCollisionTimeAndDispAndEffectiveMass(Mdouble tc, Mdouble disp, Mdouble mass) dissipation not set or has an unexpected value, (dissipation=" << disp << ")" << std::endl;
//        exit(-1);
//    }
//    if (mass <= 0)
//    {
//        std::cerr << "Error in helpers::computeRestitutionCoefficientFromCollisionTimeAndDispAndEffectiveMass(Mdouble tc, Mdouble disp, Mdouble mass) mass is not set or has an unexpected value (mass=" << mass << ")" << std::endl;
//        exit(-1);
//    }
//    return std::exp(-0.5 * disp * tc / mass);
//}
//
//Mdouble helpers::computeKFromDispAndRestitutionCoefficientAndEffectiveMass(Mdouble disp, Mdouble r, Mdouble mass)
//{
//    if (disp <= 0)
//    {
//        std::cerr << "Error in helpers::computeKFromDispAndRestitutionCoefficientAndEffectiveMass(Mdouble disp, Mdouble r, Mdouble mass) dissipation is not set or has an unexpected value, (dissipation time=" << disp << ")" << std::endl;
//        exit(-1);
//    }
//    if (r < 0)
//    {
//        std::cerr << "Error in helpers::computeKFromDispAndRestitutionCoefficientAndEffectiveMass(Mdouble disp, Mdouble r, Mdouble mass) restitution coefficient not set or has an unexpected value, (restitution coefficient=" << r << ")" << std::endl;
//        exit(-1);
//    }
//    if (mass <= 0)
//    {
//        std::cerr << "Error in helpers::computeKFromDispAndRestitutionCoefficientAndEffectiveMass(Mdouble disp, Mdouble r, Mdouble mass) mass is not set or has an unexpected value (mass=" << mass << ")" << std::endl;
//        exit(-1);
//    }
//    return 0.25 * mathsFunc::square(disp)*(constants::sqr_pi / (mathsFunc::square(std::log(r))) + 1.0) / mass;
//}
//
//Mdouble helpers::computeCollisionTimeFromDispAndRestitutionCoefficientAndEffectiveMass(Mdouble disp, Mdouble r, Mdouble mass)
//{
//    if (disp <= 0)
//    {
//        std::cerr << "Error in helpers::computeCollisionTimeFromDispAndRestitutionCoefficientAndEffectiveMass(Mdouble disp, Mdouble r, Mdouble mass) dissipation is not set or has an unexpected value, (dissipation time=" << disp << ")" << std::endl;
//        exit(-1);
//    }
//    if (r < 0)
//    {
//        std::cerr << "Error in helpers::computeCollisionTimeFromDispAndRestitutionCoefficientAndEffectiveMass(Mdouble disp, Mdouble r, Mdouble mass) restitution coefficient not set or has an unexpected value, (restitution coefficient=" << r << ")" << std::endl;
//        exit(-1);
//    }
//    if (mass <= 0)
//    {
//        std::cerr << "Error in helpers::computeCollisionTimeFromDispAndRestitutionCoefficientAndEffectiveMass(Mdouble disp, Mdouble r, Mdouble mass) mass is not set or has an unexpected value (mass=" << mass << ")" << std::endl;
//        exit(-1);
//    }
//    return -2.0 * mass * std::log(r) / disp;
//}
/////from Deen...Kuipers2006, eq. 43 and 30

//seems to be unused, consider taking out \author weinhartt
//unsigned int helpers::getSaveCountFromNumberOfSavesPerTimeUnitAndTimeStep(unsigned int numberOfSaves, Mdouble time step)
//{
//    if (numberOfSaves > 1 && time step > 0)
//    {
//        return static_cast<unsigned int>(ceil(1.0 / time step / static_cast<double>(numberOfSaves - 1)));
//    }
//    else
//    {
//        std::cerr << "Error in getSaveCountFromNumberOfSavesPerTimeUnitAndTimeStep (" << numberOfSaves << "," << time step << ")" << std::endl;
//        std::cerr << "Arguments need to be positive" << std::endl;
//        exit(-1);
//    }
//}
