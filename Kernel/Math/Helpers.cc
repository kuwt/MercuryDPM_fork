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

#include "Helpers.h"
#include <cmath>
#include <sys/stat.h>
#include <cstdio>
#include <sstream>
#include <string>
#include <Logger.h>
#include <DPMBase.h>
#include <Walls/InfiniteWall.h>
#include <Particles/BaseParticle.h>
#include <Species/BaseSpecies.h>
//#include <cassert>
#include "Math/ExtendedMath.h"
#include <numeric>
#include <chrono>
#include <sys/types.h> // required for stat.h

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else

#include <unistd.h>
#include "Particles/SphericalParticle.h"

#define GetCurrentDir getcwd
#endif


std::string helpers::lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

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
unsigned int helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(unsigned int numberOfSaves, Mdouble timeMax,
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
               "[Helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep()] numberOfSaves: %, timeMax: %, timestep: %",
               numberOfSaves, timeMax, timeStep);
        logger(ERROR, " Arguments need to be positive");
        exit(-1);
    }
}

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

/*!
 * \details This function is used to avoid errors from reading in old or manually modified restart files.
 * Instead of reading variable by variable directly from the restart stringstream,
 * a full line is read first, from which the variables are read. Thus, if a line
 * has the wrong number of arguments, it might affect the reading of the current
 * line, but correctly reads the next line.
 *
 * Example of usage:
 * > std::stringstream line;
 * > std::stringstream is = restartFile.getFStream();
 * > helpers::getLineFromStringStream(is, line);
 * > std::string dummy;
 * > line >> dummy;
 *
 * \param[in]  in the stringstream from which a line is read out should be initialized as std::stringstream(std::stringstream::out)
 * \param[out] out the stringstream into which the line is read; should be initialized as std::stringstream(std::stringstream::in | std::stringstream::out)
 */
void helpers::getLineFromStringStream(std::istream& in, std::stringstream& out)
{
    std::string line_string;
    getline(in, line_string);
    out.str(std::move(line_string));
    out.clear();
}

/*!
 * \details Provides a simple interface for writing a string to a file.
 * This function is mainly used to create ini or restart file that the code
 * later reads back in.
 *
 * Example of usage:
 * > helpers::writeToFile("RestartUnitTest.ini",
 * > "1 0 0 0 0 1 1 1\n"
 * > "0.5 0.5 0  0 0 0.5  0  0 0 0  0 0 0  0\n");
 *
 * \param[in] filename the name of the file
 * \param[in] filecontent the content
 * \returns true on success.
 */
bool helpers::writeToFile(std::string filename, std::string filecontent)
{
    std::fstream file;
    file.open(filename.c_str(), std::ios::out);
    if (file.fail())
    {
        logger(WARN, "Error in writeToFile: file could not be opened");
        return false;
    }
    file << filecontent;
    file.close();
    return true;
}

void helpers::writeCommandLineToFile(const std::string filename, const int argc, char * const argv[])
{
    std::stringstream ss;
    for (int i=0; i<argc; ++i) {
        ss << argv[i] << ' ';
    }
    writeToFile(filename,ss.str());
}


void helpers::gnuplot(std::string command)
{
#ifdef __CYGWIN__
    logger(WARN, "[helpers::gnuplot] is not supported on Cygwin");
#else
#ifdef WINDOWS
    logger(WARN, "[helpers::gnuplot] is not supported on Windows");
#else
    FILE* pipe = popen("gnuplot -persist", "w");
    fprintf(pipe, "%s", command.c_str());
    fflush(pipe);
#endif
#endif
}

bool helpers::addToFile(std::string filename, std::string filecontent)
{
    std::fstream file;
    file.open(filename.c_str(), std::ios::app);
    if (file.fail())
    {
        std::cerr << "Error in writeToFile: file could not be opened" << std::endl;
        return false;
    }
    file << filecontent;
    file.close();
    return true;
}

/*!
 * \details This is a FileExist routine, which is used to test if a run have
 * already need preformed, allows me to plug holes in parm studies.
 */
bool helpers::fileExists(std::string strFilename)
{
    struct stat stFileInfo;
    bool blnReturn;
    int intStat;

    // Attempt to get the file attributes

    intStat = stat(strFilename.c_str(), &stFileInfo);
    if (intStat == 0)
    {
        // We were able to get the file attributes
        // so the file obviously exists.
        blnReturn = true;
    }
    else
    {
        // We were not able to get the file attributes.
        // This may mean that we don't have permission to
        // access the folder which contains this file. If you
        // need to do that level of checking, lookup the
        // return values of stat which will give you
        // more details on why stat failed.
        blnReturn = false;
    }

    return blnReturn;
}

/*!
 * \details Provides a simple interface for opening a file, in order to avoid
 * that the user has to learn the syntax for opening a file.
 * \param[out] file The std::fstream object that the user can write to.
 * \param[in] filename The name of the file.
 * \param[in] mode The openmode of the file, typically std::fstream::out or std::fstream::in.
 * \return true is the file was successfully opened, false else.
 */
bool helpers::openFile(std::fstream& file, std::string filename, std::fstream::openmode mode)
{
    file.open(filename.c_str(), mode);
    if (file.fail())
        return false;
    else
        return true;
}

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

std::vector<double> helpers::readArrayFromFile(std::string filename, int& n, int& m)
{
    std::fstream file;
    file.open(filename.c_str(), std::ios::in);
    if (file.fail())
    {
        std::cerr << "Error in readArrayFromFile: file could not be opened" << std::endl;
        exit(-1);
    }
    file >> n >> m;
    std::vector<double> v;
    Mdouble val;
    for (int i = 0; i < n * m; i++)
    {
        file >> val;
        v.push_back(val);
    }
    file.close();
    return v;
}

void helpers::more(std::string filename, unsigned nLines)
{
    if (nLines != constants::unsignedMax)
        std::cout << "First " << nLines << " lines of " << filename << ":\n";
    std::fstream file;
    file.open(filename.c_str(), std::ios::in);
    if (file.fail())
        logger(ERROR, "Error in readArrayFromFile: file could not be opened");
    std::string line;
    for (unsigned i = 0; i < nLines; i++)
    {
        if (file.eof()) break;
        getline(file, line);
        std::cout << " " << line << '\n';
    }
    file.close();
}

Mdouble helpers::round(const Mdouble value, unsigned precision)
{
    const Mdouble logValue = log10(value);
    const int factor = std::pow(10, precision - std::ceil(logValue));
    return std::round(value * factor) / factor;
}

std::string helpers::to_string(const Mdouble value, unsigned precision)
{
    std::ostringstream stm;
    stm << round(value,precision);
    return stm.str();
}

/**
 * Creates a DPMBase with a particles of unit size and a flat wall and loads/unloads the particle-wall contact
 * \param[in] species      particle species specifying the contact law
 * \param[in] displacement peak displacement before unloading
 * \param[in] velocity     loading/unloading velocity
 */
void helpers::loadingTest(const ParticleSpecies* species, Mdouble displacement, Mdouble velocity, Mdouble radius,
                          std::string name)
{
    class LoadingTest : public DPMBase
    {
        const ParticleSpecies* species;
        Mdouble displacement;
        Mdouble velocity;
        Mdouble radius;
    public:
        //public variables
        LoadingTest(const ParticleSpecies* species, Mdouble displacement, Mdouble velocity, Mdouble radius)
                : species(species), displacement(displacement), velocity(velocity), radius(radius)
        {}

        void setupInitialConditions() override
        {
            //setName("LoadingTest"+species->getName());
            setTimeMax(2.0 * displacement / velocity);
            setTimeStep(2e-3 * getTimeMax());
            setSaveCount(1);
            setFileType(FileType::NO_FILE);
            fStatFile.setFileType(FileType::ONE_FILE);

            setMax({radius, radius, radius + radius});
            setMin({-radius, -radius, 0});
            setSystemDimensions(3);
            setParticleDimensions(3);

            speciesHandler.copyAndAddObject(*species);

            SphericalParticle p;
            p.setSpecies(speciesHandler.getObject(0));
            p.setRadius(radius);
            p.setPosition({0, 0, radius});
            particleHandler.copyAndAddObject(p);

            InfiniteWall w;
            w.setSpecies(speciesHandler.getObject(0));
            w.set(Vec3D(0, 0, -1), Vec3D(0.0, 0.0, 0.0));
            wallHandler.copyAndAddObject(w);
        }

        void actionsBeforeTimeStep() override
        {
            BaseParticle* p = particleHandler.getLastObject();
            logger.assert(p,"Empty particle handler");
            p->setAngularVelocity({0, 0, 0});

            //Moving particle normally into surface
            if (getTime() <= displacement / velocity)
            {
                p->setVelocity({0, 0, velocity});
                p->setPosition({0, 0, radius - velocity * getTime()});
            }
            else
            {
                p->setVelocity({0, 0, -velocity});
                p->setPosition({0, 0, radius - displacement - displacement + velocity * getTime()});
            }
        }
    } test(species, displacement, velocity, radius);
    test.setName(name);
    test.solve();
    writeToFile(test.getName() + ".gnu", "plot '" + test.getName() + ".fstat' u 7:9 w lp");
    logger(INFO, "finished loading test: run 'gnuplot %.gnu' to view output", test.getName());
}

/**
 * Creates a DPMBase with a particles of unit size and a flat wall and loads/unloads/reloads the particle-wall contact in tangential direction
 * \param[in] species      particle species specifying the contact law
 * \param[in] displacement peak displacement before unloading
 * \param[in] velocity     loading/unloading velocity
 */
void helpers::normalAndTangentialLoadingTest(const ParticleSpecies* species, Mdouble displacement,
                                             Mdouble tangentialDisplacement, Mdouble velocity, Mdouble radius,
                                             std::string name)
{
    class LoadingTest : public DPMBase
    {
        const ParticleSpecies* species;
        Mdouble displacement;
        Mdouble tangentialDisplacement;
        Mdouble velocity;
        Mdouble radius;
    public:
        //public variables
        LoadingTest(const ParticleSpecies* species, Mdouble displacement, Mdouble tangentialDisplacement,
                    Mdouble velocity, Mdouble radius)
                : species(species), displacement(displacement), tangentialDisplacement(tangentialDisplacement),
                  velocity(velocity), radius(radius)
        {}

        void setupInitialConditions() override
        {
            //setName("TangentialLoadingTest"+species->getName());
            setTimeMax(4.0 * tangentialDisplacement / velocity);
            setTimeStep(4e-4 * getTimeMax());
            setSaveCount(1);
            setFileType(FileType::NO_FILE);
            fStatFile.setFileType(FileType::ONE_FILE);

            setMax({radius, radius, radius + radius});
            setMin({-radius, -radius, 0});
            setSystemDimensions(3);
            setParticleDimensions(3);

            speciesHandler.copyAndAddObject(*species);

            SphericalParticle p;
            p.setSpecies(speciesHandler.getObject(0));
            p.setRadius(radius);
            p.setPosition({0, 0, radius - displacement});
            particleHandler.copyAndAddObject(p);

            InfiniteWall w;
            w.setSpecies(speciesHandler.getObject(0));
            w.set(Vec3D(0, 0, -1), Vec3D(0.0, 0.0, 0.0));
            wallHandler.copyAndAddObject(w);
        }

        void actionsBeforeTimeStep() override
        {
            BaseParticle* p = particleHandler.getLastObject();
            logger.assert(p,"Empty particle handler");
            p->setAngularVelocity({0, 0, 0});
    
            //Moving particle cyclically right and left between +-tangentialDisplacement
            bool moveRight = static_cast<int>(getTime() / (2.0*tangentialDisplacement / velocity) +0.5)%2==0;
            if (moveRight)
            {
                p->setVelocity({-velocity, 0, 0});
                p->setPosition({tangentialDisplacement - velocity * getTime(), 0, radius - displacement});
            }
            else
            {
                p->setVelocity({velocity, 0, 0});
                p->setPosition({-2*tangentialDisplacement + velocity * getTime(), 0, radius - displacement});
            }
        }

    } test(species, displacement, tangentialDisplacement, velocity, radius);
    test.setName(name);
    test.solve();
    writeToFile(test.getName() + ".gnu", "plot '" + test.getName() + ".fstat' u 8:($10*$14) w lp");
    logger(INFO, "finished tangential loading test: run 'gnuplot %.gnu' to view output", test.getName());
}

/**
 * Creates a DPMBase with a particles of unit size and a flat wall, loads the particle-wall contact in normal and tangential direction, then rotates.
 * \param[in] species      particle species specifying the contact law
 * \param[in] displacement peak displacement before unloading
 * \param[in] velocity     loading/unloading velocity
 */
void helpers::objectivenessTest(const ParticleSpecies* species, Mdouble displacement, Mdouble tangentialDisplacement,
                                Mdouble velocity, Mdouble radius, std::string name)
{
    class ObjectivenessTest : public DPMBase
    {
        const ParticleSpecies* species;
        Mdouble displacement;
        Mdouble tangentialDisplacement;
        Mdouble velocity;
        Mdouble radius;
    public:
        //public variables
        ObjectivenessTest(const ParticleSpecies* species, Mdouble displacement, Mdouble tangentialDisplacement,
                          Mdouble velocity, Mdouble radius)
                : species(species), displacement(displacement), tangentialDisplacement(tangentialDisplacement),
                  velocity(velocity), radius(radius)
        {}

        void setupInitialConditions() override
        {
            //setName("ObjectivenessTest"+species->getName());
            setTimeMax((tangentialDisplacement + 0.5 * constants::pi * radius) / velocity);
            setTimeStep(1e-4 * getTimeMax());
            setSaveCount(20);
            setFileType(FileType::NO_FILE);
            dataFile.setFileType(FileType::ONE_FILE);
            fStatFile.setFileType(FileType::ONE_FILE);

            setMax(radius * Vec3D(2, 2, 2));
            setMin(radius * Vec3D(-2, -2, -2));
            setSystemDimensions(2);
            setParticleDimensions(3);

            speciesHandler.copyAndAddObject(*species);
    
            SphericalParticle p;
            p.setSpecies(speciesHandler.getObject(0));
            p.setRadius(radius);
            p.setPosition({0, radius - displacement, 0});
            particleHandler.copyAndAddObject(p);
            p.setPosition({0, -radius + displacement, 0});
            particleHandler.copyAndAddObject(p);
        }

        void actionsBeforeTimeStep() override
        {
            BaseParticle* p = particleHandler.getObject(0);
            BaseParticle* q = particleHandler.getLastObject();
            logger.assert(p,"Empty particle handler");
            logger.assert(q,"Empty particle handler");

            //Moving particle normally into surface
            if (getTime() <= tangentialDisplacement / velocity)
            {
                p->setAngularVelocity({0, 0, 0});
                p->setVelocity({velocity, 0, 0});
                p->setPosition({-tangentialDisplacement + velocity * getTime(), radius - displacement, 0});
                q->setAngularVelocity({0, 0, 0});
                q->setVelocity({-velocity, 0, 0});
                q->setPosition({tangentialDisplacement - velocity * getTime(), -radius + displacement, 0});
            }
            else
            {
                Mdouble angle = velocity / (radius - displacement) * (getTime() - tangentialDisplacement / velocity);
                Mdouble s = sin(angle);
                Mdouble c = cos(angle);
                p->setAngularVelocity(velocity / (radius - displacement) * Vec3D(0, 0, -1));
                //p->setAngularVelocity(Vec3D(0,0,0));
                p->setOrientation({1, 0, 0, 0});
                p->setVelocity(velocity * Vec3D(c, -s, 0));
                //p->setVelocity(Vec3D(0,0,0));
                p->setPosition((radius - displacement) * Vec3D(s, c, 0));
                q->setAngularVelocity(-p->getAngularVelocity());
                q->setOrientation(-p->getOrientation());
                q->setVelocity(-p->getVelocity());
                q->setPosition(-p->getPosition());
            }
        }

    } test(species, displacement, tangentialDisplacement, velocity, radius);
    test.setName(name);
    test.solve();
    writeToFile(test.getName() + ".gnu", "set size ratio -1; plot '" + test.getName() + ".fstat' u 14:15 every 2 w lp");
    logger(INFO, "finished objectiveness test: run 'gnuplot %.gnu' to view output", test.getName());
}

bool helpers::compare(std::istream& is, std::string s)
{
    // Get current position
    //check if the next line starts with 'interactionFile'; otherwise, skip interaction
    int len = is.tellg();
    std::string dummy;
    is >> dummy;
    if (dummy != s)
    {
        is.seekg(len, std::ios_base::beg);
        return false;
        logger(INFO, "helpers::compare: Next stream value (%) is not %", dummy, s);
    }
    return true;
}

bool helpers::readFromCommandLine(int argc, char *argv[], std::string varName)
{
    for (unsigned i=0; i<argc; ++i) {
        if (varName == argv[i]) {
            logger(INFO, "readFromCommandLine: % set to true",varName.substr(1));
            return true;
        }
    }
    //if the variable is not found
    logger(INFO, "readFromCommandLine: % set to default value false",varName.substr(1));
    return false;
}


template<class T>
void checkTemplate(T real, T ideal, double error, std::string whatIsChecked)
{
    logger.assert_always(mathsFunc::isEqual(real, ideal, error),
                         whatIsChecked + ": % (should be %) ", real, ideal);
    logger(INFO, whatIsChecked + ": % (correct)", real);
}

void helpers::check(double real, double ideal, double error, std::string whatIsChecked)
{
    checkTemplate(real,ideal,error,whatIsChecked);
}

void helpers::check(Vec3D real, Vec3D ideal, double error, std::string whatIsChecked)
{
    checkTemplate(real,ideal,error,whatIsChecked);
}

void helpers::check(MatrixSymmetric3D real, MatrixSymmetric3D ideal, double error, std::string whatIsChecked)
{
    checkTemplate(real,ideal,error,whatIsChecked);
}

void helpers::check(Matrix3D real, Matrix3D ideal, double error, std::string whatIsChecked)
{
    checkTemplate(real,ideal,error,whatIsChecked);
}

std::string helpers::getPath()
{
    char cCurrentPath[FILENAME_MAX];
    if (not GetCurrentDir(cCurrentPath, sizeof(cCurrentPath))) {
        logger(WARN, "Get current dir failed: %", cCurrentPath);
    }
    return std::string(cCurrentPath);
}

Mdouble helpers::getRealTime()
{
    // record start time
    static auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end - start;
    return diff.count();
}

/**
 * reads next value in stream as a string and compares it with name.
 * If name is equal to string, the function outputs true.
 * If name is not equal to string, the function undoes the read by setting seekg, and outputs false.
 * @param is
 * @param name
 * @return
 */
bool helpers::isNext(std::istream& is, const std::string name) {
    std::string dummy;
    auto pos = is.tellg();
    is >> dummy;
    if (dummy != name) {
        is.seekg(pos);
        return false;
    } else {
        return true;
    }
}

template<>
std::string helpers::readFromCommandLine<std::string>(int argc, char *argv[], std::string varName, std::string value)
{
    for (unsigned i=0; i<argc-1; ++i) {
        if (varName == argv[i]) {
            value = argv[i+1];
            logger(INFO, "readFromCommandLine: % set to % ", varName.substr(1), value);
            return value;
        }
    }
    //if the variable is not found
    logger(INFO, "readFromCommandLine: % set to default value % ", varName.substr(1), value);
    return value;
}

bool helpers::createDirectory(std::string path) {
    //see https://stackoverflow.com/questions/20358455/cross-platform-way-to-make-a-directory
    mode_t nMode = 0733; // UNIX style permissions
    int nError = 0;
#if defined(_WIN32)
    nError = _mkdir(path.c_str()); // can be used on Windows
#else
    nError = mkdir(path.c_str(),nMode); // can be used on non-Windows
#endif
    if (nError != 0) {
        // handle your error here
    }
    return false;
}

Mdouble helpers::getRayleighTime(Mdouble radius, Mdouble shearModulus, Mdouble poisson, Mdouble density) {
    return constants::pi*radius*sqrt(density/shearModulus)/(0.1631*poisson+0.8766);
}