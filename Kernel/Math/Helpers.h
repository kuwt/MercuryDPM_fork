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

#ifndef HELPERS_H
#define HELPERS_H

#include "ExtendedMath.h"
#include <iostream> //std::istream and std::stringstream
#include <fstream> //std::fstream
#include <vector> //std::fstream
#include <Logger.h>

class ParticleSpecies;

namespace helpers
{
    // returns the input string after converting upper-case characters to lower case
    std::string lower(std::string s);

    /*!
 * \brief return type specifically for fuctions returning k and disp at once
 */
class KAndDisp
{
public:
    Mdouble k;
    Mdouble disp;
};

//    /*!
//     * \brief Set disp and k such that is matches a given collision time tc and restitution coefficient r
//     * for a collision of effective/reduced mass m.
//     * \deprecated use species->setCollisionTimeAndRestitutionCoefficient
//     *   (collisionTime, dissipationTimeScale, 2.0*effectiveMass) instead
//     */
//    MERCURY_DEPRECATED
//    KAndDisp computeKAndDispFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(Mdouble tc, Mdouble r, Mdouble mass);
//
//    /*!
//     * \brief Calculates the collision time for a given stiffness, dissipation, and effective mass
//     * \deprecated use species->computeCollisionTime(2.0*effectiveMass) instead
//     * \todo This does not result in the same value as the given alternative.
//     */
//    MERCURY_DEPRECATED
//    Mdouble computeCollisionTimeFromKAndDispAndEffectiveMass(Mdouble k, Mdouble disp, Mdouble mass);
//
//    /*!
//     * \brief Calculates the restitution coefficient time for a given stiffness, dissipation, and effective mass
//     * \deprecated use species->computeRestitutionCoefficient(2.0*effectiveMass) instead
//     */
//    MERCURY_DEPRECATED
//    Mdouble computeRestitutionCoefficientFromKAndDispAndEffectiveMass(Mdouble k, Mdouble disp, Mdouble mass);
//
//    /*!
//     * \brief Calculates the dissipation for a given stiffness, restitution coefficient, and effective mass
//     */
//    MERCURY_DEPRECATED
//    Mdouble computeDispFromKAndRestitutionCoefficientAndEffectiveMass(Mdouble k, Mdouble r, Mdouble mass);
//
//    /*!
//     * \brief Calculates the collision time for a given stiffness, restitution coefficient, and effective mass
//     * \deprecated use species->computeCollisionTime(2.0*effectiveMass) instead
//     */
//    MERCURY_DEPRECATED
//    Mdouble computeCollisionTimeFromKAndRestitutionCoefficientAndEffectiveMass(Mdouble k, Mdouble r, Mdouble mass);
//
//    /*!
//     * \brief Calculates the dissipation for a given stiffness, collision time, and effective mass
//     * \deprecated use species->setStiffnessAndRestitutionCoefficient(2.0*effectiveMass) instead
//     */
//    MERCURY_DEPRECATED
//    Mdouble computeDispFromKAndCollisionTimeAndEffectiveMass(Mdouble k, Mdouble tc, Mdouble mass);
//
//    /*!
//     * \brief Calculates the restitution coefficient for a given stiffness, collision time, and effective mass
//     * \deprecated use species->computeRestitutionCoefficient(2.0*effectiveMass) instead
//     */
//    MERCURY_DEPRECATED
//    Mdouble computeRestitutionCoefficientFromKAndCollisionTimeAndEffectiveMass(Mdouble k, Mdouble tc, Mdouble mass);
//
//    /*!
//     * \brief Calculates the dissipation for a given collision time, restitution coefficient, and effective mass
//     * \deprecated use species->setCollisionTimeAndRestitutionCoefficient(2.0*effectiveMass) instead
//     */
//    MERCURY_DEPRECATED
//    Mdouble computeDispFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(Mdouble tc, Mdouble r, Mdouble mass);
//
//    /*!
//     * \brief Calculates the stiffness for a given collision time, restitution coefficient, and effective mass
//     * \deprecated use species->setCollisionTimeAndRestitutionCoefficient(2.0*effectiveMass) instead
//     */
//    MERCURY_DEPRECATED
//    Mdouble computeKFromCollisionTimeAndRestitutionCoefficientAndEffectiveMass(Mdouble tc, Mdouble r, Mdouble mass);
//
//    /*!
//     * \brief Calculates the stiffness for a given collision time, dissipation, and effective mass
//     */
//    MERCURY_DEPRECATED
//    Mdouble computeKFromCollisionTimeAndDispAndEffectiveMass(Mdouble tc, Mdouble disp, Mdouble mass);
//
//    /*!
//     * \brief Calculates the resitution coefficient for a given collision time, dissipation, and effective mass
//     * \deprecated use species->computeRestitutionCoefficient(2.0*effectiveMass) instead
//     */
//    MERCURY_DEPRECATED
//    Mdouble computeRestitutionCoefficientFromCollisionTimeAndDispAndEffectiveMass(Mdouble tc, Mdouble disp, Mdouble mass);
//
//    /*!
//     * \brief Calculates the stiffness for a given dissipation, restitution coefficient, and effective mass
//     */
//    MERCURY_DEPRECATED
//    Mdouble computeKFromDispAndRestitutionCoefficientAndEffectiveMass(Mdouble disp, Mdouble r, Mdouble mass);
//
//    /*!
//     * \brief Calculates the collision time for a given dissipation, restitution coefficient, and effective mass
//     * \deprecated use species->computeCollisionTime(2.0*effectiveMass) instead
//     */
//    MERCURY_DEPRECATED
//    Mdouble computeCollisionTimeFromDispAndRestitutionCoefficientAndEffectiveMass(Mdouble disp, Mdouble r, Mdouble mass);

///return type specifically for fuctions returning k, disp, kt, dispt at once

/*!
 * \brief return type specifically for fuctions returning k, disp, kt, dispt at once
 */
class KAndDispAndKtAndDispt
{
public:
    Mdouble k;
    Mdouble disp;
    Mdouble kt;
    Mdouble dispt;
};

/*!
 * \brief Set disp, k, dispt and kt such that is matches a given collision time tc and a normal and tangential
 * restitution coefficient r, beta for a collision of effective/reduced mass m. From Deen...Kuipers2006,
 * eq. 43 and 30
 * \todo what should be used instead of this function?
 */
MERCURY_DEPRECATED
KAndDispAndKtAndDispt
computeDisptFromCollisionTimeAndRestitutionCoefficientAndTangentialRestitutionCoefficientAndEffectiveMass(
        Mdouble tc, Mdouble r, Mdouble beta, Mdouble mass);

/*!
 * \brief Calculates the maximum relative velocity allowed for a normal collision of two particles of radius r and
 * particle mass m (for higher velocities particles could pass through each other)
 * \todo what should be used instead of this function?
 */
MERCURY_DEPRECATED
Mdouble getMaximumVelocity(Mdouble k, Mdouble disp, Mdouble radius, Mdouble mass);

/*!
 * \brief Returns the correct saveCount if the total number of saves, the final time and the time step is known
 */
unsigned int getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(unsigned int numberOfSaves, Mdouble timeMax,
                                                                Mdouble timeStep);

/*!
 * \brief Reads a line from one stringstream into another, and prepares the latter for reading in.
 */
void getLineFromStringStream(std::istream& in, std::stringstream& out);

/*!
 * \brief Writes a string to a file.
 */
bool writeToFile(std::string filename, std::string filecontent);

/*!
 * \brief Writes a string to a file.
 */
void writeCommandLineToFile(const std::string filename, const int argc, char * const argv[]);

/*!
 * \brief Plots to a gnuplot window.
 */
void gnuplot(std::string command);

/*!
 * \brief Adds a string to an existing file.
 */
bool addToFile(std::string filename, std::string filecontent);

/*!
 * \brief Function to check if a file exists, is used to check if a run has already need done.
 */
bool fileExists(std::string strFilename);

/*!
 * \brief Provides a simple interface for opening a file.
 */
bool openFile(std::fstream& file, std::string filename, std::fstream::openmode mode);

/*!
 * \brief Calculates the effective mass of a particle pair, i.e. half the harmonic mean of two particle masses.
 */
Mdouble getEffectiveMass(Mdouble mass0, Mdouble mass1);

std::vector<double> readArrayFromFile(std::string filename, int& n, int& m);

void more(std::string filename, unsigned nLines = constants::unsignedMax);

template<typename T>
std::string to_string(const T& n)
{
    std::ostringstream stm;
    stm << n;
    return stm.str();
}

std::string to_string(Mdouble value, unsigned precision);

/**
 * \brief Reads optional variables in the restart file
 *
 * \details A variable is stored in the restart file by storing the variables name and value, e.g.
 *   " name value"
 * If a variable is always written to the restart file, it can be read-in like this:
 *   is >> dummy >> variable;
 * If a variable is optional, the variable name has to be checked, and should be read in like this:
 *   readOptionalVariable(is,"name",variable);
 */
template<typename T>
bool readOptionalVariable(std::istream& is, const std::string& name, T& variable)
{
    ///\todo readOptionalVariable should check the full variable name, not just the next value.
    /// However, I don't know how to put the location in the ifstream back.
    const auto pos = is.tellg();
    std::string dummy;
    is >> dummy;
    if (dummy == name)
    {
        is >> variable;
        return true;
    } else {
        is.seekg(pos);
        return false;
    }
}

/**
 * Used to test the loading/unloading properties of a contact law
 **/
void loadingTest(const ParticleSpecies* species, Mdouble displacement, Mdouble velocity, Mdouble radius,
                 std::string name);

/**
 * Used to test the tangential loading/unloading/reloading properties of a frictional contact law
 **/
void normalAndTangentialLoadingTest(const ParticleSpecies* species, Mdouble displacement,
                                    Mdouble tangentialDisplacement, Mdouble velocity, Mdouble radius,
                                    std::string name);

/**
 * Used to test the tangential loading/unloading/reloading properties of a frictional contact law
 **/
void objectivenessTest(const ParticleSpecies* species, Mdouble displacement, Mdouble tangentialDisplacement,
                       Mdouble velocity, Mdouble radius, std::string name);

/**
 * \brief Checks if the next argument in the input stream is a certain string
 */
bool compare(std::istream& is, std::string s);


/**
 * Allows a quick read-in from a parameter file.
 *
 * For example, the following code reads in the variable time from the file in:
 * \code{.cpp}
 *   double time = readFromFile("in","time",24);
 * \endcode
 *
 * The in file needs to contain the string time, followed by the value, e.g.
 * \code{.sh}
 *   time 20
 * \endcode
 *
 * If the file cannot be opened, or the parameter is not found,
 * the variable is set to the default value specified by the third argument.
 *
 * @param fileName name of input
 * @param varName  variable name as it appears in the input file
 * @param value    default value (used if the parameter could not be read)
 * @return         value of variable
 */
template<typename T>
T readFromFile(std::string fileName, std::string varName, T value)
{
    //open filestream
    std::ifstream is(fileName.c_str(), std::ios::in);
    if (is.fail())
    {
        logger(INFO, "readFromFile: file % could not be opened, variable % set to default value %",
               fileName, varName, value);
        return value;
    }
    
    //read in variables, until the right one is fount
    std::string s;
    while (!is.eof())
    {
        is >> s;
        if (!s.compare(varName))
        {
            is >> value;
            logger(INFO, "readFromFile: variable % set to % ", varName, value);
            return value;
        }
    }
    
    //if the right variable is never found
    logger(WARN, "readFromFile: variable % not set in file %, using default value % ", varName, fileName, value);
    return value;
}

/**
 * Returns true if command line arguments contain varName, false else
 * Usage example:
 *    if (readFromCommandLine(argc, argv, '-verbose')) ...
 * @param argc pass through number of command line arguments
 * @param argv pass through values of command line arguments
 * @param varName name of command line arguments that is required to return true
 * @return true or false
 */
bool readFromCommandLine(int argc, char *argv[], std::string varName);

template<typename T>
T readFromCommandLine(int argc, char *argv[], std::string varName, T value)
{
    for (unsigned i=0; i<argc-1; ++i) {
        if (varName == argv[i]) {
            value = atof(argv[i+1]);
            logger(INFO, "readFromCommandLine: % set to % ", varName.substr(1), value);
            return value;
        }
    }
    //if the variable is not found
    logger(INFO, "readFromCommandLine: % set to default value % ", varName.substr(1), value);
    return value;
}

template<typename T, size_t n>
std::array<T,n> readArrayFromCommandLine(int argc, char *argv[], std::string varName, std::array<T,n> value)
{
    for (unsigned i=0; i<argc-1; ++i) {
        if (varName == argv[i]) {
            unsigned j = i+1;
            std::stringstream out;
            for (auto& v : value) {
                v = atof(argv[j]);
                out << v << ' ';
                ++j;
            }
            logger(INFO, "readFromCommandLine: % set to % ", varName.substr(1), out.str());
            return value;
        }
    }
    //if the variable is not found
    std::stringstream out;
    for (auto& v : value) out << v << ' ';
    logger(INFO, "readFromCommandLine: % set to default value % ", varName.substr(1), out.str());
    return value;
}

template<typename T>
std::vector<T> readVectorFromCommandLine(int argc, char *argv[], std::string varName, size_t n, std::vector<T> values)
{
    for (unsigned i=0; i<argc-1; ++i) {
        if (varName == argv[i]) {
            // read until the next argument starts
            values.resize(0);
            std::stringstream out;
            for (int j = i+1; j<argc and argv[j][0]!='-'; ++j) {
                values.push_back(atof(argv[j]));
                out << values.back() << ' ';
            }
            logger(INFO, "readFromCommandLine: % set to % ", varName.substr(1), out.str());
            return values;
        }
    }
    //if the variable is not found
    std::stringstream out;
    for (auto& v : values) out << v << ' ';
    logger(INFO, "readFromCommandLine: % set to default value % ", varName.substr(1), out.str());
    return values;
}

template<>
std::string readFromCommandLine<std::string>(int argc, char *argv[], std::string varName, std::string value);

void check(double real, double ideal, double error, std::string errorMessage);

void check(Vec3D real, Vec3D ideal, double error, std::string errorMessage);

void check(Matrix3D real, Matrix3D ideal, double error, std::string errorMessage);

void check(MatrixSymmetric3D real, MatrixSymmetric3D ideal, double error, std::string errorMessage);

std::string getPath();

Mdouble getRealTime();

bool isNext(std::istream& is, const std::string name);

bool createDirectory(std::string);

Mdouble round(const Mdouble value, unsigned precision);

/*
 * \brief Returns the Rayleigh time step for a Hertz contact law.
 * \detailed An accepted time step for Hertz is 10-20% of the Rayleigh time step.
 * See \cite Marigo2015
 */
Mdouble getRayleighTime(Mdouble radius, Mdouble shearModulus, Mdouble poisson, Mdouble density);

}

#endif
