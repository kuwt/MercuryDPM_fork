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

#include "Logger.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Species/LinearPlasticViscoelasticSpecies.h"
#include "Math/ExtendedMath.h"
#include "Math/Helpers.h"
#include "iomanip"

int main(int argc UNUSED, char *argv[] UNUSED)
{
    logger(VERBOSE, "Running HelperFunctionsUnitTest\n"
                    "  Testing getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep");
    
    unsigned int saveCount = helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(21, 10, 0.001);
    if (saveCount != 501)
    {
        logger(ERROR, "save count is %, but should be %", saveCount, 501);
    }
    
    logger(VERBOSE, "  Testing getLineFromStringStream");
    std::stringstream is("1.0\n2.000000000000001\n3.0");
    std::stringstream line;
    helpers::getLineFromStringStream(is, line);
    Mdouble value;
    line >> value >> value; //the second value is not reading from the second line
    if (value != 1.0)
    {
        logger(ERROR, "value is %, but should be %", value, 0.0);
    }
    helpers::getLineFromStringStream(is, line);
    line >> value; //reads double precision
    //std::cout << std::setprecision(16) << value << std::endl;
    if (value != 2.000000000000001)
    {
        logger(ERROR, "value is %, but should be %", value, 2.000000000000001);
    }
    
    logger(VERBOSE, "  Testing writeToFile, open");
    helpers::writeToFile("HelperFunctionsUnitTest.ini",
                         "1 0 0 0 0 1 1 1\n"
                         "0.5 0.5 0  0 0 0.5  0  0 0 0  0 0 0  0\n");
    std::fstream file;
    helpers::openFile(file, "HelperFunctionsUnitTest.ini", std::fstream::in);
    value = 0.0;
    file >> value;
    if (value != 1.0)
    {
        logger(ERROR, "value is %, but should be %", value, 1.0);
    }
    
    logger(VERBOSE, "  Testing writeToFile, open");
    bool fileExists = helpers::fileExists("HelperFunctionsUnitTest.ini");
    if (fileExists != true)
    {
        logger(ERROR, "HelperFunctionsUnitTest.ini exists, but is not detected");
    }
    fileExists = helpers::fileExists("HelperFunctionsUnitTest.out");
    if (fileExists != false)
    {
        logger(ERROR, "HelperFunctionsUnitTest.out does not exist, but is detected");
    }
    
    logger(VERBOSE, "  Testing getEffectiveMass");
    Mdouble effectiveMass = helpers::getEffectiveMass(1.0, 2.0);
    if (!mathsFunc::isEqual(effectiveMass, 2.0 / 3.0, 1e-10))
    {
        logger(ERROR, "effective mass is %, but should be %", effectiveMass, 2.0 / 3.0);
    }
    
    
    
    logger(INFO, " Testing removeFromCommandline");
    char arguments[][20] = { "programname", "-test1", "-test2", "1"};
    int argcTest = 4;
    char** argvTest;
    argvTest = new char*[argcTest];
    
    for (unsigned int i=0; i<argcTest; i++)
    {
        argvTest[i] = &arguments[i][0];
    }
    
    if (!(helpers::readFromCommandLine(argcTest, argvTest, "-test1") && helpers::removeFromCommandline(argcTest, argvTest, "-test1", 0) && !helpers::readFromCommandLine(argcTest, argvTest, "-test1"))
        || !(helpers::readFromCommandLine(argcTest, argvTest, "-test2", 2) == 1 && helpers::removeFromCommandline(argcTest, argvTest, "-test2", 1) && helpers::readFromCommandLine(argcTest, argvTest, "-test2", 2) == 2)
        || argcTest>1 )
    {
        logger(ERROR, "helpers::removeFromCommandline not working as expected");
    }
    
    
    logger(VERBOSE, "Running LinearViscoelasticSpecies helper functions unit test");
    Mdouble realStiffness = 2e5;
    Mdouble realDissipation = 25.0;
    Mdouble mass = 1.0; //particle mass
    Mdouble radius = 0.5; //particle radius
    Mdouble realCollisionTime = 0.004971179385062563;
    Mdouble realRestitution = 0.883132984295725;
    Mdouble realMaximumVelocity = 316.227766016838;
    
    LinearViscoelasticSpecies species0;
    species0.setStiffness(realStiffness);
    species0.setDissipation(realDissipation);
    
    logger(VERBOSE, "  Testing getMaximumVelocity");
    Mdouble maximumVelocity = species0.getMaximumVelocity(radius, mass);
    //std::cout << std::setprecision(16) << maximumVelocity << std::endl;
    if (!mathsFunc::isEqual(maximumVelocity, realMaximumVelocity, 1e-10))
    {
        logger(ERROR, "maximum velocity is %, but should be %", maximumVelocity, realMaximumVelocity);
    }
    
    logger(VERBOSE, "  Testing getCollisionTime, getRestitutionCoefficient");
    Mdouble collisionTime = species0.getCollisionTime(mass);
    Mdouble restitution = species0.getRestitutionCoefficient(mass);
    if (!mathsFunc::isEqual(collisionTime, realCollisionTime, 1e-10))
    {
        logger(ERROR, "collision time is %, but should be %", collisionTime, realCollisionTime);
    }
    if (!mathsFunc::isEqual(restitution, realRestitution, 1e-10))
    {
        logger(ERROR, "restitution coefficient is %, but should be %", restitution, realRestitution);
    }
    
    logger(VERBOSE, "  Testing setStiffnessAndRestitutionCoefficient");
    species0.setStiffnessAndRestitutionCoefficient(realStiffness, realRestitution, mass);
    if (!mathsFunc::isEqual(species0.getStiffness(), realStiffness, 1e-10))
    {
        logger(ERROR, "stiffness is %, but should be %", species0.getStiffness(), realStiffness);
    }
    if (!mathsFunc::isEqual(species0.getDissipation(), realDissipation, 1e-10))
    {
        logger(ERROR, "dissipation is %, but should be %", species0.getDissipation(), realDissipation);
    }
    
    logger(VERBOSE, "  Testing setCollisionTimeAndRestitutionCoefficient");
    species0.setCollisionTimeAndRestitutionCoefficient(realCollisionTime, realRestitution, mass);
    if (!mathsFunc::isEqual(species0.getStiffness(), realStiffness, 1e-10))
    {
        logger(ERROR, "stiffness is %, but should be %", species0.getStiffness(), realStiffness);
    }
    if (!mathsFunc::isEqual(species0.getDissipation(), realDissipation, 1e-10))
    {
        logger(ERROR, "dissipation is %, but should be %", species0.getDissipation(), realDissipation);
    }
    
    logger(VERBOSE, "  Testing setCollisionTimeAndRestitutionCoefficient for two masses");
    species0.setCollisionTimeAndRestitutionCoefficient(realCollisionTime, realRestitution, mass, mass);
    if (!mathsFunc::isEqual(species0.getStiffness(), realStiffness, 1e-10))
    {
        logger(ERROR, "stiffness is %, but should be %", species0.getStiffness(), realStiffness);
    }
    if (!mathsFunc::isEqual(species0.getDissipation(), realDissipation, 1e-10))
    {
        logger(ERROR, "dissipation is %, but should be %", species0.getDissipation(), realDissipation);
    }
    
    logger(VERBOSE, "  Testing copy constructor");
    LinearViscoelasticSpecies species1 = species0;
    if (!mathsFunc::isEqual(species1.getStiffness(), realStiffness, 1e-10))
    {
        logger(ERROR, "stiffness is %, but should be %", species1.getStiffness(), realStiffness);
    }
    if (!mathsFunc::isEqual(species1.getDissipation(), realDissipation, 1e-10))
    {
        logger(ERROR, "dissipation is %, but should be %", species1.getDissipation(), realDissipation);
    }
    
    logger(VERBOSE, "  Testing mix");
    LinearViscoelasticMixedSpecies species2;
    species2.mixAll(&species0, &species1);
    if (!mathsFunc::isEqual(species2.getStiffness(), realStiffness, 1e-10))
    {
        logger(ERROR, "stiffness is %, but should be %", species2.getStiffness(), realStiffness);
    }
    if (!mathsFunc::isEqual(species2.getDissipation(), realDissipation, 1e-10))
    {
        logger(ERROR, "dissipation is %, but should be %", species2.getDissipation(), realDissipation);
    }
    
    logger(VERBOSE, "Running LinearPlasticViscoelasticSpecies helper functions unit test");
    LinearPlasticViscoelasticSpecies species3;
    
    logger(VERBOSE, "  Testing setPlasticParameters");
    species3.setPlasticParameters(realStiffness, 2.0 * realStiffness, 0.5 * realStiffness, 0.5);
    if (!mathsFunc::isEqual(species3.getLoadingStiffness(), realStiffness, 1e-10))
    {
        logger(ERROR, "stiffness is %, but should be %", species3.getLoadingStiffness(), realStiffness);
    }
    if (!mathsFunc::isEqual(species3.getUnloadingStiffnessMax(), 2.0 * realStiffness, 1e-10))
    {
        logger(ERROR, "max. unloading stiffness is %, but should be %", species3.getUnloadingStiffnessMax(),
               2.0 * realStiffness);
    }
    if (!mathsFunc::isEqual(species3.getCohesionStiffness(), 0.5 * realStiffness, 1e-10))
    {
        logger(ERROR, "max. unloading stiffness is %, but should be %", species3.getCohesionStiffness(),
               0.5 * realStiffness);
    }
    if (!mathsFunc::isEqual(species3.getPenetrationDepthMax(), 0.5, 1e-10))
    {
        logger(ERROR, "dissipation is %, but should be %", species3.getPenetrationDepthMax(), 0.5);
    }
    
    logger(VERBOSE, "  Testing setCollisionTimeAndRestitutionCoefficient");
    species3.setCollisionTimeAndRestitutionCoefficient(realCollisionTime, realRestitution, mass);
    if (!mathsFunc::isEqual(species3.getLoadingStiffness(), realStiffness, 1e-10))
    {
        logger(ERROR, "stiffness is %, but should be %", species3.getLoadingStiffness(), realStiffness);
    }
    if (!mathsFunc::isEqual(species3.getUnloadingStiffnessMax(), realStiffness, 1e-10))
    {
        logger(ERROR, "max. unloading stiffness is %, but should be %", species3.getUnloadingStiffnessMax(),
               realStiffness);
    }
    if (!mathsFunc::isEqual(species3.getDissipation(), realDissipation, 1e-10))
    {
        logger(ERROR, "dissipation is %, but should be %", species3.getDissipation(), realDissipation);
    }
    
    logger(VERBOSE, "  Testing computeTimeStep");
    Mdouble timeStep = species3.computeTimeStep(mass);
    if (!mathsFunc::isEqual(timeStep, 0.02 * realCollisionTime, 1e-10))
    {
        logger(ERROR, "time step is %, but should be %", timeStep, 0.02 * realCollisionTime);
    }
}
