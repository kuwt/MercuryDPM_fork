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

#ifndef MECURYDPM_RNG_H
#define MECURYDPM_RNG_H

//Used by the randomise function
#include <ctime>
#include <cmath>
#include <vector>
#include "ExtendedMath.h"

/*!
 * \brief 
 */
enum class RNGType : unsigned char
{
    LAGGED_FIBONACCI_GENERATOR = 0,
    LINEAR_CONGRUENTIAL_GENERATOR = 1
};

/*!
 * \class RNG
 * \brief This is a class that generates random numbers i.e. named the Random Number Generator (RNG). 
 * \details This is a stand-along class; but is encapsulated (used) by the MD class.
 * To make it architecture safe the both LCG and function is hard codes i.e. does not use the internal C++ one.
 * \todo (AT) implement new C++-standard RNG instead of this one (Kudos on the hard work done here though ;).
 * NB: maybe something for Mercury 2?
 */
class RNG
{
public:
    /*!
     * \brief default constructor
     */
    RNG();
    
    /*!
     * \brief This is the seed for the random number generator (note the call to seed_LFG is only required really
     * if using that type of generator, but the other one is always required)
     */
    void setRandomSeed(unsigned long int new_seed);
    
    void read(std::istream& is);

    void write(std::ostream& os) const;
    
    /*!
     * \brief This is a random generating routine can be used for initial positions
     */
    Mdouble getRandomNumber();
    
    Mdouble getRandomNumber(Mdouble min, Mdouble max);

    /*!
     * \brief Shorthand for getRandomNumber(min, max)
     */
    Mdouble operator()(Mdouble min, Mdouble max)
    {
        return getRandomNumber(min, max);
    }

    /*
     * \brief Shorthand for getRandomNumber()
     */
    Mdouble operator()()
    {
        return getRandomNumber(0.0, 1.0);
    }

    /*!
     * \brief Produces a random number according to a normal distribution with
     * mean 0 and standard deviation 1.
     */
    Mdouble getNormalVariate();

    /*! 
     * \brief Produces a random number according to a normal distribution.
     */
    Mdouble getNormalVariate(Mdouble mean, Mdouble stdev);

    /*!
     * \brief Produces a random number according to a Poisson distribution
     */
    unsigned int getPoissonVariate(Mdouble lambda);

    /*!
     * \brief This function tests the quality of random numbers, based on the chi-squared test. 
     */
    Mdouble test();

    /*!
     * \brief This functions set the parameters for the LCG random number generator. It goes multiplier, addition, mod.
     */
    void setLinearCongruentialGeneratorParmeters(const unsigned int a, const unsigned int c, unsigned int m);
    
    /*!
     * \brief sets the random variables such that they differ for each run
     */
    void randomise();
    
    /*!
     * \brief This function sets the parameters for the LFG random number generator
     */
    void setLaggedFibonacciGeneratorParameters(const unsigned int p, const unsigned int q);

    /*!
     * \brief Allows the user to set which random number generator is used
     */
    void setRandomNumberGenerator(RNGType type);

private:
    
    /*!
     * \brief This is a basic Linear Congruential Generator Random 
     */
    Mdouble getRandomNumberFromLinearCongruentialGenerator(Mdouble min, Mdouble max);

    /*!
     * \brief This is a Lagged Fibonacci Generator
     */
    Mdouble getRandomNumberFromLaggedFibonacciGenerator(Mdouble min, Mdouble max);

    /*!
     * \brief This seed the LFG
     */
    void seedLaggedFibonacciGenerator();

    /*!
     * \brief This is the initial seed of the RNG
     */
    unsigned long int randomSeedLinearCongruentialGenerator_;

    /*!
     * \brief This is the seeds required for the LFG
     */
    std::vector<Mdouble> randomSeedLaggedFibonacciGenerator_;

    /*!
     * \brief This are the two parameters that control the LCG random generated
     */
    unsigned long int a_, c_, m_;

    /*!
     * \brief This are the parameters that control the LFG random generator
     */
    unsigned long int p_, q_;

    /*!
     * \brief This is the type of random number generator
     */
    RNGType type_;

    /*
     * \todo JMFT: These Box--Muller private variables are currently not
     * read/write from .restart files. This leads to undefined (unrepeatable)
     * behaviour when restarting if you use normal variates.
     */

    /*!
     * \brief A flag that keeps track of whether or not to generate a new pair of normal variates (using Box--Muller)
     */
    bool haveSavedBoxMuller_;

    /*!
     * \brief A storage space for the so-far-unused variate from the pair generated by Box--Muller
     */
    Mdouble savedBoxMuller_;
    
};

#endif
