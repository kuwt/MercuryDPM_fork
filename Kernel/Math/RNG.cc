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

#include "RNG.h"
#include "Helpers.h"
#include <limits>

/**
 * This is a random number generator and returns a Mdouble within the range specified
 * \todo{Thomas: This code does sth. when min>max; I would prefer to throw an error.}
 * \todo{the random seed should be stored in restart}
 * */

RNG::RNG()
{
    randomSeedLinearCongruentialGenerator_ = 0;
    a_ = 1103515245;
    c_ = 12345;
    m_ = 1024 * 1024 * 1024;
    type_ = RNGType::LAGGED_FIBONACCI_GENERATOR;
    p_ = 607;
    q_ = 273;
    randomSeedLaggedFibonacciGenerator_.resize(p_);
    seedLaggedFibonacciGenerator();

    haveSavedBoxMuller_ = false;
    savedBoxMuller_ = 0;

}

void RNG::setRandomSeed(unsigned long int new_seed)
{
    randomSeedLinearCongruentialGenerator_ = new_seed;
    seedLaggedFibonacciGenerator();
}

void RNG::read(std::istream& is)
{
    std::string dummy;
    unsigned int type;
    is >> type;
    type_ = static_cast<RNGType>(type);
    is >> a_;
    is >> c_;
    is >> m_;
    is >> p_;
    is >> q_;
    is >> randomSeedLinearCongruentialGenerator_;
    //note: the seeds for the LaggedFibonacciGenerator cannot be restarted currently.
    seedLaggedFibonacciGenerator();
    //randomSeedLaggedFibonacciGenerator_.resize(p_);
    //for (auto& v : randomSeedLaggedFibonacciGenerator_) 
    //    is >> v;
}

void RNG::write(std::ostream& os) const
{
    os << " " << static_cast<unsigned int>(type_);
    os << " " << a_;
    os << " " << c_;
    os << " " << m_;
    os << " " << p_;
    os << " " << q_;
    os << " " << randomSeedLinearCongruentialGenerator_;
    //for (auto v : randomSeedLaggedFibonacciGenerator_) 
    //    os << " " << v;
}

void RNG::setLinearCongruentialGeneratorParmeters(unsigned const int a, unsigned const int c, unsigned const int m)
{
    a_ = a;
    c_ = c;
    m_ = m;
}

void RNG::randomise()
{
#ifdef MERCURY_USE_MPI
    //First set a random seed on the root
    if (PROCESSOR_ID == 0)
    {
        setRandomSeed(static_cast<unsigned long int>(time(nullptr)));
    }

    //Communicate this to the rest of the processes
    std::vector<int> values(7);
    if (PROCESSOR_ID == 0)
    {
        values[0] = static_cast<unsigned int>(type_);
        values[1] = a_;
        values[2] = c_;
        values[3] = m_;
        values[4] = p_;
        values[5] = q_;
        values[6] = randomSeedLinearCongruentialGenerator_;
    } 
    MPIContainer::Instance().broadcast(values.data(),7,0);

    //Update the generators on the other processors   
    if (PROCESSOR_ID != 0)
    {
        type_ = static_cast<RNGType>(values[0]);
        a_ = values[1];
        c_ = values[2];
        m_ = values[3];
        p_ = values[4];
        q_ = values[5];
        randomSeedLinearCongruentialGenerator_ = values[6];
    }
    seedLaggedFibonacciGenerator();
#else
    setRandomSeed(static_cast<unsigned long int>(time(nullptr)));
#endif
}

void RNG::setRandomNumberGenerator(RNGType type)
{
    type_ = type;
}

Mdouble RNG::getRandomNumber()
{
    return getRandomNumber(0, 1);
}

Mdouble RNG::getRandomNumber(Mdouble min, Mdouble max)
{
    logger.assert(min <= max, "getRandomNumber: min cannot be larger than max");
    if (type_ == RNGType::LINEAR_CONGRUENTIAL_GENERATOR) {
        return getRandomNumberFromLinearCongruentialGenerator(min, max);
    } else {
        return getRandomNumberFromLaggedFibonacciGenerator(min, max);
    }
}

/* \details This uses the Box--Muller transform
 *   https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
 * The Box--Muller algorithm generates a pair of independently distributed
 * normal random variables. One of these will be returned straight away and the
 * other will be stored and returned the next time you call this function.
 */
Mdouble RNG::getNormalVariate()
{
    static const double epsilon = std::numeric_limits<Mdouble>::min();

    if (haveSavedBoxMuller_)
    {
        /* If we have already generated a normal variate, use it. */
        haveSavedBoxMuller_ = false;
        return savedBoxMuller_;
    }
    else
    {
        /* Otherwise, generate a pair of normal variates, return one of them,
         * and save the other. */
        Mdouble radius, theta;
        do
        {
            radius = getRandomNumber(0, 1);
            theta = getRandomNumber(0, 2 * constants::pi);
        } while (radius <= epsilon);
        // make sure that the radius generated is not too small
        // (unlikely to happen, just a safety check)

        savedBoxMuller_ = sqrt(-2.0 * log(radius)) * sin(theta);
        haveSavedBoxMuller_ = true;
        return sqrt(-2.0 * log(radius)) * cos(theta);
    }
}

Mdouble RNG::getNormalVariate(Mdouble mean, Mdouble stdev)
{
    if (stdev == 0) {
        logger(WARN,
               "[RNG::getNormalVariate(Mdouble, Mdouble)] Zero stdev?");
        return mean;
    } else if (stdev < 0) {
        logger(ERROR,
               "[RNG::getNormalVariate(Mdouble, Mdouble)] Negative stdev is not allowed.");
        exit(-1);
    } else {
        return getNormalVariate() * stdev + mean;
    }
}

/*!
 * \details This uses Knuth's algorithm for generating Poisson variates. It's
 * simple but slow for large values of lambda --- beware. 
 */
unsigned int RNG::getPoissonVariate(Mdouble lambda)
{
    if (lambda > 50)
    {
        logger(WARN, "[RNG::getPoissonVariate(Mdouble)] Knuth's algorithm for Poissons may be slow for lambda = %", lambda);
    }
    unsigned int k = 0; 
    Mdouble p = 1; 
    Mdouble u;
    do 
    {
        k++;
        u = getRandomNumber(0, 1);
        p *= u;
    }
    while (u > exp(-lambda));
    return k-1;
}

/**
 * This is a basic Linear Congruential Generator Random
 * Is described by three parameters, the multiplication a, the addition c and the mod m
 */
#pragma optimize( "", off )
Mdouble RNG::getRandomNumberFromLinearCongruentialGenerator(Mdouble min, Mdouble max)
{
    //Update the random seed
    randomSeedLinearCongruentialGenerator_ = (a_ * randomSeedLinearCongruentialGenerator_ + c_) % m_;
    
    //Generate a random number in the required range
    
    Mdouble range = max - min;
    Mdouble random_num = min + range * randomSeedLinearCongruentialGenerator_ / (static_cast<Mdouble>(m_) + 1.0);
    
    return random_num;
}
#pragma optimize( "", on )

/**************************************
 * This sets the seed for LFG using LCG
 * 
 * 
 **************************************/
void RNG::seedLaggedFibonacciGenerator()
{
    for (unsigned int i = 0; i < p_; i++)
    {
        randomSeedLaggedFibonacciGenerator_[i] = getRandomNumberFromLinearCongruentialGenerator(0, 1.0);
    }
}

/**
 * This is a basic Linear Fibonacci Generator Random
 * Is described by three parameters, the multiplication a, the addition c and the mod m
 */
Mdouble RNG::getRandomNumberFromLaggedFibonacciGenerator(Mdouble min, Mdouble max)
{
#pragma optimize( "", off )
    Mdouble new_seed = fmod(randomSeedLaggedFibonacciGenerator_[0] + randomSeedLaggedFibonacciGenerator_[p_ - q_],
                            static_cast<Mdouble>(1.0));
    //Update the random seed
    randomSeedLaggedFibonacciGenerator_.erase(randomSeedLaggedFibonacciGenerator_.begin());
    randomSeedLaggedFibonacciGenerator_.emplace_back(new_seed);
    
    //Generate a random number in the required range
    
    Mdouble random_num;
    
    Mdouble range = max - min;
    random_num = min + range * new_seed;
    return random_num;
#pragma optimize( "", on )
}

/**
 * This function tests the quality of random numbers, based on the chi-squared test. 
 * It reports a probability that the random number being generated are coming from a uniform distributed. 
 * If this number is less than 0.95, it is strongly advised that you change the parameters being used
 * */
Mdouble RNG::test()
{
    //This are the fixed parameters that define the test
    static unsigned int num_of_tests = 100000;
    static Mdouble max_num = 100.0;
    static unsigned int num_of_bins = 10;
    
    //This is the generated random_number
    Mdouble rn;
    //This is the bin the random number will lie in
    unsigned int bin = 0;
    //This is a vector of bins
    std::vector<int> count;
    count.resize(num_of_bins);
    
    //Initialisation of the bins
    for (unsigned int i = 0; i < num_of_bins; i++)
    {
        count[bin] = 0;
    }
    
    //Loop over a number of tests
    for (unsigned int i = 0; i < num_of_tests; i++)
    {
        rn = getRandomNumber(0.0, max_num);
        bin = static_cast<unsigned int>(std::floor(rn * num_of_bins / max_num));
        
        //Add one to the bin count
        count[bin]++;
        
    }
    
    //Final post-process the result and report on the random number
    Mdouble chi_cum = 0.0;
    Mdouble expected = num_of_tests / num_of_bins;
    
    for (unsigned int i = 0; i < num_of_bins; i++)
    {
        chi_cum = chi_cum + (count[i] - expected) * (count[i] - expected) / expected;
        std::cout << i << " : "
                  << count[i] << " : " << (count[i] - expected) * (count[i] - expected) / expected << std::endl;
    }
    //end for loop over computing the chi-squared value.
    std::cout << "chi_cum " << chi_cum << std::endl;

    return mathsFunc::chi_squared_prob(chi_cum, num_of_bins);
}

///This function sets the parameters for the LFG random number generator
void RNG::setLaggedFibonacciGeneratorParameters(const unsigned int p, const unsigned int q)
{
    //p must be greater than q so makes sure this is true. Not sure what happens if you set p=q, in the LFG alogrithm.
    if (p < q)
    {
        p_ = q;
        q_ = p;
    }
    
    randomSeedLaggedFibonacciGenerator_.resize(p_);
    seedLaggedFibonacciGenerator();
}
