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

#ifndef MERCURY_PSD_H
#define MERCURY_PSD_H

#include <fstream>
#include <vector>

/**
 * Stores a radius and a cumulative number density:
 * To be used as a vector, std::vector<PSD> psd
 *
 * Cumulative number density: nc_i is the number fraction of particles whose radius is less than r_i.
 * It requires vc_0=0, vc_end=1. This variable is used as psd in MercuryDPM, as it can be interpreted as the probability that a particle's radius is is less than r_i:
 *     p(r<r_i) = nc_i.
 *
 * Cumulative volume density: vc_i is the volume fraction of particles whose radius is less than r_i. It requires vc_0=0, vc_end=1
 *
 * Subtractive volume density: v_i is the volume fraction of particles whose radius is between r_i-1 and r_i. It requires v_0=0, sum(v_i)=1
 *
 * Sieve data: s_i is the volume fraction of particles whose radius is between r_i and r_i+1. It requires sum(s_i)=1. Property: s_i = v_i-1
 */
struct PSD
{
//    // sets the psd from the probability values given in the file batch_materials_PSD.xlsx
//    static void checkPSD (std::vector<PSD>& psd) {
//        if (psd.empty()) {
//            logger(ERROR,"psd is empty");
//        }
//        if (psd[0].probability!=0) {
//            logger(WARN,"psd should start with a zero probability");
//            psd.emplace(psd.begin(),PSD{psd[0].radius,0});
//        }
//        if (psd.back().probability!=1) {
//            logger(WARN,"psd should end with a unit probability");
//            psd.emplace_back(PSD{psd.back().radius,1});
//        }
//    }

    // validate whether psd is a cumulative distribution
    static void validateCumulativeDistribution (std::vector<PSD>& psd);

    // combine vectors of radii and probabilities into psd vector
    static std::vector<PSD> createPSDFromRadiiAndProbabilities (const std::vector<double>& radius,
            const std::vector<double>& probability);

    // converts a subtractive to a cumulative psd
    static void convertSubtractiveToCumulative (std::vector<PSD>& psd);

    // converst a single to a cumulative psd
    static void convertCumulativeToSubtractive (std::vector<PSD>& psd);

    // converts a subtractive number psd to a subtractive volume psd
    static void convertSubtractiveVolumeToNumber (std::vector<PSD>& psd);

    // converts a cumulative number psd to a cumulative volume psd
    static void convertCumulativeVolumeToNumber (std::vector<PSD>& psd);

    // converts a subtractive number psd to a subtractive volume psd
    static void convertSubtractiveNumberToVolume (std::vector<PSD>& psd);

    // converts a cumulative number psd to a cumulative volume psd
    static void convertCumulativeNumberToVolume (std::vector<PSD>& psd);

    // returns a cumulative number particle size distribution from a given volume particle size distribution
    static std::vector<PSD> setFromVolumePSD (const std::vector<double>& radius,
                                                          const std::vector<double>& probability);

    // get quantile size
    static double getD0 (const std::vector<PSD>& psd){
        return 2.0*psd[0].radius;
    }

    static double getD10 (const std::vector<PSD>& psd){
        return 2.0*getQuantile(psd,0.1);
    }

    static double getD50 (const std::vector<PSD>& psd){
        return 2.0*getQuantile(psd,0.5);
    }

    static double getD90 (const std::vector<PSD>& psd){
        return 2.0*getQuantile(psd,0.9);
    }

    static double getD100 (const std::vector<PSD>& psd){
        return 2.0*psd.back().radius;
    }


    // get median size
    static double getQuantile(std::vector<PSD> psd, double quantile);

    // get radius r such that a monodisperse system has the same number of particles as a polydisperse system
    static double getVolumetricMean(std::vector<PSD> psd);

    /**
     * required to use std::lower_bound for finding when the probability is higher than a certain value
     *   std::vector<PSD> psd;
     *   double probability;
     *   std::lower_bound(psd.begin(),psd.end(),probability);
     * @param probability
     * @return
     */
    bool operator<(const double probability) const
    {
        return this->probability < probability;
    }
    
    /*!
     * \brief Writes to output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const PSD& psd) {
        os << psd.radius << ' ' << psd.probability;
        return os;
    }
    
    /*!
     * \brief Reads from input stream
     */
    friend std::istream& operator>>(std::istream& is, PSD& psd) {
        is >> psd.radius >> psd.probability;
        return is;
    }
    
    double radius;
    double probability;
};



#endif //MERCURY_PSD_H
