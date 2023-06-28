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

#include "PSDContinuous.h"
#include <Logger.h>
#include "ExtendedMath.h"
using mathsFunc::square;
using mathsFunc::cubic;

void PSDContinuous::print (std::vector<PSDContinuous>& psd) {
    if (psd.front().radius>1e-1) {
        for (const auto p : psd) {
            std::cout << p.radius << "m\t" << p.probability * 100 << "\%\n";
        }
    } else if (psd.front().radius>1e-4) {
        for (const auto p : psd) {
            std::cout << p.radius * 1000 << "mm\t" << p.probability * 100 << "\%\n";
        }
    } else {
        for (const auto p : psd) {
            std::cout << p.radius * 1e6 << "um\t" << p.probability * 100 << "\%\n";
        }
    }
}

// validate whether psd is a cumulative distribution
void PSDContinuous::validateCumulativeDistribution (std::vector<PSDContinuous>& psd) {
    logger.assert_always(psd.size()>0,"psd vector is empty");
    //check whether the distribution is cumulative
    for (auto it=psd.begin()+1; it!=psd.end(); ++it) {
        logger.assert_always(it->probability >= (it-1)->probability,"psd is not cumulative");
    }
    //psd needs to start with a probability of zero
    if (psd[0].probability!=0) {
//        logger(INFO,"adding a zero at teh beginning of the psd");
        psd.insert(psd.begin(),{psd[0].radius,0});
    }
    //psd needs to end with a probability of one
    if (psd.back().probability<1) {
        logger(INFO,"adding a one at the end of the psd");
        psd.push_back({psd.back().radius,1});
    }
}

// combine vectors of sieve data (radii and probabilities) into subtractive volume distribution
std::vector<PSDContinuous> PSDContinuous::createPSDFromRadiiAndProbabilities (const std::vector<double>& radius,
                                                          const std::vector<double>& probability) {
    logger.assert_always(radius.size()==probability.size(),"Number of radius and probability values don't match");
    std::vector<PSDContinuous> psd;
//    psd.reserve(radius.size());
//    for (int i = 0; i < radius.size(); ++i) {
//        psd.push_back({radius[i],probability[i]});
//    }
    psd.reserve(radius.size()+1);
    psd.push_back({radius.front(),0});
    for (int i = 1; i < radius.size(); ++i) {
        psd.push_back({radius[i],probability[i-1]});
    }
    psd.push_back({radius.back(),probability.back()});
    return psd;
}

// converts a subtractive to a cumulative psd: vc_0 = 0, vc_i=v_i + vc_i-1, vc_end=1
void PSDContinuous::convertSubtractiveToCumulative (std::vector<PSDContinuous>& psd) {
    //add up probabilities
    for (auto it=psd.begin()+1; it!=psd.end(); ++it) {
        it->probability = std::min(1.0, it->probability + (it-1)->probability);
    }
    //check whether cumulative distribution adds up to one
    logger.assert_debug(fabs(psd.back().probability-1)<1e-12,"cumulative distribution needs to add up to one");
    //remove rounding errors
    psd.back().probability=1;
}

// converst a single to a cumulative psd: v_i= vc_i - vc_i-1
void PSDContinuous::convertCumulativeToSubtractive (std::vector<PSDContinuous>& psd) {
    //subtract cumulative probabilities
    //double v = psd[1].probability - psd[0].probability;
    double vcOld=0, vc;
    for (auto it=psd.begin(); it!=psd.end(); ++it) {
        vc = it->probability;
        it->probability -= vcOld;
        vcOld = vc;
    }
}

// converts a subtractive number psd to a subtractive volume psd
void PSDContinuous::convertSubtractiveVolumeToNumber (std::vector<PSDContinuous>& psd) {
    logger.assert_always(psd[0].probability==0,"psd has to start with zero");
    //subtract cumulative probabilities
    double sum=0;
    for (auto it=psd.begin()+1; it!=psd.end(); ++it) {
        ///\todo should this be:
        //         const Mdouble R = it->radius;
        //        const Mdouble r = (it - 1)->radius;
        //        it->probability = R==r ? 0 : it->probability/(R*R*R*R-r*r*r*r);
        it->probability /= 0.25*(square(it->radius)+square((it-1)->radius))*(it->radius+(it-1)->radius);
        sum += it->probability;
    }
    for (auto& p : psd) {
        p.probability /= sum;
    }
}

// converts a cumulative number psd to a cumulative volume psd
void PSDContinuous::convertCumulativeVolumeToNumber (std::vector<PSDContinuous>& psd) {
    convertCumulativeToSubtractive(psd);
    convertSubtractiveVolumeToNumber(psd);
    convertSubtractiveToCumulative(psd);
}

void PSDContinuous::interpolateCSD (std::vector<PSDContinuous>& psd, unsigned n) {
    std::vector<PSDContinuous> psd2;
    for (int i = 0; i < psd.size()-1; ++i) {
        const double r0 = psd[i].radius;
        const double dr = psd[i+1].radius-psd[i].radius;
        const double p0 = psd[i].probability;
        const double dp = psd[i+1].probability-psd[i].probability;
        for (double j = 0; j < n; ++j) {
            psd2.push_back({r0+j/n*dr,p0+j/n*dp});
        }
    }
    psd2.push_back(psd.back());
    psd = psd2;
}


// converts a subtractive number psd to a subtractive volume psd
void PSDContinuous::convertSubtractiveNumberToVolume (std::vector<PSDContinuous>& psd) {
    logger.assert_always(psd[0].probability==0,"psd has to start with zero");
    //subtract cumulative probabilities
    double sum=0;
    for (auto it=psd.begin()+1; it!=psd.end(); ++it) {
        ///\todo should this be:
        //const Mdouble R = it->radius;
        //const Mdouble r = (it - 1)->radius;
        //it->probability *= R*R*R*R-r*r*r*r;
        it->probability *= 0.25 * (square(it->radius) + square((it - 1)->radius)) * (it->radius + (it - 1)->radius);
        sum += it->probability;
    }
    for (auto& p : psd) {
        p.probability /= sum;
    }
}

// converts a cumulative number psd to a cumulative volume psd
void PSDContinuous::convertCumulativeNumberToVolume (std::vector<PSDContinuous>& psd) {
    convertCumulativeToSubtractive(psd);
    convertSubtractiveNumberToVolume(psd);
    convertSubtractiveToCumulative(psd);
}
//// returns a cumulative number particle size distribution from a given volume particle size distribution
//std::vector<PSDContinuous> PSDContinuous::setFromVolumePSD (const std::vector<double>& radius,
//                                                      const std::vector<double>& probability) {
//    logger.assert_always(radius.size()==probability.size(),"Number of radius and probability values don't match");
//    std::vector<PSDContinuous> psd;
//    psd.reserve(radius.size());
//    for (int i = 0; i < radius.size(); ++i) {
//        psd.push_back({radius[i],probability[i]});
//    }
//    return psd;
//    ///\todo why is there no conversion in this function?
//}

// get median size
double PSDContinuous::getQuantile(std::vector<PSDContinuous> psd, double quantile) {
    logger.assert_always(quantile<=1 && quantile>=0,"quantile is not between 0 and 1");
    convertCumulativeNumberToVolume(psd);
    auto high = std::lower_bound(psd.begin(),psd.end(),quantile);
    auto low = std::max(psd.begin(),high-1);
    if (high->probability==low->probability)
        return high->radius;
    else
        return low->radius + (high->radius-low->radius)*(quantile-low->probability)/(high->probability-low->probability);
}


double PSDContinuous::getVolumetricMean(std::vector<PSDContinuous> psd) {
    convertCumulativeToSubtractive(psd);
    double mean =0;
    for (auto it=psd.begin()+1; it!=psd.end(); ++it)
        mean += it->probability * 0.5*(cubic(it->radius)+cubic((it-1)->radius));
    mean = pow(mean,1./3.);
    return mean;
}

const PSD convertPSD2ToPSD(const std::vector<PSDContinuous>& psd2) {
    std::vector<DistributionElements> rp;
    for (const auto p : psd2) {
        rp.push_back({p.radius, p.probability});
    }
    PSD psd;
    psd.setParticleSizeDistribution(rp);
    return psd;
}