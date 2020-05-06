#include "PSD.h"
#include <Logger.h>
#include "ExtendedMath.h"
using mathsFunc::square;
using mathsFunc::cubic;

// validate whether psd is a cumulative distribution
void PSD::validateCumulativeDistribution (std::vector<PSD>& psd) {
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
std::vector<PSD> PSD::createPSDFromRadiiAndProbabilities (const std::vector<double>& radius,
        const std::vector<double>& probability) {
    logger.assert_always(radius.size()==probability.size(),"Number of radius and probability values don't match");
    std::vector<PSD> psd;
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
void PSD::convertSubtractiveToCumulative (std::vector<PSD>& psd) {
    //add up probabilities
    for (auto it=psd.begin()+1; it!=psd.end(); ++it) {
        it->probability = std::min(1.0, it->probability + (it-1)->probability);
    }
    //check whether cumulative distribution adds up to one
    logger.assert(fabs(psd.back().probability-1)<1e-12,"cumulative distribution needs to add up to one");
    //remove rounding errors
    psd.back().probability=1;
}

// converst a single to a cumulative psd: v_i= vc_i - vc_i-1
void PSD::convertCumulativeToSubtractive (std::vector<PSD>& psd) {
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
void PSD::convertSubtractiveVolumeToNumber (std::vector<PSD>& psd) {
    logger.assert_always(psd[0].probability==0,"psd has to start with zero");
    //subtract cumulative probabilities
    double sum=0;
    for (auto it=psd.begin()+1; it!=psd.end(); ++it) {
//        it->probability /= 0.25*(pow(it->radius,4)-pow((it-1)->radius,4))/(it->radius-(it-1)->radius);
        it->probability /= 0.25*(square(it->radius)+square((it-1)->radius))*(it->radius+(it-1)->radius);
        //p.probability /= cubic(p.radius);
        sum += it->probability;
    }
    for (auto& p : psd) {
        p.probability /= sum;
    }
}

// converts a cumulative number psd to a cumulative volume psd
void PSD::convertCumulativeVolumeToNumber (std::vector<PSD>& psd) {
    convertCumulativeToSubtractive(psd);
    convertSubtractiveVolumeToNumber(psd);
    convertSubtractiveToCumulative(psd);
}

// converts a subtractive number psd to a subtractive volume psd
void PSD::convertSubtractiveNumberToVolume (std::vector<PSD>& psd) {
    logger.assert_always(psd[0].probability==0,"psd has to start with zero");
    //subtract cumulative probabilities
    double sum=0;
    for (auto it=psd.begin()+1; it!=psd.end(); ++it) {
        //it->probability *= 0.25 * (pow(it->radius, 4) - pow((it - 1)->radius, 4)) / (it->radius - (it - 1)->radius);
        it->probability *= 0.25 * (square(it->radius) + square((it - 1)->radius)) * (it->radius + (it - 1)->radius);
        //p.probability /= cubic(p.radius);
        sum += it->probability;
    }
    for (auto& p : psd) {
        p.probability /= sum;
    }
}

// converts a cumulative number psd to a cumulative volume psd
void PSD::convertCumulativeNumberToVolume (std::vector<PSD>& psd) {
    convertCumulativeToSubtractive(psd);
    convertSubtractiveNumberToVolume(psd);
    convertSubtractiveToCumulative(psd);
}
// returns a cumulative number particle size distribution from a given volume particle size distribution
std::vector<PSD> PSD::setFromVolumePSD (const std::vector<double>& radius,
                                                      const std::vector<double>& probability) {
    logger.assert_always(radius.size()==probability.size(),"Number of radius and probability values don't match");
    std::vector<PSD> psd;
    psd.reserve(radius.size());
    for (int i = 0; i < radius.size(); ++i) {
        psd.push_back({radius[i],probability[i]});
    }
    return psd;
}

// get median size
double PSD::getQuantile(std::vector<PSD> psd, double quantile) {
    logger.assert_always(quantile<=1 && quantile>=0,"quantile is not between 0 and 1");
    convertCumulativeNumberToVolume(psd);
    auto high = std::lower_bound(psd.begin(),psd.end(),quantile);
    auto low = std::max(psd.begin(),high-1);
    return low->radius + (high->radius-low->radius)*(quantile-low->probability)/(high->probability-low->probability);
}


double PSD::getVolumetricMean(std::vector<PSD> psd) {
    convertCumulativeToSubtractive(psd);
    double mean =0;
    for (auto it=psd.begin()+1; it!=psd.end(); ++it)
        mean += it->probability * 0.5*(cubic(it->radius)+cubic((it-1)->radius));
    mean = pow(mean,1./3.);
    return mean;
}
