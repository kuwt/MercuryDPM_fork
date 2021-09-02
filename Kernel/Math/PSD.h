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

#ifndef PSD_H
#define PSD_H

#include <fstream>
#include <Logger.h>
#include <vector>
#include <Tools/csvReader.h>
#include <random>
#include "iostream"
#include "ExtendedMath.h"

using mathsFunc::square;
using mathsFunc::cubic;

/*!
 * \brief Contains a vector with radii and probabilities of a user defined particle size distribution (PSD)
 *
 * \details Stores radii and probabilities of a particle size distribution (PSD) in a vector of type PSD::RadiusAndProbability and
 * converts them to a PSD which can be used by insertionBoundaries which insert particles into a simulation:
 *
 * Cumulative distribution function (CDF): gives percentage of particles p_i whose radius r is less than
 * a certain radius r_i. It requires p_0 = 0 and p_end = 1. The CDF's probabilities can be number, length, area or
 * volume based. The cumulative number distribution function (CUMULATIVE_NUMBER_DISTRIBUTION) is used as PSD in MercuryDPM, as it can be
 * interpreted as the probability that a particle's radius is less than r_i:
 *     CDF(r<r_i) = p_i.
 *
 * Probability density function (PDF): p_i is the percentage of particles whose radius is between r_i-1 and
 * r_i. It requires p_0=0 and sum(p_i)=1. The PDF's probabilities can also be number, length, area or volume based.
 * PDF's are utilized to convert any type of PDF to a probability number density function (PROBABILITYDENSITY_NUMBER_DISTRIBUTION) which in a next step
 * are converted to the default CUMULATIVE_NUMBER_DISTRIBUTION.
 *
 * Sieve data: p_i is the percentage of particles whose radius is between r_i and r_i+1. It requires sum(p_i)=1.
 * relation to PDF: p_i = p(PDF)_i-1. Sieve data is not yet used in this class.
 *
 * Default distribution is the CUMULATIVE_NUMBER_DISTRIBUTION.
 */
class PSD
{

public:
    /*!
     * \brief Enum class which stores the possible types of CDFs and PDFs.
     * Particle size distributions can be represented by different probabilities based on number of particles, length
     * of particles, surface area of particles or volume of particles.
     */
    enum class TYPE
    {
        CUMULATIVE_NUMBER_DISTRIBUTION,
        CUMULATIVE_LENGTH_DISTRIBUTION,
        CUMULATIVE_VOLUME_DISTRIBUTION,
        CUMULATIVE_AREA_DISTRIBUTION,
        PROBABILITYDENSITY_NUMBER_DISTRIBUTION,
        PROBABILITYDENSITY_LENGTH_DISTRIBUTION,
        PROBABILITYDENSITY_AREA_DISTRIBUTION,
        PROBABILITYDENSITY_VOLUME_DISTRIBUTION
    };
    
    /*!
     * \brief Class which stores radii and probabilities of a PSD. This class should be used as a vector<PSD::RadiusAndProbability>.
     */
    class RadiusAndProbability
    {
    public:
        Mdouble radius;
        Mdouble probability;
    };
    
    /*!
     * \brief Constructor; sets everything to 0 or default.
     */
    PSD();
    
    /*!
     * \brief Copy constructor with deep copy.
     */
    PSD(const PSD& other);
    
    /*!
     * \brief Destructor; default destructor.
     */
    ~PSD();
    
    /*!
     * \brief Creates a copy on the heap and returns a pointer.
     */
    PSD* copy() const;
    
    /*!
     * \brief Prints radii and probabilities of the PSD vector.
     */
    void printPSD();
    
    /*!
     * \brief Draw a sample radius from a CUMULATIVE_NUMBER_DISTRIBUTION.
     */
    Mdouble drawSample();
    
    /*!
     * \brief Draw sample radius manually per size class and check the volumeAllowed of each size class to insert the
     * PSD as accurate as possible.
     */
    Mdouble insertManuallyByVolume(Mdouble volume);
    
    /*!
     * \brief Validates if a CDF starts with zero and adds up to unity.
    */
    void validateCumulativeDistribution();
    
    /*!
     * \brief Validates if the integral of the PDF equals to unity.
     */
    void validateProbabilityDensityDistribution();
    
    /*!
     * \brief Deprecated version of reading in PSDs from a vector.
     */
    MERCURY_DEPRECATED
    void setPSDFromVector(std::vector<RadiusAndProbability> psd, TYPE PSDType);
    
    /*!
     * \brief read in the PSD vector with probabilities and radii saved in a .csv file.
     */
    void setPSDFromCSV(const std::string& fileName, TYPE PSDType, bool headings = false, Mdouble
    unitScalingFactorRadii = 1.0);
    
    /*!
     * \brief Converts a PDF to a CDF by integration.
     */
    void convertProbabilityDensityToCumulative();
    
    /*!
     * \brief Converts a CDF to a PDF by derivation.
     */
    void convertCumulativeToProbabilityDensity();
    
    /*!
     * \brief convert any PDF to a PROBABILITYDENSITY_NUMBER_DISTRIBUTION.
     */
    void convertProbabilityDensityToProbabilityDensityNumberDistribution(TYPE PDFType);
    
    /*!
     * \brief convert a PROBABILITYDENSITY_NUMBER_DISTRIBUTION to a PROBABILITYDENSITY_VOLUME_DISTRIBUTION.
     */
    void convertProbabilityDensityNumberDistributionToProbabilityDensityVolumeDistribution();
    
    /*!
     * \brief convert any other CDF to a CUMULATIVE_NUMBER_DISTRIBUTION.
     */
    void convertCumulativeToCumulativeNumberDistribution(TYPE CDFType);
    
    /*!
     * \brief cutoff the PSD at given percentiles.
     */
    void cutoffCumulativeNumber(Mdouble percentileMin, Mdouble percentileMax, Mdouble minPolydispersity = 0.1);
    
    /*!
     * \brief cutoff the PSD at given percentiles and make it less polydisperse by squeezing it.
     */
    void cutoffAndSqueezeCumulative(Mdouble percentileMin, Mdouble percentileMax, Mdouble squeeze,
                                    Mdouble minPolydispersity = 0.1);
    
    /*!
     * \brief Get smallest radius of the PSD.
     */
    Mdouble getMinRadius();
    
    /*!
     * \brief Get largest radius of the PSD.
     */
    Mdouble getMaxRadius();
    
    /*!
     * \brief Get the PSD vector.
     */
    const std::vector<RadiusAndProbability> getParticleSizeDistribution() const;
    
    /*!
     * \brief Get the number of particles already inserted into the simulation.
     */
    int getInsertedParticleNumber();
    
    /*!
     * \brief Calculate a certain diameter (e.g. D10, D50, D90, etc.) from a percentile of the PSD.
     */
    Mdouble getDx(Mdouble x);
    
    /*!
     * \brief Calculate the percentile of the PSD.
     */
    Mdouble getRadiusByPercentile(Mdouble percentile);
    
    /*!
     * \brief get a volumetric mean radius of the PSD.
     */
    Mdouble getVolumetricMeanRadius();
    
    /*!
     * \brief compute raw momenta of the user defined PSD.
     */
    void computeRawMomenta();
    
    /*!
     * \brief compute central momenta of the user defined PSD.
     */
    void computeCentralMomenta();
    
    /*!
     * \brief compute standardised momenta of the user defined PSD.
     */
    void computeStandardisedMomenta();
    
    /*!
     * \brief get momenta of the user defined PSD.
     */
    std::array<Mdouble, 6> getMomenta();
    
    /*!
     * \brief determines if a certain value of the PSD vector is lower than another one. Used for std::lower_bound()
     */
    friend bool operator<(const PSD::RadiusAndProbability& l, const PSD::RadiusAndProbability& r);
    
    /*!
     * \brief determines if a certain value of the PSD vector is lower than a double.
     */
    friend bool operator<(const PSD::RadiusAndProbability& l, Mdouble r);
    
    /*!
     * \brief Writes to output stream.
     */
    friend std::ostream& operator<<(std::ostream& os, PSD::RadiusAndProbability& p);
    
    /*!
     * \brief Reads from input stream.
     */
    friend std::istream& operator>>(std::istream& is, PSD::RadiusAndProbability& p);
    
    /*!
     * \brief Determines if a certain value of the PSD vector is equal to a double.
     */
    friend Mdouble operator==(PSD::RadiusAndProbability l, Mdouble r);


private:
    /*!
     * Vector of the PSD::RadiusAndProbability class which stores radii and probabilities of the PSD.
     */
    std::vector<RadiusAndProbability> particleSizeDistribution_;
    
    /*!
     * Array of doubles which stores the moments of a user defined discrete PROBABILITYDENSITY_NUMBER_DISTRIBUTION.
     */
    std::array<Mdouble, 6> momenta_{};
    
    /*!
     * Vector of integers which represents the number of inserted particles in each size class.
     * The classes in this vector are defined to contain the particles between size r_i and r_i-1. (e.g. size class
     * 12 consists of particles between size class 12 and 11 of the PDF)
     */
    std::vector<int> nParticlesPerClass_;
    
    /*!
     * Vector of doubles which stores the volume of inserted particles for each size class. This vector is used in
     * the insertManuallyByVolume() function to check if the volumeAllowed per class is exceeded and thus no further
     * particles should be added to a certain class.
     * The classes in this vector are defined to contain the volume of particles between size r_i and r_i-1. (e.g. size
     * class 12 consists of the particles' volume between size class 12 and 11 of the PDF)
     */
    std::vector<Mdouble> volumePerClass_;
};


#endif //PSD_H
