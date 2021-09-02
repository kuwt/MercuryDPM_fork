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

#include <DPMBase.h>
#include "PSD.h"

/*!
 * \details Default constructor; sets every data member to 0 or default.
 */
PSD::PSD()
{
    for (auto& p : momenta_)
        p = 0;
}

/*!
 * \details Copy constructor
 */
PSD::PSD(const PSD& other)
{
    particleSizeDistribution_ = other.particleSizeDistribution_;
    momenta_ = other.momenta_;
}

/*!
 * \details Destructor. Since there are no pointers in this class, there is no
 *          need for any actions here.
 */
PSD::~PSD()
= default;


/*!
 * \details Copy method; creates a copy on the heap and returns its pointer.
 * \return pointer to the copy on the heap
 */
PSD* PSD::copy() const
{
#ifdef DEBUG_CONSTRUCTOR
    std::cout << "PSD::copy() const finished" << std::endl;
#endif
    return new PSD(*this);
}

/*!
 * \details Prints the radii [m] and probabilities [%] of the psd vector. It currently supports nanometers, micrometers,
 *          milimeters and meters as size units.
 */
void PSD::printPSD()
{
    if (particleSizeDistribution_.front().radius > 1e-1)
    {
        for (const auto p : particleSizeDistribution_)
        {
            std::cout << p.radius << "m\t" << p.probability * 100 << "\%\n";
        }
    }
    else if (particleSizeDistribution_.front().radius > 1e-4)
    {
        for (const auto p : particleSizeDistribution_)
        {
            std::cout << p.radius * 1000 << "mm\t" << p.probability * 100 << "\%\n";
        }
    }
    else if (particleSizeDistribution_.front().radius > 1e-7)
    {
        for (const auto p : particleSizeDistribution_)
        {
            std::cout << p.radius * 1e6 << "um\t" << p.probability * 100 << "\%\n";
        }
    }
    else
    {
        for (const auto p : particleSizeDistribution_)
        {
            std::cout << p.radius * 1e9 << "nm\t" << p.probability * 100 << "\%\n";
        }
    }
}

/*!
 * \details Draws a sample probability of a real uniform distribution. A random number is
 * generated in range [0,1] and by linear interpolation a suitable radius is returned. This function is only valid
 * for CumulativeNumberDistributions as particles are drawn and not volumes, areas or lengths.
 * \return A Radius for a randomly assigned probability.
 */
Mdouble PSD::drawSample()
{
    // draw a number between 0 and 1, uniformly distributed
    /// \todo TP: We should add a variable seed. Maybe make it definable by the user?
    static std::mt19937 gen(0);
    static std::uniform_real_distribution<Mdouble> dist(0, 1);
    Mdouble prob = dist(gen);
    // find the interval [low,high] of the psd in which the number sits
    auto high = std::lower_bound(particleSizeDistribution_.begin(), particleSizeDistribution_.end(), prob);
    auto low = std::max(particleSizeDistribution_.begin(), high - 1);
    Mdouble rMin = low->radius;
    Mdouble rMax = high->radius;
    Mdouble pMin = low->probability;
    Mdouble pMax = high->probability;
    // interpolate linearly between [low.radius,high.radius] (assumption: CDF is piecewise linear)
    Mdouble a = (prob - pMin) / (pMax - pMin);
    return a * rMin + (1 - a) * rMax;
}

/*!
 * \details Draws a sample probability of a real uniform distribution within a given size class. A random number is
 * generated in the size class range [r_i,r_i+1] and by linear interpolation a suitable radius is returned. from a
 * CUMULATIVE_VOLUME_DISTRIBUTION the volumeAllowed of each class is checked to insert the PSD as accurate as possible.
 * This function is only valid for cumulativeNumberDistributions as particles are drawn and not volumes, areas or
 * lengths. Furthermore this insertion routine is most accurate for non-continuous particle insertion.
 * \param[in] volume                         volume of the geometry to be filled.
 * \return A Radius for a randomly assigned probability of a specific size class.
 */
Mdouble PSD::insertManuallyByVolume(Mdouble volume)
{
    // initialize the particleNumberPerClass vector if empty.
    if (nParticlesPerClass_.empty())
        nParticlesPerClass_.resize(particleSizeDistribution_.size());
    // initialize the particleNumberPerClass vector if empty.
    if (volumePerClass_.empty())
        volumePerClass_.resize(particleSizeDistribution_.size());
    
    for (auto it = particleSizeDistribution_.end() - 1; it != particleSizeDistribution_.begin(); --it)
    {
        static std::mt19937 gen(0);
        static std::uniform_real_distribution<Mdouble> dist(0, 1);
        Mdouble prob = dist(gen);
        // map the probability to the current class interval
        prob = (it - 1)->probability + (it->probability - (it - 1)->probability) * prob;
        // find the interval [low,high] of the psd in which the number sits
        auto high = std::lower_bound(particleSizeDistribution_.begin(), particleSizeDistribution_.end(), prob);
        auto low = std::max(particleSizeDistribution_.begin(), high - 1);
        Mdouble rMin = low->radius;
        Mdouble rMax = high->radius;
        Mdouble pMin = low->probability;
        Mdouble pMax = high->probability;
        // interpolate linearly between [low.radius,high.radius] (assumption: CDF is piecewise linear)
        int index = std::distance(particleSizeDistribution_.begin(), high);
        Mdouble a = (prob - pMin) / (pMax - pMin);
        Mdouble rad = a * rMin + (1 - a) * rMax;
        // Compare inserted volume against volumeAllowed
        Mdouble radVolume = 4.0 / 3.0 * constants::pi * rad * rad * rad;
        convertCumulativeToProbabilityDensity();
        convertProbabilityDensityNumberDistributionToProbabilityDensityVolumeDistribution();
        Mdouble volumeAllowed = particleSizeDistribution_[index].probability * volume;
        convertProbabilityDensityToProbabilityDensityNumberDistribution(TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION);
        convertProbabilityDensityToCumulative();
        if (volumePerClass_[index] > volumeAllowed)
        {
            continue;
        }
        else if (radVolume + volumePerClass_[index] > volumeAllowed)
        {
            Mdouble differenceVolumeLow = -(volumePerClass_[index] - volumeAllowed);
            Mdouble differenceVolumeHigh = radVolume + volumePerClass_[index] - volumeAllowed;
            Mdouble volumeRatio = differenceVolumeLow / differenceVolumeHigh;
            // If volumeRatio > 1 it will insert anyways, because it should be no problem for the distribution
            prob = dist(gen);
            if (prob <= volumeRatio)
            {
                ++nParticlesPerClass_[index];
                volumePerClass_[index] += radVolume;
                return rad;
            }
            else
            {
                continue;
            }
        }
        else
        {
            ++nParticlesPerClass_[index];
            volumePerClass_[index] += radVolume;
            return rad;
        }
    }
    return 0;
}


/*!
 * \details Validates if a distribution is cumulative by first checking if the psd vector is empty and that the
 * scaling of probabilities is in the range [0,1]. Also it checks if each consecutive value is higher than the latter
 * value (i.e. assuming piecewise linear CDF: p_i > p_i-1). Further it replaces probabilities if the CDF does not
 * start with zero or does not end with unity (i.e. p_0=0 and p_end=1)
 */
void PSD::validateCumulativeDistribution()
{
    // ensure interval of probabilities to be [0,1]
    for (auto p = 0; p < particleSizeDistribution_.size(); ++p)
        particleSizeDistribution_[p].probability =
                particleSizeDistribution_[p].probability / particleSizeDistribution_.back().probability;
    // check whether the distribution is cumulative
    for (auto it = particleSizeDistribution_.begin() + 1; it != particleSizeDistribution_.end(); ++it)
    {
        logger.assert_always(it->probability >= (it - 1)->probability, "psd is not cumulative");
    }
    // cdf needs to start with a probability of zero
    if (particleSizeDistribution_[0].probability != 0)
    {
        logger(INFO, "adding a zero at the beginning of the psd");
        particleSizeDistribution_.insert(particleSizeDistribution_.begin(),
                                         {std::max(1e-3 * particleSizeDistribution_[0].radius,
                                                   2 * particleSizeDistribution_[0].radius -
                                                   particleSizeDistribution_[1].radius), 0});
    }
    // cdf needs to end with a probability of one
    if (particleSizeDistribution_.back().probability < 1)
    {
        logger(INFO, "adding a one at the end of the psd");
        particleSizeDistribution_.push_back({particleSizeDistribution_.back().radius, 1});
    }
    // cut off equal subsequent values on the end of the cdf to remove zero size classes.
    for (auto i = particleSizeDistribution_.end() - 1; i != particleSizeDistribution_.begin(); i--)
    {
        if (i->probability == (i - 1)->probability)
        {
            particleSizeDistribution_.erase(particleSizeDistribution_.end() + std::distance(particleSizeDistribution_
                                                                                                    .end(), i));
        }
    }
}

/*!
 * \details Validates if a PDF is valid by first checking if the PDF starts with zero (i.e. p_0=0). Further it checks
 * if the integral is equal to unity (i.e. assuming piecewise constant PDF: sum(p_i)=1).
 */
void PSD::validateProbabilityDensityDistribution()
{
    Mdouble sum = 0;
    // check if the psd starts with zero probability (i.e. p_0=0)
    if (particleSizeDistribution_[0].probability != 0)
    {
        logger(INFO, "adding a zero at the beginning of PDF");
        particleSizeDistribution_.insert(particleSizeDistribution_.begin(),
                                         {std::max(1e-3 * particleSizeDistribution_[0].radius,
                                                   2 * particleSizeDistribution_[0].radius -
                                                   particleSizeDistribution_[1].radius), 0});
    }
    // sum up probabilities to check if it equals to unity (sum(p_i)=1)
    for (auto& it : particleSizeDistribution_)
    {
        sum += it.probability;
    }
    // ensure interval of probabilities to be [0,1] by normalization.
    for (auto& it : particleSizeDistribution_)
        it.probability /= sum;
    // if the sum of probabilities is not equal to unity, sum up the normalized probabilities again
    if (sum - 1 > 1e-6)
    {
        sum = 0;
        for (auto& it : particleSizeDistribution_)
        {
            sum += it.probability;
        }
    }
    logger.assert(sum - 1 < 1e-6, "PDF is not valid: Integral of PDF is not equal to unity");
}

/*!
 * \details creates the psd vector from radii and probabilities filled in by hand. The Type of PSD will be converted
 * to the default cumulative number distribution function (CUMULATIVE_NUMBER_DISTRIBUTION) for further processing.
 * \param[in] psdVector                         Vector containing radii and probabilities ([0,1]).
 * \param[in] PSDType                           Type of the PSD: CUMULATIVE_VOLUME_DISTRIBUTION,
 *                                              CUMULATIVE_NUMBER_DISTRIBUTION, CUMULATIVE_LENGTH_DISTRIBUTION,
 *                                              CUMULATIVE_AREA_DISTRIBUTION, PROBABILITYDENSITY_VOLUME_DISTRIBUTION,
 *                                              PROBABILITYDENSITY_NUMBER_DISTRIBUTION,
 *                                              PROBABILITYDENSITY_LENGTH_DISTRIBUTION,
 *                                              PROBABILITYDENSITY_AREA_DISTRIBUTION.
 * \deprecated This is the old way of inserting PSDs. In the future use setPSDFromCSV().
 */
void PSD::setPSDFromVector(std::vector<RadiusAndProbability> psdVector, TYPE PSDType)
{
    particleSizeDistribution_ = psdVector;
    switch (PSDType)
    {
        // default case is a cumulative number distribution (CND)
        default:
            validateCumulativeDistribution();
            break;
        case TYPE::CUMULATIVE_VOLUME_DISTRIBUTION:
            validateCumulativeDistribution();
            convertCumulativeToProbabilityDensity();
            convertProbabilityDensityToProbabilityDensityNumberDistribution(
                    TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION);
            convertProbabilityDensityToCumulative();
            break;
        case TYPE::CUMULATIVE_LENGTH_DISTRIBUTION:
            validateCumulativeDistribution();
            convertCumulativeToProbabilityDensity();
            convertProbabilityDensityToProbabilityDensityNumberDistribution(
                    TYPE::PROBABILITYDENSITY_LENGTH_DISTRIBUTION);
            convertProbabilityDensityToCumulative();
            break;
        case TYPE::CUMULATIVE_AREA_DISTRIBUTION:
            validateCumulativeDistribution();
            convertCumulativeToProbabilityDensity();
            convertProbabilityDensityToProbabilityDensityNumberDistribution(TYPE::PROBABILITYDENSITY_AREA_DISTRIBUTION);
            convertProbabilityDensityToCumulative();
            break;
        case TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION:
            validateProbabilityDensityDistribution();
            convertProbabilityDensityToProbabilityDensityNumberDistribution(PSDType);
            convertProbabilityDensityToCumulative();
            break;
        case TYPE::PROBABILITYDENSITY_NUMBER_DISTRIBUTION:
            validateProbabilityDensityDistribution();
            convertProbabilityDensityToCumulative();
            break;
        case TYPE::PROBABILITYDENSITY_LENGTH_DISTRIBUTION:
            validateProbabilityDensityDistribution();
            convertProbabilityDensityToProbabilityDensityNumberDistribution(PSDType);
            convertProbabilityDensityToCumulative();
            break;
        case TYPE::PROBABILITYDENSITY_AREA_DISTRIBUTION:
            validateProbabilityDensityDistribution();
            convertProbabilityDensityToProbabilityDensityNumberDistribution(PSDType);
            convertProbabilityDensityToCumulative();
            break;
    }
}

/*!
 * \details creates the PSD vector from probabilities and radii saved in a .csv file. Radii should be located at
 *          the first column and probabilities at the second column of the file. The Type of PSD will be converted to the
 *          default cumulative number distribution function (CUMULATIVE_NUMBER_DISTRIBUTION) for further processing.
 * \param[in] fileName                          Name of the .csv file containing the radii and probabilities of the PSD.
 * \param[in] PSDType                           Type of the PSD: CUMULATIVE_VOLUME_DISTRIBUTION,
 *                                              CUMULATIVE_NUMBER_DISTRIBUTION, CUMULATIVE_LENGTH_DISTRIBUTION,
 *                                              CUMULATIVE_AREA_DISTRIBUTION, PROBABILITYDENSITY_VOLUME_DISTRIBUTION,
 *                                              PROBABILITYDENSITY_NUMBER_DISTRIBUTION,
 *                                              PROBABILITYDENSITY_LENGTH_DISTRIBUTION PROBABILITYDENSITY_AREA_DISTRIBUTION.
 * \param[in] headings                          If TRUE the file is assumed to have headings and the first row will
 *                                              be skipped. If FALSE the file has no headings and the file will be
 *                                              read in as is. Default is FALSE.
 * \param[in] unitScalingFactorRadii            Scaling factor of radii to match SI-units.
 */
void PSD::setPSDFromCSV(const std::string& fileName, TYPE PSDType, bool headings, Mdouble unitScalingFactorRadii)
{
    // read in csv file using the csvReader class
    csvReader csv;
    csv.setHeader(headings);
    csv.read(fileName);
    std::vector<Mdouble> radii = csv.getFirstColumn(unitScalingFactorRadii);
    std::vector<Mdouble> probabilities = csv.getSecondColumn(1.0);
    logger.assert_always(radii.size() == probabilities.size(), "The radii and probabilities vector have to be the "
                                                               "same size");
    // combine radii and probabilities
    for (int i = 0; i < radii.size(); ++i)
    {
        particleSizeDistribution_.push_back({radii[i], probabilities[i]});
    }
    logger.assert_always(!particleSizeDistribution_.empty(), "PSD cannot be empty");
    switch (PSDType)
    {
        // default case is a CUMULATIVE_NUMBER_DISTRIBUTION
        default:
            validateCumulativeDistribution();
            break;
        case TYPE::CUMULATIVE_LENGTH_DISTRIBUTION:
            validateCumulativeDistribution();
            convertCumulativeToProbabilityDensity();
            convertProbabilityDensityToProbabilityDensityNumberDistribution(
                    TYPE::PROBABILITYDENSITY_LENGTH_DISTRIBUTION);
            convertProbabilityDensityToCumulative();
            break;
        case TYPE::CUMULATIVE_AREA_DISTRIBUTION:
            validateCumulativeDistribution();
            convertCumulativeToProbabilityDensity();
            convertProbabilityDensityToProbabilityDensityNumberDistribution(TYPE::PROBABILITYDENSITY_AREA_DISTRIBUTION);
            convertProbabilityDensityToCumulative();
            break;
        case TYPE::CUMULATIVE_VOLUME_DISTRIBUTION:
            validateCumulativeDistribution();
            convertCumulativeToProbabilityDensity();
            convertProbabilityDensityToProbabilityDensityNumberDistribution(
                    TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION);
            convertProbabilityDensityToCumulative();
            break;
        case TYPE::PROBABILITYDENSITY_NUMBER_DISTRIBUTION:
            validateProbabilityDensityDistribution();
            convertProbabilityDensityToCumulative();
            break;
        case TYPE::PROBABILITYDENSITY_LENGTH_DISTRIBUTION:
            validateProbabilityDensityDistribution();
            convertProbabilityDensityToProbabilityDensityNumberDistribution(PSDType);
            convertProbabilityDensityToCumulative();
            break;
        case TYPE::PROBABILITYDENSITY_AREA_DISTRIBUTION:
            validateProbabilityDensityDistribution();
            convertProbabilityDensityToProbabilityDensityNumberDistribution(PSDType);
            convertProbabilityDensityToCumulative();
            break;
        case TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION:
            validateProbabilityDensityDistribution();
            convertProbabilityDensityToProbabilityDensityNumberDistribution(PSDType);
            convertProbabilityDensityToCumulative();
            break;
    }
}

/*!
 * \details Convert any type of PDF to a CDF. Probabilities are integrated for each radius by cumulatively summing
 * them up (i.e. p_i = p_i + p_i-1).
 */
void PSD::convertProbabilityDensityToCumulative()
{
    // add up probabilities
    for (auto it = particleSizeDistribution_.begin() + 1; it != particleSizeDistribution_.end(); ++it)
    {
        it->probability = std::min(1.0, it->probability + (it - 1)->probability);
    }
    // check whether cumulative distribution is valid
    validateCumulativeDistribution();
}

/*!
 * \details Convert any type of CDF to a PDF. Probabilities are derivated for each radius by substracting the CDF
 * probabilities (i.e. pPDF_i = pCDF_i - pCDF_i-1).
 */
void PSD::convertCumulativeToProbabilityDensity()
{
    // subtract cumulative probabilities
    Mdouble probabilityOld = 0, probabilityCDF;
    for (auto& it : particleSizeDistribution_)
    {
        probabilityCDF = it.probability;
        it.probability -= probabilityOld;
        probabilityOld = probabilityCDF;
    }
    // check whether probability density distribution is valid
    validateProbabilityDensityDistribution();
}

/*!
 * \details converts any of the PDFTypes to a PROBABILITYDENSITY_NUMBER_DISTRIBUTION based on their TYPE. This is a helper function
 * which enables the convertProbabilityDensityToCumulative function to convert the psd into the default TYPE (CUMULATIVE_NUMBER_DISTRIBUTION).
 * \param[in] PDFType                            Type of the PDF: PROBABILITYDENSITY_LENGTH_DISTRIBUTION, PROBABILITYDENSITY_AREA_DISTRIBUTION or PROBABILITYDENSITY_VOLUME_DISTRIBUTION,
 *                                               Where L = Length, A = Area and V = Volume.
 */
void PSD::convertProbabilityDensityToProbabilityDensityNumberDistribution(TYPE PDFType)
{
    Mdouble sum = 0;
    switch (PDFType)
    {
        default:
            logger(ERROR, "Wrong PDFType");
            break;
        case TYPE::PROBABILITYDENSITY_LENGTH_DISTRIBUTION:
            for (auto it = particleSizeDistribution_.begin() + 1; it != particleSizeDistribution_.end(); ++it)
            {
                it->probability /= 0.5 * (it->radius + (it - 1)->radius);
                // old conversion
                //p.probability /= p.radius;
                // sum up probabilities
                sum += it->probability;
            }
            // normalize
            for (auto& p : particleSizeDistribution_)
            {
                p.probability /= sum;
            }
            break;
        case TYPE::PROBABILITYDENSITY_AREA_DISTRIBUTION:
            for (auto it = particleSizeDistribution_.begin() + 1; it != particleSizeDistribution_.end(); ++it)
            {
                it->probability /=
                        1.0 / 3.0 * (square(it->radius) + it->radius * (it - 1)->radius + square((it - 1)->radius));
                // old conversion
                //p.probability /= square(p.radius);
                // sum up probabilities
                sum += it->probability;
            }
            // normalize
            for (auto& p : particleSizeDistribution_)
            {
                p.probability /= sum;
            }
            break;
        case TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION:
            for (auto it = particleSizeDistribution_.begin() + 1; it != particleSizeDistribution_.end(); ++it)
            {
                it->probability /=
                        0.25 * (square(it->radius) + square((it - 1)->radius)) * (it->radius + (it - 1)->radius);
                // old conversion
                //p.probability /= cubic(p.radius);
                // sum up probabilities
                sum += it->probability;
            }
            // normalize
            for (auto& p : particleSizeDistribution_)
            {
                p.probability /= sum;
            }
            break;
    }
}

/*!
 * \details converts a PROBABILITYDENSITY_NUMBER_DISTRIBUTION to a PROBABILITYDENSITY_VOLUME_DISTRIBUTION. Used for
 * getVolumetricMeanRadius() and insertManuallyByVolume().
 */
void PSD::convertProbabilityDensityNumberDistributionToProbabilityDensityVolumeDistribution()
{
    Mdouble sum = 0;
    for (auto it = particleSizeDistribution_.begin() + 1; it != particleSizeDistribution_.end(); ++it)
    {
        it->probability *=
                0.25 * (square(it->radius) + square((it - 1)->radius)) * (it->radius + (it - 1)->radius);
        // old conversion
        //p.probability /= cubic(p.radius);
        // sum up probabilities
        sum += it->probability;
    }
    // normalize
    for (auto& p : particleSizeDistribution_)
    {
        p.probability /= sum;
    }
}

/*!
 * \details converts any of the CDFTypes to a the default CUMULATIVE_NUMBER_DISTRIBUTION based on their TYPE.
 * \param[in] PDFType                            Type of the PDF: PROBABILITYDENSITY_LENGTH_DISTRIBUTION, PROBABILITYDENSITY_AREA_DISTRIBUTION or PROBABILITYDENSITY_VOLUME_DISTRIBUTION,
 *                                               Where L = Length, A = Area and V = Volume.
 */
/// \todo TP: NOT WORKING! If anyone knows how to do it feel free to add
void PSD::convertCumulativeToCumulativeNumberDistribution(TYPE CDFType)
{
    Mdouble sum = 0;
    switch (CDFType)
    {
        default:
            logger(ERROR, "Wrong CDFType");
            break;
        case TYPE::CUMULATIVE_LENGTH_DISTRIBUTION:
            for (auto it = particleSizeDistribution_.begin() + 1; it != particleSizeDistribution_.end(); ++it)
            {
                // add conversion here
            
                // sum up probabilities
                sum += it->probability;
            }
            // normalize
            for (auto& p : particleSizeDistribution_)
            {
                p.probability /= sum;
            }
            break;
        case TYPE::CUMULATIVE_AREA_DISTRIBUTION:
            for (auto it = particleSizeDistribution_.begin() + 1; it != particleSizeDistribution_.end(); ++it)
            {
                // add conversion here
            
                // sum up probabilities
                sum += it->probability;
            }
            // normalize
            for (auto& p : particleSizeDistribution_)
            {
                p.probability /= sum;
            }
            break;
        case TYPE::CUMULATIVE_VOLUME_DISTRIBUTION:
            for (auto it = particleSizeDistribution_.begin() + 1; it != particleSizeDistribution_.end(); ++it)
            {
                // add conversion here
            
                // sum up probabilities
                sum += it->probability;
            }
            // normalize
            for (auto& p : particleSizeDistribution_)
            {
                p.probability /= sum;
            }
            break;
    }
}

/*!
 * \details Cuts off the CDF at given minimum and maximum percentiles and applies a minimum polydispersity at the base.
 * \param[in] percentileMin                       undersize percentile to cut off the lower part of the CDF.
 * \param[in] percentileMax                       oversize percentile to cut off the upper part of the CDF.
 * \param[in] minPolydispersity                 Applies a minimum of polydispersity ([0,1]) at the base of the CDF.
 */
void PSD::cutoffCumulativeNumber(Mdouble percentileMin, Mdouble percentileMax, Mdouble minPolydispersity)
{
    Mdouble radiusMin = getRadiusByPercentile(percentileMin);
    Mdouble radiusMax = getRadiusByPercentile(percentileMax);
    // to get a minimum polydispersity at the base
    Mdouble radiusMinCut = std::min(radiusMin * (1 + minPolydispersity), radiusMax);
    // cut off min
    while (particleSizeDistribution_.front().radius <= radiusMinCut)
        particleSizeDistribution_.erase(particleSizeDistribution_.begin());
    particleSizeDistribution_.insert(particleSizeDistribution_.begin(),
                                     particleSizeDistribution_[radiusMinCut, percentileMin]);
    particleSizeDistribution_.insert(particleSizeDistribution_.begin(), particleSizeDistribution_[radiusMin, 0]);
    // cut off max
    while (particleSizeDistribution_.back().radius >= radiusMax) particleSizeDistribution_.pop_back();
    particleSizeDistribution_.push_back(particleSizeDistribution_[radiusMax, percentileMax]);
    particleSizeDistribution_.push_back(particleSizeDistribution_[radiusMax, 1]);
}

/*!
 * \details Cuts off the CDF at given minimum and maximum percentiles, applies a minimum polydispersity at the base and
 * squeezes the distribution to make it less polydisperse.
 * \param[in] percentileMin                       undersize percentile to cut off the lower part of the CDF.
 * \param[in] percentileMax                       oversize percentile to cut off the upper part of the CDF.
 * \param[in] squeeze                           applies a squeezing factor ([0,1]) which determines the degree the
 *                                              PDF gets squeezed.
 * \param[in] minPolydispersity                 applies a minimum of polydispersity ([0,1]) at the base of the CDF.
 */
void PSD::cutoffAndSqueezeCumulative(Mdouble percentileMin, Mdouble percentileMax, Mdouble squeeze, Mdouble
minPolydispersity)
{
    Mdouble r50 = 0.5 * PSD::getDx(0.5);
    // cut off
    cutoffCumulativeNumber(percentileMin, percentileMax, minPolydispersity);
    // squeeze psd
    for (auto& p : particleSizeDistribution_)
    {
        p.radius = r50 + (p.radius - r50) * squeeze;
    }
}

/*!
 * \details Gets the diameter from a certain percentile of the PSD.
 * \return A double which is the diameter corresponding to a certain percentile.
 * \param[in] x                 double which determines the obtained diameter as a percentile of the PSD.
 */
Mdouble PSD::getDx(Mdouble x)
{
    return 2.0 * getRadiusByPercentile(x / 100);
}

/*!
 * \details gets the radius from a certain percentile of the PSD
 * \return A double which is the radius corresponding to a certain percentile of the PSD.
 * \param[in] percentile                 double which determines the returned radius as a percentile of the PSD.
 */
Mdouble PSD::getRadiusByPercentile(Mdouble percentile)
{
    logger.assert_always(percentile <= 1 && percentile >= 0, "percentile is not between 0 and 1");
    // find the percentile corresponding to the PSD
    auto high = std::lower_bound(particleSizeDistribution_.begin(), particleSizeDistribution_.end(), percentile);
    auto low = std::max(particleSizeDistribution_.begin(), high - 1);
    if (high->probability == low->probability)
        return high->radius;
    else
        return low->radius + (high->radius - low->radius) * (percentile - low->probability) /
                             (high->probability - low->probability);
}

/*!
 * \details Gets a radius such that a monodisperse system has the same number of particles as a polydisperse system.
 * (i.e. mean += p_i * 0.5*(r_i^3 + r_i-1^3)
 * \return A double which corresponds to the volumetric mean radius.
 */
Mdouble PSD::getVolumetricMeanRadius()
{
    Mdouble mean = 0;
    convertCumulativeToProbabilityDensity();
    convertProbabilityDensityNumberDistributionToProbabilityDensityVolumeDistribution();
    for (auto it = particleSizeDistribution_.begin() + 1; it != particleSizeDistribution_.end(); ++it)
        mean += it->probability * 0.5 * (cubic(it->radius) + cubic((it - 1)->radius));
    mean = pow(mean, 1. / 3.);
    convertProbabilityDensityToProbabilityDensityNumberDistribution(TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION);
    convertProbabilityDensityToCumulative();
    return mean;
}

/*!
 * \details Gets the minimal radius of the PSD.
 * \return A double which corresponds to the minimal radius of the PSD.
 */
Mdouble PSD::getMinRadius()
{
    return particleSizeDistribution_[0].radius;
}

/*!
 * \details Gets the maximum radius of the PSD.
 * \return A double which corresponds to the maximum radius of the PSD.
 */
Mdouble PSD::getMaxRadius()
{
    return particleSizeDistribution_.back().radius;
}

/*!
 * \details Gets the vector containing radii and probabilities of the PSD.
 * \return A vector containing radii and probabilities of the PSD.
 */
const std::vector<PSD::RadiusAndProbability> PSD::getParticleSizeDistribution() const
{
    return particleSizeDistribution_;
}

/*!
 * \details Gets the number of particles already inserted into the simulation by summing up the particles inserted in
 * each class.
 * \return An integer containing the number of particles already inserted into the simulation.
 */
int PSD::getInsertedParticleNumber()
{
    int sum = 0;
    for (auto& it : nParticlesPerClass_)
        sum += it;
    return sum;
}

/*!
 * \details Compute the raw momenta of the inserted PSD by converting it to a PDF, calculating the moments m_i according
 * to \f$ m_i = \sum_{j=0}^n \f$ scaling it by the number
 * of particles inserted into the simulation (moments are computed from a PROBABILITYDENSITY_NUMBER_DISTRIBUTION) and
 * converting it back to a CDF.
 * See <a href="https://en.wikipedia.org/wiki/Moment_(mathematics)#Significance_of_the_moments">Wikipedia</a> for
 * details.
 */
void PSD::computeRawMomenta()
{
    convertCumulativeToProbabilityDensity();
    for (size_t im = 0; im < momenta_.size(); ++im)
    {
        // prevent summing up of moments for each time step of the simulation
        momenta_[im] = 0;
        for (auto& it : particleSizeDistribution_)
        {
            momenta_[im] += std::pow(it.radius, im) * it.probability;
        }
    }
    // zeroth-moment equals particle number for a PROBABILITYDENSITY_NUMBER_DISTRIBUTION
    momenta_[0] *= getInsertedParticleNumber();
    convertProbabilityDensityToCumulative();
}

/*!
 * \details Compute the central momenta of the PSD from their respective raw momenta.
 */
void PSD::computeCentralMomenta()
{
    computeRawMomenta();
    Mdouble mean = momenta_[1];
    momenta_[5] += -5 * mean * momenta_[4] + 10 * mean * mean * momenta_[3]
                   - 10 * mean * mean * mean * momenta_[2] + 4 * mean * mean * mean * mean * mean;
    momenta_[4] += -4 * mean * momenta_[3] + 6 * mean * mean * momenta_[2] - 3 * mean * mean * mean * mean;
    momenta_[3] += -3 * mean * momenta_[2] + 2 * mean * mean * mean;
    momenta_[2] += -mean * mean;
}

/*!
 * \details Compute the standardised momenta of the PSD from their respective central momenta.
 */
void PSD::computeStandardisedMomenta()
{
    computeCentralMomenta();
    Mdouble std = std::sqrt(momenta_[2]);
    momenta_[3] /= std * std * std;
    momenta_[4] /= std * std * std * std;
    momenta_[5] /= std * std * std * std * std;
}

/*!
 * \details Get momenta of the user defined particleSizeDistribution_ vector.
 * \return Array of momenta corresponding to the first six moments of the user defined PSD.
 */
std::array<Mdouble, 6> PSD::getMomenta()
{
    return momenta_;
}

/*!
 * \details Required to use std::lower_bound for finding when the probability is higher than a certain value.
 * \return TRUE if probability from a vector of type PSD::RadiusAndProbability is higher than a certain value from a vector of type
 * PSD::RadiusAndProbability and FALSE in the opposite case.
 */
bool operator<(const PSD::RadiusAndProbability& l, const PSD::RadiusAndProbability& r)
{
    return l.probability < r.probability;
}

/*!
 * \details required to use std::lower_bound for finding when the probability provided as a double is higher than a
 * certain value.
 * \return TRUE if probability as double is higher than a certain value from a PSD::RadiusAndProbability vector and FALSE in the opposite
 * case.
 */
bool operator<(const PSD::RadiusAndProbability& l, const Mdouble prob)
{
    return l.probability < prob;
}

/*!
 * \details Required to use std::distance to find the index of the PSD size class in which a particle has to be inserted
 * \return A double which determines the size class (radius) a particle will be inserted to.
 */
Mdouble operator==(PSD::RadiusAndProbability l, const Mdouble r)
{
    return l.radius == r;
}

/*!
 * \details reads from input stream. This function is used for restart files.
 * \return a reference to an input stream.
 */
std::istream& operator>>(std::istream& is, PSD::RadiusAndProbability& p)
{
    is >> p.radius;
    is >> p.probability;
    return is;
}

/*!
 * \details Writes to output stream. This function is used for restart files.
 * \return a reference to an output stream.
 */
std::ostream& operator<<(std::ostream& os, PSD::RadiusAndProbability& p)
{
    os << p.radius << ' ' << p.probability << ' ';
    return os;
}