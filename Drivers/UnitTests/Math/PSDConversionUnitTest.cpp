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

#include "Math/PSD.h"

std::string getPSDTypeName(PSD::TYPE psdType)
{
    if (psdType == PSD::TYPE::CUMULATIVE_NUMBER_DISTRIBUTION)
        return "CUMULATIVE_NUMBER_DISTRIBUTION";
    else if (psdType == PSD::TYPE::CUMULATIVE_LENGTH_DISTRIBUTION)
        return "CUMULATIVE_LENGTH_DISTRIBUTION";
    else if (psdType == PSD::TYPE::CUMULATIVE_AREA_DISTRIBUTION)
        return "CUMULATIVE_AREA_DISTRIBUTION";
    else if (psdType == PSD::TYPE::CUMULATIVE_VOLUME_DISTRIBUTION)
        return "CUMULATIVE_VOLUME_DISTRIBUTION";
    else if (psdType == PSD::TYPE::PROBABILITYDENSITY_NUMBER_DISTRIBUTION)
        return "PROBABILITYDENSITY_NUMBER_DISTRIBUTION";
    else if (psdType == PSD::TYPE::PROBABILITYDENSITY_LENGTH_DISTRIBUTION)
        return "PROBABILITYDENSITY_LENGTH_DISTRIBUTION";
    else if (psdType == PSD::TYPE::PROBABILITYDENSITY_AREA_DISTRIBUTION)
        return "PROBABILITYDENSITY_AREA_DISTRIBUTION";
    else if (psdType == PSD::TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION)
        return "PROBABILITYDENSITY_VOLUME_DISTRIBUTION";
    else
        logger(ERROR, "Unknown PSD type");
    return "";
}

bool isEqual(const std::vector<DistributionElements>& l, const std::vector<DistributionElements>& r)
{
    if (l.size() != r.size())
        return false;

    Mdouble tol = 1.0e-15;
    for (int i = 0; i < l.size(); i++)
        if (std::fabs(l[i].internalVariable - r[i].internalVariable) > tol || std::fabs(l[i].probability - r[i].probability) > tol)
            return false;

    return true;
}

/*!
 * \brief Sets a PSD by vector with a given type and gets the PSD vector back by the same type to confirm it's the same.
 * \details The PSD class internally converts any PSD type to CUMULATIVE_NUMBER_DISTRIBUTION. It is also possible to get
 * the PSD vector by any type, where internally it is converted back from CUMULATIVE_NUMBER_DISTRIBUTION. Doing this for
 * the same type should of course give the same PSD vector back, provided the PSD vector was perfect to start with, so
 * that internal validation didn't alter it in any way.
 * \note This merely checks if converting to and converting back from CUMULATIVE_NUMBER_DISTRIBUTION will yield the same
 * PSD vector. It does not check for the correctness of the conversions itself. It is therefore possible that forwards
 * and backwards conversion both contain the same but opposite mistake. This test would then still pass.
 * @param psdType The PSD type the test.
 */
void testPSD(PSD::TYPE psdType)
{
    // Set a default PSD vector, differentiating between cumulative and probability density.
    // Note: it is important that this is a proper normalised PSD vector, so that the PSD class does not alter it in any
    // way when validating it. Because then converting back would result in a different PSD vector.
    std::vector<DistributionElements> psdVector;
    switch (psdType)
    {
        case PSD::TYPE::CUMULATIVE_NUMBER_DISTRIBUTION:
        case PSD::TYPE::CUMULATIVE_LENGTH_DISTRIBUTION:
        case PSD::TYPE::CUMULATIVE_AREA_DISTRIBUTION:
        case PSD::TYPE::CUMULATIVE_VOLUME_DISTRIBUTION:
            psdVector = { {1.0, 0.0}, {2.0, 0.5}, {3.0, 0.6}, {4.0, 1.0} };
            break;
        case PSD::TYPE::PROBABILITYDENSITY_NUMBER_DISTRIBUTION:
        case PSD::TYPE::PROBABILITYDENSITY_LENGTH_DISTRIBUTION:
        case PSD::TYPE::PROBABILITYDENSITY_AREA_DISTRIBUTION:
        case PSD::TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION:
            psdVector = { {1.0, 0.0}, {2.0, 0.5}, {3.0, 0.1}, {4.0, 0.4} };
            break;
        default:
            logger(ERROR, "Unknown PSD type");
    }

    PSD psd;
    // Set the PSD by the correct type (internally converted to CUMULATIVE_NUMBER_DISTRIBUTION).
    psd.setPSDFromVector(psdVector, psdType);
    // Get the PSD vector by the same type (internally converted from CUMULATIVE_NUMBER_DISTRIBUTION).
    std::vector<DistributionElements> psdVectorBack = psd.getParticleSizeDistributionByType(psdType);

    // The original and converted back PSD vectors should be equal.
    logger.assert_always(isEqual(psdVector, psdVectorBack), "Original and converted back PSD vectors of PSD type % are not equal!", getPSDTypeName(psdType));
}

void runTest()
{
    testPSD(PSD::TYPE::CUMULATIVE_NUMBER_DISTRIBUTION);
    testPSD(PSD::TYPE::CUMULATIVE_LENGTH_DISTRIBUTION);
    testPSD(PSD::TYPE::CUMULATIVE_AREA_DISTRIBUTION);
    testPSD(PSD::TYPE::CUMULATIVE_VOLUME_DISTRIBUTION);
    testPSD(PSD::TYPE::PROBABILITYDENSITY_NUMBER_DISTRIBUTION);
    testPSD(PSD::TYPE::PROBABILITYDENSITY_LENGTH_DISTRIBUTION);
    testPSD(PSD::TYPE::PROBABILITYDENSITY_AREA_DISTRIBUTION);
    testPSD(PSD::TYPE::PROBABILITYDENSITY_VOLUME_DISTRIBUTION);
}

int main()
{
    runTest();
}