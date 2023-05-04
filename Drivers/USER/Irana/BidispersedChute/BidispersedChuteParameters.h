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

#ifndef MERCURYDPM_BIDISPERSEDCHUTEPARAMETERS_H
#define MERCURYDPM_BIDISPERSEDCHUTEPARAMETERS_H


class BidispersedChuteParameters
{
public:
    BidispersedChuteParameters()
    {
        inflowHeight = 10;
        angle = 24;
        concentrationSmall = 0.5;
        largeParticleRadius = 0.5;
        smallParticleRadius = 0.5;
    }

    /*!
     * Constructor; sets the variable properties of the simulation
     * @param h the estimated height; determines the volume of particles that will be inserted into the domain
     * @param a inclincation of the periodic box
     * @param phi concentration of small particles
     * @param s ratio of large-to-small particle diameter
     */
    BidispersedChuteParameters(Mdouble h, Mdouble a, Mdouble phi) :
            inflowHeight(h), angle(a), concentrationSmall(phi)
    {
        logger.assert_always(h > -1e-10, "inflow height should be positive");
        logger.assert_always(a > -1e-10 && a < 90, "angle should be between 0 and 90 degrees");
        logger.assert_always(phi > -1e-10 && (phi - 1) < 1e-10, "concentration of small particles must be between "
                "0 and 1");
        if (a < constants::pi / 2)
            logger(WARN, "Angle should be in degrees");
        largeParticleRadius = 0.5;
        smallParticleRadius = 0.5/std::cbrt(2);
        fixedParticleRadius = 0.5;
        logger(DEBUG, "large particle radius: %, small particle radius: %", getLargeParticleRadius(), getSmallParticleRadius());
    }
    
    void setLargeParticleRadius(Mdouble radius)
    {
        largeParticleRadius = radius;
    }
    
    void setSmallParticleRadius(Mdouble radius)
    {
        smallParticleRadius = radius;
    }
    
    Mdouble getAngleInDegrees() const
    {
        return angle;
    }
    
    Mdouble getLargeParticleRadius() const
    {
        return largeParticleRadius;
    }
    
    //choose d^S s.t. (d^L phi^L + d^S phi^S) = 1
    Mdouble getSmallParticleRadius() const
    {
        return smallParticleRadius;
    }
    
    Mdouble getFixedParticleRadius() const
    {
        return fixedParticleRadius;
    }
    
    void setFixedParticleRadius(Mdouble radius)
    {
        fixedParticleRadius = radius;
    }

    /*!
     * Returns the number of large particles in terms of the surface area (Dx*Dy) of the domain
     */
    unsigned int getNumberOfLargeParticles(Mdouble area) const
    {
        logger.assert_always(area > 0, "chute should have positive area");
        const Mdouble totalVolume = area * inflowHeight;
        return (1-concentrationSmall) * totalVolume / std::pow(2*getLargeParticleRadius(), 3.0);
        const Mdouble volumeForLarge = (1 - concentrationSmall) * totalVolume;
        const Mdouble volumeOneLarge = 4./3 * constants::pi * std::pow(getLargeParticleRadius(), 3.0);
        const Mdouble packingFraction = 0.57;
        return static_cast<unsigned int>(volumeForLarge / volumeOneLarge * packingFraction);
    }

    unsigned int getNumberOfLargeParticlesExtraAndUpdate(Mdouble area, Mdouble height)
    {
        logger.assert_always(area > 0 && height > 0, "extra volume should be positive");
        const unsigned int currentNumberOfLarge = getNumberOfLargeParticles(area);
        inflowHeight += height * (1-concentrationSmall);
        const unsigned int targetNumberOfLarge = getNumberOfLargeParticles(area);
        return targetNumberOfLarge - currentNumberOfLarge;
    }
    
    unsigned int getNumberOfSmallParticles(Mdouble area) const
    {
        logger.assert_always(area > 0, "chute should have positive area");
        const Mdouble totalVolume = area * inflowHeight;
        return concentrationSmall * totalVolume / std::pow(2*getSmallParticleRadius(), 3.0);
        const Mdouble volumeForSmall = concentrationSmall * totalVolume;
        const Mdouble volumeOneSmall = 4./3 * constants::pi * std::pow(getSmallParticleRadius(), 3.0);
        const Mdouble packingFraction = 0.57;
        return static_cast<unsigned int>(volumeForSmall / volumeOneSmall * packingFraction);
    }
    
    unsigned int getNumberOfSmallParticlesExtraAndUpdate(Mdouble area, Mdouble height)
    {
        logger.assert_always(area > 0 && height > 0, "extra volume should be positive");
        const unsigned int currentNumberOfSmall = getNumberOfSmallParticles(area);
        inflowHeight += height * concentrationSmall;
        const unsigned int targetNumberOfSmall = getNumberOfSmallParticles(area);
        return targetNumberOfSmall - currentNumberOfSmall;
    }
    
    Mdouble getInflowHeight() const
    {
        return inflowHeight;
    }
    
    void setAngle(Mdouble newAngle)
    {
        angle = newAngle;
    }
    
    Mdouble getConcentrationSmall() const
    {
        return concentrationSmall;
    }

private:
    //total height of the particles, in large particle diameters
    Mdouble inflowHeight;
    //angle in degrees
    Mdouble angle;
    //concentration of small particles, in [0,1]
    Mdouble concentrationSmall;
    //radius of large particles
    Mdouble largeParticleRadius;
    //radius of small particles
    Mdouble smallParticleRadius;
    Mdouble fixedParticleRadius;
};


#endif //MERCURYDPM_BIDISPERSEDCHUTEPARAMETERS_H
