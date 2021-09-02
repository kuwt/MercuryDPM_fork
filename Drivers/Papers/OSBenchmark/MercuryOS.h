//Copyright (c) 2013-2021, The MercuryDPM Developers Team. All rights reserved.
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
#ifndef MERCURY_MERCURYOS_H
#define MERCURY_MERCURYOS_H

#include <Mercury3D.h>
#include <Species/HertzianViscoelasticMindlinSpecies.h>
#include <Logger.h>

/**
 * This class contains the code shared between is various OS*.cpp drivers.
 * In particular, it provides a function to define the material properties and provides two booleans that determine (a) whether to write a full set of output files and (b) whether to use smooth walls instead of the triangulated walls
 */
class MercuryOS : public Mercury3D
{
    // if true, additional output files are written
    bool writeOutput_ = false;
    // if true, the simulation time is shortened
    bool test_ = false;
    // if true, smooth walls are used instead of the stl walls
    bool useMercuryWalls_ = false;
    // if true, the simulation time is shortened
    bool soft_ = false;
    
public:
    
    // sets the variable writeOutput_
    void writeOutput(bool writeOutput)
    {
        writeOutput_ = writeOutput;
    }
    
    // returns the variable writeOutput_
    bool writeOutput() const
    {
        return writeOutput_;
    }

    // sets the variable test_
    void test(bool test)
    {
        test_ = test;
    }
    
    // returns the variable test_
    bool test() const
    {
        return test_;
    }
    
    // sets the variable soft_
    void soft(bool soft)
    {
        soft_ = soft;
    }
    
    // returns the variable soft_
    bool soft() const
    {
        return soft_;
    }
    
    // sets the variable useMercuryWalls_
    void useMercuryWalls(bool useMercuryWalls)
    {
        useMercuryWalls_ = useMercuryWalls;
    }
    
    // returns the variable useMercuryWalls_
    bool useMercuryWalls() const
    {
        return useMercuryWalls_;
    }

private:
    
    /**
     * Computes the effective elastic modulus from the Young's modulus and Poisson ratio of the two particles.
     */
    static double getEffectiveElasticModulus(double E1, double v1, double E2, double v2)
    {
        return 1. / ((1 - v1 * v1) / E1 + (1 - v2 * v2) / E2);
    }
    
    /**
     * Computes the effective shear modulus from the Young's modulus and Poisson ratio of the two particles.
     */
    static double getEffectiveShearModulus(double E1, double v1, double E2, double v2)
    {
        double G1 = E1 / 2 / (1 + v1), G2 = E2 / 2 / (1 + v2);
        return 1. / ((2 - v1) / G1 + (2 - v2) / G2);
    }
    
protected:
    // pointers to the species
    HertzianViscoelasticMindlinSpecies *m1, *m2, *steel;

public:
    /**
     * Defines the material properties of M1, M2, steel.
     */
    void setMaterialProperties()
    {
        // Young's modulus and Poisson ratios of M1, M2, steel
        double E1 = 1e9, E2 = 0.5e9, ES = 210e9;
        double v1 = 0.2, v2 = 0.2, vS = 0.2;
        
        if (soft()) {
            E1 = 2e6;
            E2 = 1e6;
        }
        
        // effective elastic and shear moduli
        double E11 = getEffectiveElasticModulus(E1, v1, E1, v1);
        double E12 = getEffectiveElasticModulus(E1, v1, E2, v2);
        double E1S = getEffectiveElasticModulus(E1, v1, ES, vS);
        double E22 = getEffectiveElasticModulus(E2, v2, E2, v2);
        double E2S = getEffectiveElasticModulus(E2, v2, ES, vS);
        double ESS = getEffectiveElasticModulus(ES, vS, ES, vS);
        double G11 = getEffectiveShearModulus(E1, v1, E1, v1);
        double G12 = getEffectiveShearModulus(E1, v1, E2, v2);
        double G1S = getEffectiveShearModulus(E1, v1, ES, vS);
        double G22 = getEffectiveShearModulus(E2, v2, E2, v2);
        double G2S = getEffectiveShearModulus(E2, v2, ES, vS);
        double GSS = getEffectiveShearModulus(ES, vS, ES, vS);
        
        // restitution coefficients
        double r11 = 0.5, r12 = 0.45, r1S = 0.4, r22 = 0.4, r2S = 0.4, rSS = 0.6;
        // sliding friction coefficients
        double mu11 = 0.3, mu12 = 0.2, mu1S = 0.2, mu22 = 0.4, mu2S = 0.2, muSS = 0.5;
        
        HertzianViscoelasticMindlinSpecies species;
        species.setEffectiveElasticModulusAndRestitutionCoefficient(E11, r11);
        species.setEffectiveShearModulus(G11);
        species.setSlidingFrictionCoefficient(mu11);
        species.setDensity(2500);
        m1 = speciesHandler.copyAndAddObject(species);
    
        species.setEffectiveElasticModulusAndRestitutionCoefficient(E22, r22);
        species.setEffectiveShearModulus(G22);
        species.setSlidingFrictionCoefficient(mu22);
        species.setDensity(2000);
        m2 = speciesHandler.copyAndAddObject(species);
    
        species.setEffectiveElasticModulusAndRestitutionCoefficient(ESS, rSS);
        species.setEffectiveShearModulus(GSS);
        species.setSlidingFrictionCoefficient(muSS);
        species.setDensity(7200);
        steel = speciesHandler.copyAndAddObject(species);
        
        auto m12 = speciesHandler.getMixedObject(m1, m2);
        m12->setEffectiveElasticModulusAndRestitutionCoefficient(E12, r12);
        m12->setEffectiveShearModulus(G12);
        m12->setSlidingFrictionCoefficient(mu12);
        
        auto m1S = speciesHandler.getMixedObject(m1, steel);
        m1S->setEffectiveElasticModulusAndRestitutionCoefficient(E1S, r1S);
        m1S->setEffectiveShearModulus(G1S);
        m1S->setSlidingFrictionCoefficient(mu1S);
        
        auto m2S = speciesHandler.getMixedObject(m2, steel);
        m2S->setEffectiveElasticModulusAndRestitutionCoefficient(E2S, r2S);
        m2S->setEffectiveShearModulus(G2S);
        m2S->setSlidingFrictionCoefficient(mu2S);
    }
    
    
};


#endif //MERCURY_MERCURYOS_H
