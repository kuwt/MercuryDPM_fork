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
#include "Box.h"
#include "Math/Vector.h"
#include "Panel.h"
#include "Sphere.h"

Box::Box(int maxLevel, int nTerms) : maxLevel_(maxLevel), p_(nTerms)
{
    for (int iL = 0; iL <= maxLevel; iL++)
    {
        std::vector<Panel*>* boxLevel = new std::vector<Panel*>;
        levels_.push_back(*boxLevel);
    }
}

void Box::addPanel(int level, Panel* panel)
{
    levels_[level].push_back(panel);
}

void Box::upwardPass()
{
    std::cout << "============================" << std::endl;
    std::cout << "Starting the upward pass." << std::endl;
    // Perform a multipole expansion about the centre of all panels located on the finest mesh
    std::cout << "Computing multipole expansion on finest level." << std::endl;
    for (Panel* panel : levels_[maxLevel_])
    {
        panel->computeMultipoleExpansion();
    }
    
    // Perform an upward pass, each level combine the previous level multipole expansions by shifting them to the centre of the current panel
    for (int iL = maxLevel_ - 1; iL >= 0; iL--)
    {
        std::cout << "Shifting multipole expansions on level " << iL << "." << std::endl;
        for (Panel* panel : levels_[iL])
        {
            panel->translateMultipoleExpansion();
        }
    }
    std::cout << "Finished upward pass." << std::endl;
    
}

void Box::downwardPass()
{
    std::cout << "============================" << std::endl;
    std::cout << "Starting the downward pass." << std::endl;
    
    //Set the first part of the local expansion to zero: everything is in the interaction list or closer
    NumericalVector<std::complex<Mdouble>> localExpansionCoefficients;
    for (Panel* panel : levels_[1])
    {
        panel->setLocalExpansionZero();
    }
    
    // Compute all local expansions for intermediate levels
    for (int iL = 1; iL < maxLevel_; iL++)
    {
        std::cout << "Computing nearby Local Expansions on level " << iL << "." << std::endl;
        for (Panel* panel : levels_[iL])
        {
            panel->computePartialLocalExpansion();
            panel->computeLocalExpansion();
            
        }
        
        std::cout << "Translating local expansions to children" << std::endl;
        for (Panel* panel : levels_[iL])
        {
            panel->translateLocalExpansion();
        }
    }
    
    // Compute local expansions for the remaining  interaction list on the finest level
    std::cout << "Computing nearby local expansions on finest level " << std::endl;
    for (Panel* panel : levels_[maxLevel_])
    {
        panel->computePartialLocalExpansion();
        panel->computeLocalExpansion();
    }
    
    std::cout << "Finished downward pass." << std::endl;
}

void Box::computeFlow(int k)
{
    //***************************************
    //* 		Initialise spheres			*
    //***************************************
    
    // The first step is to add dipoles with correct location and strengths in the domain
    // All dipoles in the domain are now expanded into multipoles
    std::cout << "Initialising spheres" << std::endl;
    upwardPass();
    downwardPass();
    
    //For all panels on the finest level
    for (Panel* panel : levels_[maxLevel_])
    {
        //For all dipoles on the finest level
        for (Dipole* iD : panel->getDipoles())
        {
            // Add a multipole in the domain
            //Multipole* multipole = new Multipole(iD->getP, iD->getSquaredFactorials, iD->getLocation());
            //panel->multipoles_.push_back(multipole);
            
            // Construct a sphere
            //Sphere* sphere = new Sphere(panel, iD->getLocation(), iD, multipole);
            
            // Add the sphere to the list of spheres and to the panel list
            //spheres_.push_back(sphere);
            //spheres_.push_back(sphere);
        }
    }
    
    //***********************************
    //* 	Compue Mirror Multipoles	*
    //***********************************
    
    
    // Compute k order mirrors
    for (int i = 1; k <= i; i++)
    {
        //For all panels on the finest level
        for (Panel* panel : levels_[maxLevel_])
        {
            // For all spheres: add a new multipole to the centre of a sphere to compensate
            for (Sphere* sphere : panel->getSpheres())
            {
                std::vector<std::complex<Mdouble>> localExpansionAroundSphere;
                Vec3D sphereCentre = sphere->getLocation();
                
                // Compute panel local expansion around sphere
                //localExpansionAroundSphere = panel->localExpansionAroundCentre_->translateLocalExpansion(sphereCentre);
                
                // Compute other sphere local expansion around sphere
                for (Sphere* sphereOther: spheres_)
                {
                    if (sphereOther != sphere)
                    {
                        // todo: implement this shizzle
                        //localExpansionAroundSphere += sphereOther->multipole_->convertMultipoleToLocal(sphereCentre);
                    }
                }
                // Compute multipole
                size_t nTerms = (p_ + 1) * (p_ + 1);
                NumericalVector<std::complex<Mdouble>> multipoleExpansionCoefficients(nTerms);
                for (int n = 0; n <= p_; n++)
                {
                    for (int m = -n; m <= n; m++)
                    {
                        int location = n * n + (m + n);
                        //multipoleExpansionCoefficients[location] = localExpansionAroundSphere[0]*n/(n+1)*std::pow(sphere->getRadius,2n+1);
                    }
                }
                //sphere->multipole_->setExpansionCoefficients(multipoleExpansionCoefficients);
            }
        }
        
        //Perform and upward and downward pass to compute the new coefficients
        upwardPass();
        downwardPass();
    }
    
    
}








