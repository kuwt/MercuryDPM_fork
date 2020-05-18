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
#ifndef PANEL_H_
#define PANEL_H_

#include <vector>
#include "Dipole.h"
#include "LocalExpansion.h"
#include "Math/NumericalVector.h"
#include "Math/Vector.h"
#include "Multipole.h"
#include "Source.h"

class Panel;

#include "Sphere.h"

class Box;

class Panel
{
public:
    Panel(Panel* root,
          int maximumPanelLevel,
          Vec3D leftBound,
          Vec3D rightBound,
          std::vector<Source*> sources,
          std::vector<Dipole*> dipoles,
          NumericalVector<>* squaredFactorials,
          Box* box);
    
    //Data structure functions
    void initialise();
    
    void computeCoefficients();
    
    void createPanels(int dim, std::vector<Source*>& sources, std::vector<Dipole*>& dipoles, Vec3D& leftBoundChild,
                      Vec3D& rightBoundChild, NumericalVector<>* squaredFactorials);
    
    void findPanelInteractions();
    
    void setPanelInteractions();
    
    //Functions used for the upward pass
    void computeMultipoleExpansion();
    
    void translateMultipoleExpansion();
    
    //Functions used for the downward pass
    void setLocalExpansionZero(); //In the current implementation this function is not used
    void computePartialLocalExpansion();
    
    void translateLocalExpansion();
    
    
    //Functions that need to be removed or shift to another topic
    void computeLocalExpansion();
    
    
    //Getters
    Vec3D getCentre()
    {
        return centre_;
    }
    
    Panel* getRoot()
    {
        return root_;
    }
    
    std::vector<Panel*> getChilderen()
    {
        return childeren_;
    }
    
    std::vector<Panel*> getNeighbours()
    {
        return neighbours_;
    }
    
    std::vector<Panel*> getSecondNeighbours()
    {
        return secondNeighbours_;
    }
    
    std::vector<Panel*> getInteractionList()
    {
        return interactionList_;
    }
    
    std::vector<Source*> getSources()
    {
        return sources_;
    }
    
    Source* getSource(int index)
    {
        return sources_[index];
    }
    
    std::vector<Dipole*> getDipoles()
    {
        return dipoles_;
    }
    
    std::vector<Multipole*> getMultipoles()
    {
        return multipoles_;
    }
    
    std::vector<Sphere*> getSpheres()
    {
        return spheres_;
    }
    
    int getPanelLevel()
    {
        return panelLevel_;
    }
    
    NumericalVector<std::complex<Mdouble>> getPartialLocalExpansion()
    {
        return partialLocalExpansion_;
    }
    
    NumericalVector<std::complex<Mdouble>> getLocalExpansion()
    {
        return localExpansion_;
    }

private:
    //Panel characteristics
    const int panelLevel_;    //The level at which the panel is living. This panel has no children when panelLevel = 1
    Mdouble dim_;            //dimension of the problem space
    Vec3D leftBound_;        // Left and bottom bounds
    Vec3D rightBound_;        // Right and top bounds
    double size_;            // the half of the width and length of the panel
    Vec3D centre_;            // Centre of the panel
    
    //Data structures
    Panel* root_ = nullptr;                    //Root of the current panel
    std::vector<Panel*> childeren_;            //List of children which have the current panel as root.
    std::vector<Panel*> neighbours_;        //List of neighbours
    std::vector<Panel*> secondNeighbours_;  //List of second neighbours.
    std::vector<Panel*> interactionList_;    //List of interaction panels
    
    //Multipole structures
    std::vector<Source*> sources_;            //List of sources contained within this panel.
    std::vector<Dipole*> dipoles_;            //List of dipoles contained within this panel.
    std::vector<Multipole*> multipoles_;    //List of multipoles contained within this panel.
    
    //For fluid calculation
    std::vector<Sphere*> spheres_;            //For the finest level of panels this vector contains the spheres within that panel
    
    // Mother of all pointers
    Box* box_;                                //Pointer to the data structure sorted per level
    
    //Computational values
    // todo: sort out this mess
    Multipole* multipoleAroundCentre_;
    LocalExpansion* partialLocalExpansionAroundCentre_;
    LocalExpansion* localExpansionAroundCentre_;
    NumericalVector<std::complex<Mdouble>> partialLocalExpansion_;
    NumericalVector<std::complex<Mdouble>> localExpansion_;
};


#endif /* PANEL_H_ */
