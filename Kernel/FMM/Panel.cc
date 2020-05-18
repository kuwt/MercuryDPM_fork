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
#include <cstddef>
#include <cmath>
#include <vector>
#include "Dipole.h"
#include "Panel.h"
#include "Math/NumericalVector.h"
#include "Math/Vector.h"
#include "Multipole.h"
#include "Source.h"
#include "Sphere.h"
#include "Box.h"

Panel::Panel(Panel* root,
             int panelLevel,
             Vec3D leftBound,
             Vec3D rightBound,
             std::vector<Source*> sources,
             std::vector<Dipole*> dipoles,
             NumericalVector<>* squaredFactorials,
             Box* box) :
        panelLevel_(panelLevel),
        leftBound_(leftBound),
        rightBound_(rightBound),
        root_(root),
        sources_(sources),
        dipoles_(dipoles),
        box_(box)
{
    //Initialise the panel characteristics
    // I am using Vec3D, so the size is always constant
/*	if (leftBound.size() != rightBound.size())
	{
		std::cout << "Bounds are not of correct dimensions" << std::endl;
		std::exit(-1);
	}*/
    
    // This code only works for 3D problems, since the maths is significantly different from 2D and 1D
    // The datastructure however allows for 1D and 2D structures,
    // In future 2D support might be added (by a master student?)
    dim_ = 3;
    
    // The panel is a cube with same size in all directions
    // Note this is half the length of a cube side
    size_ = 0.5 * (rightBound.getComponent(0) - leftBound.getComponent(0));
    
    // Determine the centre
    for (int iD = 0; iD < dim_; iD++)
    {
        centre_.setComponent(iD, (leftBound.getComponent(iD) + size_));
    }
    
    //Add panel to the complete box
    box->addPanel(panelLevel, this);
    
    //Set multipole and local expansion coefficients
    multipoleAroundCentre_ = new Multipole(box_->getNumberOfTerms(), squaredFactorials, centre_);
    multipoleAroundCentre_->computeMultipoleExpansion(); // Set vector with zeros to correct length
    partialLocalExpansionAroundCentre_ = new LocalExpansion(box_->getNumberOfTerms(), squaredFactorials, centre_);
    partialLocalExpansionAroundCentre_->initialiseLocalExpansion(); // Set vector with zeros to correct length
    localExpansionAroundCentre_ = new LocalExpansion(box_->getNumberOfTerms(), squaredFactorials, centre_);
    localExpansionAroundCentre_->initialiseLocalExpansion(); // Set vector with zeros to correct length
    
    //Create childeren based on the given dimensions;
    if (panelLevel < box->getMaxLevel()) // final level is obtained if this value is 0
    {
        //Note: For an adaptive grid, add a check to see if this panel contains enough sources to have childeren.
        
        Vec3D leftBoundChild;
        Vec3D rightBoundChild;
        
        createPanels(1, sources, dipoles, leftBoundChild, rightBoundChild, squaredFactorials);
    }
    
}

void Panel::createPanels(int dim,
                         std::vector<Source*>& sources,
                         std::vector<Dipole*>& dipoles,
                         Vec3D& leftBoundChild,
                         Vec3D& rightBoundChild,
                         NumericalVector<>* squaredFactorials)
{
    std::vector<Source*> sourcesChild;
    std::vector<Dipole*> dipolesChild;
    
    
    for (int iSide = 0; iSide < 2; ++iSide)
    {
        //Determine the x boundaries
        leftBoundChild.setComponent(dim - 1, (leftBound_.getComponent(dim - 1) + size_ * iSide));
        rightBoundChild.setComponent(dim - 1, (leftBoundChild.getComponent(dim - 1) + size_));
        
        //Determine the sources in this domain
        for (Source* iS : sources)
        {
            if ((iS->getLocation().getComponent(dim - 1) > leftBoundChild.getComponent(dim - 1)) &&
                (iS->getLocation().getComponent(dim - 1) <= rightBoundChild.getComponent(dim - 1)))
            {
                sourcesChild.push_back(iS);
            }
        }
        
        //Determine the dipoles in this domain
        for (Dipole* iD : dipoles)
        {
            if ((iD->getLocation().getComponent(dim - 1) > leftBoundChild.getComponent(dim - 1)) &&
                (iD->getLocation().getComponent(dim - 1) <= rightBoundChild.getComponent(dim - 1)))
            {
                dipolesChild.push_back(iD);
            }
        }
        
        if (dim != dim_)
        {
            createPanels(dim + 1, sourcesChild, dipolesChild, leftBoundChild, rightBoundChild, squaredFactorials);
        }
        else
        {
            //Create childeren
/*			std::cout << "Creating a child with bounds: " << std::endl;
			std::cout << "In x-direction: " << std::endl;
			std::cout << "Left bound: " << leftBoundChild[0] << " and right bound: " << rightBoundChild[0] << std::endl;
			std::cout << "In y-direction: " << std::endl;
			std::cout << "Left bound: " << leftBoundChild[1] << " and right bound: " << rightBoundChild[1] << std::endl;
			std::cout << "In z-direction: " << std::endl;
			std::cout << "Left bound: " << leftBoundChild[2] << " and right bound: " << rightBoundChild[2] << std::endl;
			std::cout << "Child contains " << sourcesChild.size() << " sources" << std::endl;
			std::cout << "==============" << std::endl;*/
            Panel* child = new Panel(this, panelLevel_ + 1, leftBoundChild, rightBoundChild, sourcesChild, dipolesChild,
                                     squaredFactorials, box_);
            childeren_.push_back(child);
            
            //reset leftBoundChild and rightBoundChild
        }
        sourcesChild.clear();
    }
}

void Panel::initialise()
{
    //Create neighbour list and interaction list
    findPanelInteractions();
    
}

void Panel::computeCoefficients()
{
    //Start upward pass creating multipole methods and shifting them on parent panels
    box_->upwardPass();
    
    //Start the downward pass
    box_->downwardPass();
}

void Panel::findPanelInteractions()
{
    if (getRoot() == nullptr)
    {
        //nothing happens
    }
    else
    {
        //Crate the links between panels
        this->setPanelInteractions();
    }
    
    //Initialise child panels
    for (Panel* iP : childeren_)
    {
        iP->findPanelInteractions();
    }
}

void Panel::setPanelInteractions()
{
    //Childeren of a parent are always a neighbour of eachother
    //Note: by definition, a panel is always a neighbour of itself
    for (Panel* iC : root_->getChilderen())
    {
        neighbours_.push_back(iC);
    }
    
    //For all the neighbours of the parent of this panel
    for (Panel* iN : root_->getNeighbours())
    {
        //For all their childeren, check if the childeren are a neighbour or a second neighbour of this panel
        // Note that all second neighbour panels are also childeren of these paretn neighbours
        for (Panel* iC : iN->getChilderen())
        {
            //Compute distance between the centre of the current panel and possible neighbour panel
            Mdouble distance = Vec3D::getDistance(iC->getCentre(), this->getCentre());
            
            //Check if the possible neighbour panel is close enough
            if (distance < 3.0 * size_) //warning: slight hack. I assume all sides are equally long
            {
                neighbours_.push_back(iC);
                //std::cout << "A neighbour has been created with a centre around: (" << iC->getCentre()[0] << "," << iC->getCentre()[1] << ")" << std::endl;
                //std::cout << "Current panel has centre: (" << this->getCentre()[0] << "," << this->getCentre()[1] << ")" << std::endl;
                //std::cout << "===============" << std::endl;
            }
            else if (distance < 6.0 * size_)
            {
                secondNeighbours_.push_back(iC);
            }
            else
            {
                // If it is not a neighbour or nearest neighbour, it belongs to the interactionlist group
                interactionList_.push_back(iC);
                //std::cout << "A panel is added to the interactionlist with a centre around: (" << iC->getCentre()[0] << "," << iC->getCentre()[1] << ")" << std::endl;
                //std::cout << "Current panel has centre: (" << this->getCentre()[0] << "," << this->getCentre()[1] << ")" << std::endl;
                //std::cout << "===============" << std::endl;
            }
        }
    }
    
    // All childeren of the second nearest neighbour of the current panel's parent, are in the interactionlist
    for (Panel* iN : root_->getSecondNeighbours())
    {
        for (Panel* iC : iN->getChilderen())
        {
            secondNeighbours_.push_back(iC);
        }
    }
}

void Panel::computeMultipoleExpansion()
{
    
    /// Source: perform a multipole expansion around the centre of the box
    // todo: write this part of the code
    // A source can be expanded around the centre of the panel without translations
    //add(iS->computeMultipoleExpansion();)
    
    /// Dipole: perform a multipole expansion around the centre of the box
    for (Dipole* iD : dipoles_)
    {
        //Convert to a multipole
        //todo: This must be done when the actual dipole is made.
        iD->computeMultipoleExpansion();
        
        //Translate multipole to centre of the panel
        multipoleAroundCentre_->addMultipoleCoefficients(iD->TranslateMultipoleExpansionTo(centre_));
        
    }
    
    /// Multipole: Transfer a multipole in the box to the centre of the box
    for (Multipole* iM : multipoles_)
    {
        multipoleAroundCentre_->addMultipoleCoefficients(iM->TranslateMultipoleExpansionTo(centre_));
    }
    
}

void Panel::translateMultipoleExpansion()
{
    // Extract the multipole of all the childeren
    for (Panel* iC : childeren_)
    {
        //Translate the multipole in the center of the child to the centre of the current panel
        //And then add it to the multipole in the centre of the current panel
        multipoleAroundCentre_->addMultipoleCoefficients(
                iC->multipoleAroundCentre_->TranslateMultipoleExpansionTo(centre_)
        );
    }
}

void Panel::setLocalExpansionZero()
{
    //Note: nothing has to be done here, the current implementation initialises all local expansions with 0.
}

void Panel::computePartialLocalExpansion()
{
    //For all panels in the interactionlist, add a local expansion of the multipole to the partialLocalExpansion
    for (Panel* panel : interactionList_)
    {
        partialLocalExpansionAroundCentre_->addLocalExpansionCoefficients(
                panel->multipoleAroundCentre_->convertMultipoleToLocal(centre_)
        );
    }
}

void Panel::computeLocalExpansion()
{
    localExpansionAroundCentre_->addLocalExpansionCoefficients(
            partialLocalExpansionAroundCentre_->getExpansionCoefficients()
    );
}

void Panel::translateLocalExpansion()
{
    for (Panel* iC : childeren_)
    {
        Vec3D centreChild = iC->getCentre();
        iC->localExpansionAroundCentre_->addLocalExpansionCoefficients(
                localExpansionAroundCentre_->translateLocalExpansion(centreChild)
        );
    }
}










