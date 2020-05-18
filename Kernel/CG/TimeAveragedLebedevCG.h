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

#ifndef TimeAveragedLebedevCG_H
#define TimeAveragedLebedevCG_H

#include <CG/TimeAveragedCG.h>

template<template<class> class BaseFunction, class Fields = CGFields::StandardFields>
class TimeAveragedLebedevCG : public TimeAveragedCGXYZ<BaseFunction, Fields>
{
public:
    
    typedef CGPoint<CGCoordinates::XYZ, Fields> Point;
    
    /*!
     * \brief Default constructor. Only sets the evaluation functions,
        no points are created initially.
     */
    TimeAveragedLebedevCG() = default;
    
    /*!
     * \brief Default copy Constructor; copies all the member variables.
     */
    TimeAveragedLebedevCG(const TimeAveragedLebedevCG& p) = default;
    
    /*!
     * \brief Default destructor; does nothing
     */
    virtual ~TimeAveragedLebedevCG();
    
    /*!
     * \brief Creates a mesh based on Lebedev quadrature points
     */
    void createMesh() override;
    
    
    /*!
     * \brief Creates a copy of the current instance 
     */
    TimeAveragedLebedevCG<BaseFunction, Fields>* copy() const;
    
    /*!
     * \brief Returns the inner radius of the grid.
     */
    double getRadiusInner();
    
    /*!
     * \brief Returns the outer radius of the grid.
     */
    double getRadiusOuter();
    
    /*!
     * \brief Returns the number of grid points in r-direction.
     */
    int getNR();
    
    /*!
     * \brief Sets the inner and outer radius of the mesh
     */
    void setR(double radiusInner, double radiusOuter);
    
    /*!
     * \brief Sets the number of mesh points in the r-direction.
     */
    void setNR(int nR);
    
    /*!
     * \brief Evaluates CG fields
     */
    void evaluate();
    
    /*!
     * \brief Evaluates the contributions of a particle.
     */
    void evaluateParticle(BaseParticle& p);
    
    /*!
     * \brief Evaluates the contribution of an interaction.
     */
    void evaluateContact(BaseInteraction& i);

private:
    /*!
     * \TODO document
     */
    double radiusInner_;
    double radiusOuter_;
    int nR_;
};

#include "TimeAveragedLebedevCG.hcc"

#endif
