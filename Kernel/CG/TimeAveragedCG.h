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
#ifndef TimeAveragedCG_H
#define TimeAveragedCG_H

#include <GeneralDefine.h>
#include "CG/CG.h"
#include <iostream>
#include <vector>

class BaseParticle;

class BaseInteraction;

class DPMBase;

class Function;

/*!
 * \brief Evaluates time-averaged continuum fields and writes the data into a 
 * stat file.
 * \details
 * Evaluates statistical output that is averaged over all time steps in the time 
 * interval [timeMin_, timeMax_].
 * 
 * This class should only be used to evaluate steady (time-independent) fields. 
 * For time-dependent fields, see CG and TimeSmoothedCG.
 * 
 * The class is derived from CG, which already contains the basic 
 * functionality to obtain time-dependent statistics. This class then averages 
 * the field values over time, with a uniform weight given to all time steps.
 */
template<class Coordinates = CGCoordinates::O,
        template<class> class BaseFunction=CGFunctions::Lucy,
        class Fields=CGFields::StandardFields>
class TimeAveragedCG : public CG<Coordinates, BaseFunction, Fields>
{
public:
    typedef BaseFunction<Coordinates> Function;
    
    /*!
     * \brief Default constructor. 
     */
    TimeAveragedCG();
    
    /*!
     * \brief Copy constructor. It copies the TimeAveragedCGFunction and all objects it contains.
     * \param[in] p the TimeAveragedCGFunction that has to be copied
     */
    TimeAveragedCG(const TimeAveragedCG& p) = default;
    
    /*!
     * \brief Destructor, it simply destructs the TimeAveragedCGFunction and all the objects it contains.
     */
    virtual ~TimeAveragedCG() = default;
    
    /*!
     * \brief
     */
    TimeAveragedCG<Coordinates, BaseFunction, Fields>* copy() const override;
    
    /*!
     * \brief
     */
    void write(std::ostream& os) const override;
    
    /*!
     * \brief
     */
    std::string getName() const override;
    
    /*!
     * \brief Called at the beginning of the DPM simulation to initialise the cg 
     * evaluation and to open the statFile.
     */
    void initialise() override;
    
    /*!
     * \brief Called after a given number of time steps (statFile::saveCount_)
     * to evaluate the CG fields.
     */
    void evaluate() override;
    
    /*!
     * \brief Called at the end of the DPM simulation to finish the cg 
     * evaluation and to close the statFile.
     */
    void finish() override;

protected:
    
    /*!
     * Stores the number of time steps over which the fields are averaged. 
     * Needed, as the average is computed by summing up all values, then 
     * dividing by the number of values, i.e. nTime_.
     * Internal variable that cannot be set by the user.
     */
    unsigned int nTime_;
};


/*!
 * \brief Specialisation of TimeAveragedCG with coordinates XYZ used for LebedevCG
 */
template<template<class> class BaseFunction, class Fields = CGFields::StandardFields>
//class test : public TimeAveragedCG_Lebedev<CGCoordinates::XYZ, BaseFunction, Fields>
class TimeAveragedCGXYZ : public TimeAveragedCG<CGCoordinates::XYZ, BaseFunction, Fields>
{
public:
    /*!
     * \brief Default constructor. 
     */
    TimeAveragedCGXYZ();
    
    /*!
     * \brief Copy constructor. It copies the TimeAveragedCGFunction and all objects it contains.
     * \param[in] p the TimeAveragedCGFunction that has to be copied
     */
    TimeAveragedCGXYZ(const TimeAveragedCGXYZ& p);
    
    /*!
     * \brief Destructor, it simply destructs the TimeAveragedCGFunction and all the objects it contains.
     */
    virtual ~TimeAveragedCGXYZ();
    
    /*!
     * \brief Copy
     */
    TimeAveragedCGXYZ<BaseFunction, Fields>* copy() const override;
    
};

#include "TimeAveragedCG.hcc"

#endif
