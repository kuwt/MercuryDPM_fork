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
#ifndef CGHandler_H
#define CGHandler_H

#include "GeneralDefine.h"
#include "BaseHandler.h"
#include "CG/BaseCG.h"

class DPMBase;

/*!
 * \todo get function that by default can distinguish species, but also density
 * \todo make tests
 * \todo change output format
 * \todo command line arguments
 * \todo make more readable definitions
 * \todo read (a) restart, (b) data/fstat files
 * \todo add speed by using the mesh
 * \todo introduce standardDev 
 * \todo add 2D support
 * \todo do we need to store the BaseObject::index_ anymore?
 * \todo Can Interaction inherit directly from BaseInteraction?
 * \todo make Files::statFile_, ... public, like the handlers; remove get functions for File's and handlers.
 * \todo take out dependence on DPMBase::statFile (i.e. the savecount) 
 * \todo TW note, to keep the code working on Windows:
 * - std::exit requires correct header cstdlib, 
 * - don't use to_string
 * (thanks to Silvia for debugging)
 */

/*!
 * \class CGHandler
 * \brief Container that stores all CG objects.
 * \details The CGHandler stores all CG objects. Its member functions initialise, 
 * evaluate, and finish are called by the DPMBase before, during and after 
 * the time loop, respectively, and call the member functions initialise, 
 * evaluate, and finish of each CG object.
 * 
 * See also: \ref MercuryCG.
 */
class CGHandler : public BaseHandler<BaseCG>
{
public:
    /*!
     * \brief Default constructor, creates an empty CGHandler.
     */
    CGHandler() = default;;
    
    /*!
     * \brief Copy constructor, copies the CGHandler and all BaseCGPoint's it contains.
     */
    CGHandler(const CGHandler& BH);
    
    /*!
     * \brief Assignment operator that copies the pointer to the DPMBase and all objects.
     */
    CGHandler& operator=(const CGHandler& rhs);
    
    /*!
     * \brief Destructor, destructs the CGHandler and all BaseCGPoint's it contains.
     */
    ~CGHandler() final = default;
    
    void addObject(BaseCG* cg) final;
    
    std::string getName() const final;
    
    /*!
     * \brief Reads objects into the CGHandler from an istream (currently not 
     * implemented).
     */
    void readAndAddObject(std::istream& is) final;
    
    /*!
     * \brief Writes objects into the CGHandler to an ostream (currently not 
     * implemented).
     */
    void write(std::ostream& os) const;
    
    /*!
     * \brief Contains the code executed before the first time step.
     */
    void initialise();
    
    /*!
     * \brief Contains the code executed at each time step.
     */
    void evaluate();
    
    /*!
     * \brief Contains the code executed after the last time step.
     */
    void finish();
    
    /*!
     * \brief loads restart file, before evaluateDataFiles is run
     */
    void restart(std::string name);
    
    void restartAndEvaluateRestartFiles(const std::string& name);
    
    void restartAndEvaluateDataFiles(const std::string& name, bool evaluateFStatFiles = true);
    
    /*!
     * \brief does the same as StatisticsVector::statistics_from_fstat_and_data:
     * loads a restart file (if existing), then several data files, and fstat files (if existing)
     */
    bool evaluateDataFiles(bool evaluateFStatFiles = true);
    
    bool evaluateRestartFiles();
    
    void computeContactPoints();
    
    Mdouble getTimeMin();
    
    Mdouble getTimeMax();

    void setInitialFileCounter(unsigned initialFileCounter) {
        this->initialFileCounter = initialFileCounter;
    }

    unsigned getInitialFileCounter() const {
        return initialFileCounter;
    }

    unsigned initialFileCounter = 0;
};

#endif
