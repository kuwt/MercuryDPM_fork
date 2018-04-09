//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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

#ifndef CG_H
#define CG_H

#include "CG/BaseCG.h"
#include <vector>

class CGHandler;

class Function;

class BaseParticle;

class BaseInteraction;

#include "CG/CGPoint.h"
#include "DPMBase.h"
#include "CG/CGHandler.h"
#include "Species/ParticleSpecies.h"
#include "Logger.h"
#include "CG/Functions/Lucy.h"
#include "CG/Functions/Heaviside.h"

/*!
 * \class CG
 * \brief Evaluates time-resolved continuum fields and writes the data into a
 * stat file.
 * \details
 * Derived from BaseCG, which already contains the basic functionality of a CG
 * object. This class adds points_, a set of CGPoints,
 * which are templated with a certain CGFunction.
 *
 * \details The CGHandler stores all CG objects. Its member functions initialise,
 * evaluate, and finish are called by the CGHandler before, during and after
 * the time loop, respectively, and contain the main functionality:
 * - open, write to, and close the statFile.
 * - for each time step, calls the functions in CGPoint that cause the
 *   evaluation of the fields.
 *
 * For time-averaged or time-smoothed fields, see TimeAverageCG and
 * TimeSmoothedCG.
 *
 * See also: \ref MercuryCG.
 *
 * \class Function
 * \brief Template argument; use a member class of CGFunctions to instantiate.
 */
template<class Coordinates = CGCoordinates::O,
        template<class> class BaseFunction=CGFunctions::Lucy,
        class Fields=CGFields::StandardFields>
class CG : public BaseCG
{
public:
    /*!
     * \brief Because of this typedefs, Point can be used instead of
     * CGPoint<Function> and Function can be used instead of
     * BaseFunction<Coordinates> in this class.
     * This was done to avoid the overly use of template notation.
     */
    typedef BaseFunction<Coordinates> Function;
    typedef CGPoint<Coordinates, Fields> Point;

    /*!
     * \brief Default constructor; does nothing, i.e. no points are created
     * initially.
     */
    CG();

    /*!
     * \brief Default copy Constructor; copies all member variables.
     */
    CG(const CG& p);

    /*!
     * \brief Default destructor; does nothing
     */
    virtual ~CG();

    /*!
     * \brief Copy operator; creates a new'ed CG object.
     */
    CG<Coordinates, BaseFunction, Fields>* copy() const override;

    /*!
     * \brief Writes class content, except for the points, into an output stream.
     * \todo TW write should be renamed writeHeader, writeAll should be renamed write.
     */
    void write(std::ostream& os) const override;

    /*!
     * \brief Writes class content, including the points_, into an output stream, 
     * usually a stat file.
     */
    void writeAll(std::ostream& os) const;

    /*!
     * \brief returns the name of the class, which is required by write.
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

    Point evaluateTotal();

    /*!
     * \brief Contains the basic for loop over all CGPoints, required to do
     * particle statistics.
     * \todo evaluateParticle and evaluateContact can be optimized by using the 
     * grid properties for smart neighborhood search.
     */
    void evaluateParticle(const BaseParticle& p);

    /*!
     * \brief Contains the basic for loop over all CGPoints, required to do
     * contact statistics.
     * \todo the check for contact statistics should be done here, not per 
     * particle.
     */
    void evaluateContact(const BaseInteraction& i);

    IntegralType getIntegralType(const BaseInteraction& c);

    /*!
     * \brief Called at the end of the DPM simulation to finish the cg
     * evaluation and to close the statFile.
     */
    void finish() override;

    const Point& getPoint(size_t i) const
    {
        return points_[i];
    }

    const std::vector<Point>& getPoints() const
    {
        return points_;
    }


protected:

    /*!
     * \brief set all variables to zero
     */
    void resetVariables();

    /*!
     * \brief divide each variable by the domain volume
     */
    void volumeAverageVariables();

    /*!
     * \brief write variables to the stat file
     */
    void writeVariables();

    /*!
     * \brief The part of evaluate that is used for CG, timeAveragedCG and timeSmoothedCG
     */
    void evaluateCommon();

    /*!
     * \brief plot total to console
     */
    void outputSumOfVariables();

    /*!
     * \brief Contains the CGPoint's, i.e. the positions at which the StandardFields
     * are evaluated.
     */
    std::vector<Point> points_;

    Function function_;


};

#include "CG/CG.hcc"

#endif
