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
#ifndef TimeSmoothedFields_H
#define TimeSmoothedFields_H

#include <GeneralDefine.h>
#include <vector>
#include "CG/Fields/StandardFields.h"

/*!
 * \brief A helper class for TimeSmoothedCG containing the time-smoothed 
 * variables.
 * \details
 * Time-smoothing of the Fields is achieved by averaging 
 * the field values F over time,
 * \f[\bar F(t,w_t) = \frac{\sum_i w(t-t_i) F(t_i)}{\sum_i w(t-t_i)},\f]
 * with the weight given by a Gaussian distribution,
 * \f[\bar w(t-t_i) = exp(-(t-t_i)^2/2w_t^2).\f]
 * Thus, the smoothed fields have to be stored over several time steps until the 
 * full sum has been computed. This is done in TimeSmoothedCG::averages_, which 
 * is a vector of TimeSmoothedFields values; for each time step $t$ at which 
 * smoothed values are computed, a TimeSmoothedFields object is added at 
 * \f$t-t_{cutoff}\f$, then values are continually added until the smoothed 
 * field is fully computed at time \f$t+t_{cutoff}\f$; the data is then written 
 * and removed from the vector.
 */

template<class Fields>
class TimeSmoothedFields
{
public:
    
    /*!
     * \brief Constructor; sets the size of the fields_ vector and time_; sets sumWeights_ to zero.
     */
    TimeSmoothedFields(std::size_t n, Mdouble time) : time_(time), sumWeights_(0), fields_(n)
    {}
    
    /*!
     * \brief Default copy Constructor; copies all member variables.
     */
    TimeSmoothedFields(const TimeSmoothedFields& a) = default;
    
    /*!
     * \brief Default destructor; does nothing
     */
    ~TimeSmoothedFields() = default;
    
    /*!
     * A vector of StandardFields values whose length is equal to the number of CGPoints.
     * Is used to compute the time-smoothed fields. 
     */
    std::vector<Fields> fields_;
    
    /*!
     * The time for which smoothed-time fields are evaluated
     */
    Mdouble time_;
    
    /*!
     * The sum of weights of all time steps that are already computed
     */
    Mdouble sumWeights_;
};

#endif
