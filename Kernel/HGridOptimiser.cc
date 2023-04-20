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

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include "HGridOptimiser.h"
#include "Particles/BaseParticle.h"

void HGridOptimiser::initialise(const MercuryBase& problem, unsigned int numberOfCells, int verbosity)
{
    //assign the number of cells
    numCells_ = numberOfCells;
    cellCheckOverContactCheckRatio_ = 0.2;
    
    cellN_.resize(numCells_);
    for (double& it : cellN_)
    {
        it = 0.0;
    }
    
    //Set the minimum and maximum radius of the particles.
    rMin_ = nextafter(problem.particleHandler.getSmallestParticle()->getMaxInteractionRadius(), 0.0);
    rMax_ = nextafter(nextafter(problem.particleHandler.getLargestParticle()->getMaxInteractionRadius(),
                                std::numeric_limits<double>::max()), std::numeric_limits<double>::max());
    
    //for all particles, add one to the cell in cellN_ which has is associated with
    //the radius of that particle.
    for (std::vector<BaseParticle*>::const_iterator it = problem.particleHandler.begin();
         it != problem.particleHandler.end(); ++it)
    {
        cellN_[radius2Cell((*it)->getMaxInteractionRadius())]++;
    }
    
    //assign what the number of particles is.
    unsigned int numParticles = problem.particleHandler.getSize();
    
    intCellN.push_back(1.5 * cellN_[0] - 0.5 * cellN_[1]);
    for (double& it : cellN_)
    {
        intCellN.push_back(it);
    }
    intCellN.push_back(1.5 * cellN_[numCells_ - 1] - 0.5 * cellN_[numCells_ - 2]);
    
    if (verbosity > 0)
    {
        logger(INFO, "rMin=% rMax=% NParticles=%\n", Flusher::NO_FLUSH, rMin_, rMax_, numParticles);
        for (unsigned int i = 0; i < cellN_.size(); i++)
        {
            logger(INFO, "From %.8 to %.8 number=%.5\n", Flusher::NO_FLUSH, cell2Min(i), cell2Max(i), cellN_[i]);
        }
        logger(INFO, "Domain size [%,%]x[%,%]x[%,%]\n",
               Flusher::NO_FLUSH, problem.getXMin(), problem.getXMax(), problem.getYMin(), problem.getYMax(), problem
                       .getZMin(), problem.getZMax());
        for (unsigned int i = 0; i < intCellN.size() - 1; i++)
        {
            logger(INFO, "From %.8 to %.8 linear from %.8 to %.8",
                   intCell2Min(i), intCell2Max(i), intCellN[i], intCellN[i + 1]);
        }
        /*std::cout<<"["<<cellN_[0];
         for (unsigned int i=1;i<cellN_.size();i++)
         {
         std::cout<<","<<cellN_[i];
         }       
         std::cout<<"]"<<std::endl;*/
    }
    
    dimension_ = problem.getSystemDimensions();
    switch (dimension_)
    {
        case 1:
            length_ = pow((problem.getXMax() - problem.getXMin()) / numParticles, 1);
            break;
        case 2:
            length_ = pow(
                    (problem.getXMax() - problem.getXMin()) * (problem.getYMax() - problem.getYMin()) / numParticles,
                    1.0 / 2.0);
            break;
        case 3:
            length_ = pow((problem.getXMax() - problem.getXMin()) * (problem.getYMax() - problem.getYMin()) *
                          (problem.getZMax() - problem.getZMin()) / numParticles, 1.0 / 3.0);
            break;
    }
}

void
HGridOptimiser::initialisePolyFunc(double omega, std::vector<double>& coeff, unsigned int numberOfCells, int verbosity)
{
    numCells_ = numberOfCells;
    cellCheckOverContactCheckRatio_ = 0.2;
    cellN_.resize(numCells_);
    for (double& it : cellN_)
    {
        it = 0;
    }
    
    rMin_ = nextafter(1, 0.0);
    rMax_ = nextafter(omega, std::numeric_limits<double>::max());
    for (unsigned int i = 0; i < cellN_.size(); i++)
    {
        double start = cell2Min(i);
        double end = cell2Max(i);
        for (unsigned int j = 0; j < coeff.size(); j++)
        {
            cellN_[i] += coeff[j] / (1.0 + j) * (std::pow(end, j + 1) - std::pow(start, j + 1));
        }
    }
    
    intCellN.push_back(1.5 * cellN_[0] - 0.5 * cellN_[1]);
    for (double& it : cellN_)
    {
        intCellN.push_back(it);
    }
    intCellN.push_back(1.5 * cellN_[numCells_ - 1] - 0.5 * cellN_[numCells_ - 2]);
    
    dimension_ = 1;
    length_ = 1.0;
    
    if (verbosity > 0)
    {
        logger(INFO, "rMin=% rMax=%\n", Flusher::NO_FLUSH, rMin_, rMax_);
        for (unsigned int i = 0; i < cellN_.size(); i++)
        {
            logger(INFO, "From %.8 to %.8 number=%.5\n", Flusher::NO_FLUSH, cell2Min(i), cell2Max(i), cellN_[i]);
        }
        
        for (unsigned int i = 0; i < intCellN.size() - 1; i++)
        {
            logger(INFO, "From %.8 to %.8 linear from %.8 to %.8",
                   intCell2Min(i), intCell2Max(i), intCellN[i], intCellN[i + 1]);
        }
    }
}

/*!
 * \details Assigns a cell to a BaseParticle with the given radius. Note that the
 *          index of the cells are linear in the radius of a particle. For example,
 *          if numCells_ = 10 and we are looking in the radius range [0,10], than
 *          the cell number of the particle with radius r is the number r rounded
 *          down to an integer.
 */
unsigned int HGridOptimiser::radius2Cell(double r)
{
    if (r < rMin_ || r >= rMax_)
        logger(ERROR, "The radius % is not in the range [%,%]", r, rMin_, rMax_);
    unsigned int y = static_cast<unsigned int> (floor(numCells_ * (r - rMin_) / (rMax_ - rMin_)));
    return y;
}

unsigned int HGridOptimiser::radius2IntCell(double r)
{
    if (r < rMin_ || r >= rMax_)
        throw;
    unsigned int y = static_cast<unsigned int> (floor(numCells_ * (r - rMin_) / (rMax_ - rMin_) + 0.5));
    return y;
}

double HGridOptimiser::intCell2Min(unsigned int i)
{
    if (i == 0)
        return rMin_;
    else
        return rMin_ + (rMax_ - rMin_) * (i - 0.5) / numCells_;
}

double HGridOptimiser::intCell2Max(unsigned int i)
{
    if (i == numCells_)
        return rMax_;
    else
        return rMin_ + (rMax_ - rMin_) * (i + 0.5) / numCells_;
}

/*!
 * Computes the left bound of the cell with given ordinal number.
 */
double HGridOptimiser::cell2Min(unsigned int i)
{
    return rMin_ + (rMax_ - rMin_) * i / numCells_;
}

/*!
 * Computes the right bound of the cell with given ordinal number.
 */
double HGridOptimiser::cell2Max(unsigned int i)
{
    return rMin_ + (rMax_ - rMin_) * (i + 1) / numCells_;
}

double HGridOptimiser::pdfIntCell(double start, double end, unsigned int i, int p)
{
    double a = (intCellN[i + 1] - intCellN[i]) / (intCell2Max(i) - intCell2Min(i));
    double b = (intCellN[i] * intCell2Max(i) - intCellN[i + 1] * intCell2Min(i)) / (intCell2Max(i) - intCell2Min(i));
    return (a / (p + 2) * (std::pow(end, p + 2) - std::pow(start, p + 2)) +
            b / (p + 1) * (std::pow(end, p + 1) - std::pow(start, p + 1)));
}

///This function calculates:
///int(f(r)*r^power*dr,r=start..end)/int(f(r)*dr,r=0..omega)
///with r=a*r+b
double HGridOptimiser::pdfInt(double start, double end, int p)
{
    //std::cout<<"pdfInit("<<start<<","<<end<<","<<p<<");"<<std::endl;
    unsigned int startCell = radius2IntCell(start);
    unsigned int endCell = radius2IntCell(end);
    
    double numerator = 0;
    if (startCell == endCell)
    {
        numerator = pdfIntCell(start, end, startCell, p);
    }
    else if (endCell < startCell)
    {
        numerator = 0;
    }
    else
    {
        numerator = pdfIntCell(start, intCell2Max(startCell), startCell, p);
        //std::cout<<"Teller+="<<pdfIntCell(start,intCell2Max(startCell),startCell,p)<<std::endl;
        for (unsigned int i = startCell + 1; i < endCell; i++)
        {
            numerator += pdfIntCell(intCell2Min(i), intCell2Max(i), i, p);
            //std::cout<<"Teller+="<<pdfIntCell(intCell2Min(i),intCell2Max(i),i,p)<<std::endl;
        }
        numerator += pdfIntCell(intCell2Min(endCell), end, endCell, p);
        //std::cout<<"Teller+="<<pdfIntCell(intCell2Min(endCell),end,endCell,p)<<std::endl;
    }
    double denominator = 0;
    for (unsigned int i = 0; i <= numCells_; i++)
    {
        denominator += pdfIntCell(intCell2Min(i), intCell2Max(i), i, 0);
        //std::cout<<"Noemer+="<<pdfIntCell(intCell2Min(i),intCell2Max(i),i,0)<<std::endl;
    }
    return numerator / denominator;
}

///diff(int(f(r)*r^power*dr,r=s..e)/int(f(r)*dr,r=0..omega),e)=f(e)*e^power/int(f(r)*dr,r=0..omega)
double HGridOptimiser::diffPdfInt(double x, int p)
{
    unsigned int xCell = radius2IntCell(x);
    double denominator = 0;
    for (unsigned int i = 0; i <= numCells_; i++)
    {
        denominator += pdfIntCell(intCell2Min(i), intCell2Max(i), i, 0);
    }
    double a = (intCellN[xCell + 1] - intCellN[xCell]) / (intCell2Max(xCell) - intCell2Min(xCell));
    double b = (intCellN[xCell] * intCell2Max(xCell) - intCellN[xCell + 1] * intCell2Min(xCell)) /
               (intCell2Max(xCell) - intCell2Min(xCell));
    return (a * std::pow(x, p + 1) + b * std::pow(x, p)) / denominator;
}

double HGridOptimiser::expectedCellsIntegralCellNumerator(double start, double end, unsigned int i, int p, double h)
{
    double a = (intCellN[i + 1] - intCellN[i]) / (intCell2Max(i) - intCell2Min(i));
    double b = (intCellN[i] * intCell2Max(i) - intCellN[i + 1] * intCell2Min(i)) / (intCell2Max(i) - intCell2Min(i));
    double r;
    r = start;
    double min = std::pow(2.0 * (r + h) / h, p) * (a * p * r - a * h + a * r + b * p + 2 * b) * (r + h) /
                 ((p + 1) * (p + 2));
    r = end;
    double plus = std::pow(2.0 * (r + h) / h, p) * (a * p * r - a * h + a * r + b * p + 2 * b) * (r + h) /
                  ((p + 1) * (p + 2));
    return plus - min;
}

double
HGridOptimiser::diffHExpectedCellsIntegralCellNumerator(double start, double end, unsigned int i, int p, double h)
{
    double a = (intCellN[i + 1] - intCellN[i]) / (intCell2Max(i) - intCell2Min(i));
    double b = (intCellN[i] * intCell2Max(i) - intCellN[i + 1] * intCell2Min(i)) / (intCell2Max(i) - intCell2Min(i));
    double r;
    r = start;
    //double min =std::pow(2.0*r/h+2.0,p)*(a*p*r-a*h+a*r+b*p+2*b)*(r+h)/((p+1)*(p+2));
    double min =
            (-2.0 * r / h / h * p * std::pow(2.0 * r / h + 2.0, p - 1) * (a * p * r - a * h + a * r + b * p + 2 * b) *
             (r + h) +
             -std::pow(2.0 * r / h + 2.0, p) * a * (r + h) +
             std::pow(2.0 * r / h + 2.0, p) * (a * p * r - a * h + a * r + b * p + 2 * b)) / ((p + 1) * (p + 2));
    r = end;
    //double plus=std::pow(2.0*r/h+2.0,p)*(a*p*r-a*h+a*r+b*p+2*b)*(r+h)/((p+1)*(p+2));
    double plus =
            (-2.0 * r / h / h * p * std::pow(2.0 * r / h + 2.0, p - 1) * (a * p * r - a * h + a * r + b * p + 2 * b) *
             (r + h) +
             -std::pow(2.0 * r / h + 2.0, p) * a * (r + h) +
             std::pow(2.0 * r / h + 2.0, p) * (a * p * r - a * h + a * r + b * p + 2 * b)) / ((p + 1) * (p + 2));
    return plus - min;
}

double HGridOptimiser::expectedCellsIntegralCellDenominator(double start, double end, unsigned int i)
{
    double a = (intCellN[i + 1] - intCellN[i]) / (intCell2Max(i) - intCell2Min(i));
    double b = (intCellN[i] * intCell2Max(i) - intCellN[i + 1] * intCell2Min(i)) / (intCell2Max(i) - intCell2Min(i));
    return (a / 2.0 * (pow(end, 2) - pow(start, 2)) + b * (pow(end, 1) - pow(start, 1)));
}

///This function calculates:
///int((2*r/h+2)^d f(r) dr,r=s..e)/int(f(r) dr,r=s..e)+
///Used to calculated the expected number of cells to check at the level with 
///maximum size h for particle radius between start and end
double HGridOptimiser::expectedCellsIntegral(double start, double end, int p, double h)
{
    unsigned int startCell = radius2IntCell(start);
    unsigned int endCell = radius2IntCell(end);
    
    double numerator = 0;
    double denominator = 0;
    if (startCell == endCell)
    {
        numerator = expectedCellsIntegralCellNumerator(start, end, startCell, p, h);
        denominator = expectedCellsIntegralCellDenominator(start, end, startCell);
    }
    else
    {
        numerator = expectedCellsIntegralCellNumerator(start, intCell2Max(startCell), startCell, p, h);
        denominator = expectedCellsIntegralCellDenominator(start, intCell2Max(startCell), startCell);
        for (unsigned int i = startCell + 1; i < endCell; i++)
        {
            numerator += expectedCellsIntegralCellNumerator(intCell2Min(i), intCell2Max(i), i, p, h);
            denominator += expectedCellsIntegralCellDenominator(intCell2Min(i), intCell2Max(i), i);
        }
        numerator += expectedCellsIntegralCellNumerator(intCell2Min(endCell), end, endCell, p, h);
        denominator += expectedCellsIntegralCellDenominator(intCell2Min(endCell), end, endCell);
        
    }
    //    std::cout<<"New: numerator="<<numerator<<" denominator="<<denominator<<std::endl;
    return numerator / denominator;
}

double
HGridOptimiser::diffStartExpectedCellsIntegral(double start, double end, int p, double h)
{
    unsigned int startCell = radius2IntCell(start);
    unsigned int endCell = radius2IntCell(end);
    
    double numerator = 0;
    double denominator = 0;
    
    double a = (intCellN[startCell + 1] - intCellN[startCell]) / (intCell2Max(startCell) - intCell2Min(startCell));
    double b = (intCellN[startCell] * intCell2Max(startCell) - intCellN[startCell + 1] * intCell2Min(startCell)) /
               (intCell2Max(startCell) - intCell2Min(startCell));
    
    double diffTeller = -(
            2.0 / h * p * std::pow(2.0 * (start + h) / h, p - 1) * (a * p * start - a * h + a * start + b * p + 2 * b) *
            (start + h) / ((p + 1) * (p + 2)) +
            std::pow(2.0 * (start + h) / h, p) * (a * p + a) * (start + h) / ((p + 1) * (p + 2)) +
            std::pow(2.0 * (start + h) / h, p) * (a * p * start - a * h + a * start + b * p + 2 * b) /
            ((p + 1) * (p + 2)));
    double diffNoemer = -(a * start + b);
    if (startCell == endCell)
    {
        numerator = expectedCellsIntegralCellNumerator(start, end, startCell, p, h);
        denominator = expectedCellsIntegralCellDenominator(start, end, startCell);
    }
    else
    {
        numerator = expectedCellsIntegralCellNumerator(start, intCell2Max(startCell), startCell, p, h);
        denominator = expectedCellsIntegralCellDenominator(start, intCell2Max(startCell), startCell);
        for (unsigned int i = startCell + 1; i < endCell; i++)
        {
            numerator += expectedCellsIntegralCellNumerator(intCell2Min(i), intCell2Max(i), i, p, h);
            denominator += expectedCellsIntegralCellDenominator(intCell2Min(i), intCell2Max(i), i);
        }
        numerator += expectedCellsIntegralCellNumerator(intCell2Min(endCell), end, endCell, p, h);
        denominator += expectedCellsIntegralCellDenominator(intCell2Min(endCell), end, endCell);
        
    }
    //std::cout<<"new: numerator="<<numerator<<" denominator="<<denominator<<" diffTeller="<<diffTeller<<" diffNoemer="<<diffNoemer<<std::endl;
    return (diffTeller * denominator - numerator * diffNoemer) / (denominator * denominator);
}

double HGridOptimiser::diffEndExpectedCellsIntegral(double start, double end, int p, double h)
{
    unsigned int startCell = radius2IntCell(start);
    unsigned int endCell = radius2IntCell(end);
    
    double numerator = 0;
    double denominator = 0;
    
    double a = (intCellN[endCell + 1] - intCellN[endCell]) / (intCell2Max(endCell) - intCell2Min(endCell));
    double b = (intCellN[endCell] * intCell2Max(endCell) - intCellN[endCell + 1] * intCell2Min(endCell)) /
               (intCell2Max(endCell) - intCell2Min(endCell));
    
    double diffTeller = (
            2.0 / h * p * std::pow(2.0 * (end + h) / h, p - 1) * (a * p * end - a * h + a * end + b * p + 2 * b) *
            (end + h) / ((p + 1) * (p + 2)) +
            std::pow(2.0 * (end + h) / h, p) * (a * p + a) * (end + h) / ((p + 1) * (p + 2)) +
            std::pow(2.0 * (end + h) / h, p) * (a * p * end - a * h + a * end + b * p + 2 * b) / ((p + 1) * (p + 2)));
    double diffNoemer = (a * end + b);
    if (startCell == endCell)
    {
        numerator = expectedCellsIntegralCellNumerator(start, end, startCell, p, h);
        denominator = expectedCellsIntegralCellDenominator(start, end, startCell);
    }
    else
    {
        numerator = expectedCellsIntegralCellNumerator(start, intCell2Max(startCell), startCell, p, h);
        denominator = expectedCellsIntegralCellDenominator(start, intCell2Max(startCell), startCell);
        for (unsigned int i = startCell + 1; i < endCell; i++)
        {
            numerator += expectedCellsIntegralCellNumerator(intCell2Min(i), intCell2Max(i), i, p, h);
            denominator += expectedCellsIntegralCellDenominator(intCell2Min(i), intCell2Max(i), i);
        }
        numerator += expectedCellsIntegralCellNumerator(intCell2Min(endCell), end, endCell, p, h);
        denominator += expectedCellsIntegralCellDenominator(intCell2Min(endCell), end, endCell);
        
    }
    //std::cout<<"new: numerator="<<numerator<<" denominator="<<denominator<<" diffTeller="<<diffTeller<<" diffNoemer="<<diffNoemer<<std::endl;
    return (diffTeller * denominator - numerator * diffNoemer) / (denominator * denominator);
}

double HGridOptimiser::diffHExpectedCellsIntegral(double start, double end, int p, double h)
{
    unsigned int startCell = radius2IntCell(start);
    unsigned int endCell = radius2IntCell(end);
    
    double denominator = 0;
    double diffTeller = 0;
    
    if (startCell == endCell)
    {
        diffTeller = diffHExpectedCellsIntegralCellNumerator(start, end, startCell, p, h);
        denominator = expectedCellsIntegralCellDenominator(start, end, startCell);
    }
    else
    {
        diffTeller = diffHExpectedCellsIntegralCellNumerator(start, intCell2Max(startCell), startCell, p, h);
        denominator = expectedCellsIntegralCellDenominator(start, intCell2Max(startCell), startCell);
        for (unsigned int i = startCell + 1; i < endCell; i++)
        {
            diffTeller += diffHExpectedCellsIntegralCellNumerator(intCell2Min(i), intCell2Max(i), i, p, h);
            denominator += expectedCellsIntegralCellDenominator(intCell2Min(i), intCell2Max(i), i);
        }
        diffTeller += diffHExpectedCellsIntegralCellNumerator(intCell2Min(endCell), end, endCell, p, h);
        denominator += expectedCellsIntegralCellDenominator(intCell2Min(endCell), end, endCell);
    }
    return diffTeller / denominator;
}

void
HGridOptimiser::calculateDiffWork(std::vector<double>& hGridCellSizes, std::vector<double>& dfdx, HGridMethod method,
                                  int verbosity)
{
    unsigned int NLevels = hGridCellSizes.size() - 1;
    unsigned int NParams = hGridCellSizes.size();
    
    if (verbosity > 0)
    {
        logger(INFO, "Length scale=%\n", Flusher::NO_FLUSH, length_);
        logger(INFO, "Level sizes:", Flusher::NO_FLUSH);
        for (unsigned int i = 0; i < NParams; i++)
        {
            logger(INFO, " %", hGridCellSizes[i], Flusher::NO_FLUSH);
        }
        logger(INFO, "");
    }
    
    std::vector<double> particlesPerLevel;
    for (unsigned int i = 0; i < NLevels; i++)
    {
        particlesPerLevel.push_back(pdfInt(0.5 * hGridCellSizes[i], 0.5 * hGridCellSizes[i + 1], 0));
        if (particlesPerLevel.back() < 0)
        {
            particlesPerLevel.back() = 0;
        }
    }
    if (verbosity > 0)
    {
        logger(INFO, "Particles per level:\n", Flusher::NO_FLUSH);
        for (unsigned int i = 0; i < NLevels; i++)
        {
            logger(INFO, " %", particlesPerLevel[i], Flusher::NO_FLUSH);
        }
        logger(INFO, "");
    }
    
    //How changes particlesPerLeve[i] when hGridCellSize[j] is changed
    std::vector<std::vector<double> > diffParticlesPerLevel;
    diffParticlesPerLevel.resize(NLevels);
    for (unsigned int i = 0; i < NLevels; i++)
    {
        for (unsigned int j = 0; j < NParams; j++)
        {
            diffParticlesPerLevel[i].push_back(0.0);
            if (j == i)
            {
                diffParticlesPerLevel[i][j] += -0.5 * diffPdfInt(0.5 * hGridCellSizes[i], 0);
            }
            if (j == i + 1)
            {
                diffParticlesPerLevel[i][j] += 0.5 * diffPdfInt(0.5 * hGridCellSizes[i + 1], 0);
            }
        }
    }
    if (verbosity > 0)
    {
        logger(INFO, "Diff Particles per level: \n", Flusher::NO_FLUSH);
        for (unsigned int i = 0; i < NLevels; i++)
        {
            for (unsigned int j = 0; j < NParams; j++)
            {
                logger(INFO, " %.12", diffParticlesPerLevel[i][j], Flusher::NO_FLUSH);
            }
            logger(INFO, "");
        }
    }
    
    std::vector<double> cellsPerLevel;
    for (unsigned int i = 0; i < NLevels; i++)
    {
        cellsPerLevel.push_back(pow(length_ / hGridCellSizes[i + 1], dimension_));
    }
    if (verbosity > 0)
    {
        logger(INFO, "Cells per level:", Flusher::NO_FLUSH);
        for (unsigned int i = 0; i < NLevels; i++)
        {
            logger(INFO, " %", cellsPerLevel[i], Flusher::NO_FLUSH);
        }
        logger(INFO, "");
    }
    
    //How changes cellsPerLevel[i] when hGridCellSize[j] is changed
    std::vector<std::vector<double> > diffCellsPerLevel;
    diffCellsPerLevel.resize(hGridCellSizes.size());
    for (unsigned int i = 0; i < NLevels; i++)
    {
        for (unsigned int j = 0; j < NParams; j++)
        {
            if (j == i + 1)
            {
                diffCellsPerLevel[i].push_back(
                        -pow(length_ / hGridCellSizes[i + 1], dimension_) * dimension_ / hGridCellSizes[i + 1]);
            }
            else
            {
                diffCellsPerLevel[i].push_back(0.0);
            }
        }
    }
    if (verbosity > 0)
    {
        logger(INFO, "Diff Cells per level: \n", Flusher::NO_FLUSH);
        for (unsigned int i = 0; i < NLevels; i++)
        {
            for (unsigned int j = 0; j < NParams; j++)
            {
                logger(INFO, " %.12", diffCellsPerLevel[i][j], Flusher::NO_FLUSH);
            }
            logger(INFO, "");
        }
    }
    
    std::vector<double> particlesPerCell;
    for (unsigned int i = 0; i < NLevels; i++)
    {
        particlesPerCell.push_back(particlesPerLevel[i] / cellsPerLevel[i]);
    }
    if (verbosity > 0)
    {
        logger(INFO, "Particles per cell:", Flusher::NO_FLUSH);
        for (double i : particlesPerCell)
        {
            logger(INFO, " %", i, Flusher::NO_FLUSH);
        }
        logger(INFO, "");
    }
    
    //How changes particlesPerCell[i] when hGridCellSize[j] is changed
    std::vector<std::vector<double> > diffParticlesPerCell;
    diffParticlesPerCell.resize(NLevels);
    for (unsigned int i = 0; i < NLevels; i++)
    {
        for (unsigned int j = 0; j < NParams; j++)
        {
            diffParticlesPerCell[i].push_back(
                    (diffParticlesPerLevel[i][j] * cellsPerLevel[i] - particlesPerLevel[i] * diffCellsPerLevel[i][j]) /
                    (cellsPerLevel[i] * cellsPerLevel[i]));
        }
    }
    if (verbosity > 0)
    {
        logger(INFO, "Diff Particles per Cell: \n", Flusher::NO_FLUSH);
        for (unsigned int i = 0; i < NLevels; i++)
        {
            for (unsigned int j = 0; j < NParams; j++)
            {
                logger(INFO, " %.12", diffParticlesPerCell[i][j], Flusher::NO_FLUSH);
            }
            logger(INFO, "");
        }
    }
    
    double contactWork = 0, overheadWork = 0;
    std::vector<double> diffContactWork;
    std::vector<double> diffOverheadWork;
    diffContactWork.resize(NParams);
    diffOverheadWork.resize(NParams);
    for (unsigned int j = 0; j < NParams; j++)
    {
        diffContactWork[j] = 0;
        diffOverheadWork[j] = 0;
    }
    
    for (unsigned int i = 0; i < NLevels; i++)
    {
        contactWork += 0.5 * pow(3, dimension_) * particlesPerCell[i] * particlesPerLevel[i];
        overheadWork += 0.5 * (pow(3, dimension_) + 1) * particlesPerLevel[i];
        for (unsigned int j = 0; j < NParams; j++)
        {
            diffContactWork[j] += 0.5 * pow(3, dimension_) * (diffParticlesPerCell[i][j] * particlesPerLevel[i] +
                                                              particlesPerCell[i] * diffParticlesPerLevel[i][j]);
            diffOverheadWork[j] += 0.5 * (pow(3, dimension_) + 1) * diffParticlesPerLevel[i][j];
        }
        
        unsigned int startJ, endJ;
        switch (method)
        {
            case BOTTOMUP:
            {
                startJ = i + 1;
                endJ = NLevels;
                break;
            }
            case TOPDOWN:
            {
                startJ = 0;
                endJ = i;
                break;
            }
        }
        for (unsigned int j = startJ; j < endJ; j++)
        {
            double numberOfCellToCheck;
            numberOfCellToCheck = expectedCellsIntegral(0.5 * hGridCellSizes[i], 0.5 * hGridCellSizes[i + 1],
                                                        dimension_, hGridCellSizes[j + 1]);
            
            std::vector<double> diffNumberOfCellToCheck;
            for (unsigned int k = 0; k < NParams; k++)
            {
                diffNumberOfCellToCheck.push_back(0.0);
                if (k == i)
                {
                    diffNumberOfCellToCheck[k] += 0.5 * diffStartExpectedCellsIntegral(0.5 * hGridCellSizes[i],
                                                                                       0.5 * hGridCellSizes[i + 1],
                                                                                       dimension_,
                                                                                       hGridCellSizes[j + 1]);
                }
                if (k == i + 1)
                {
                    diffNumberOfCellToCheck[k] += 0.5 * diffEndExpectedCellsIntegral(0.5 * hGridCellSizes[i],
                                                                                     0.5 * hGridCellSizes[i + 1],
                                                                                     dimension_, hGridCellSizes[j + 1]);
                }
                if (k == j + 1)
                {
                    diffNumberOfCellToCheck[k] += diffHExpectedCellsIntegral(0.5 * hGridCellSizes[i],
                                                                             0.5 * hGridCellSizes[i + 1], dimension_,
                                                                             hGridCellSizes[j + 1]);
                }
            }
            contactWork += particlesPerLevel[i] * numberOfCellToCheck * particlesPerCell[j];
            overheadWork += particlesPerLevel[i] * numberOfCellToCheck;
            //std::cout<<"i= "<<i<<" j= "<<j<<" NumberOfCellToCheck= "<<numberOfCellToCheck<<" diffNumberOfCellToCheck[2]= "<<diffNumberOfCellToCheck[2]<<std::endl;
            for (unsigned int k = 0; k < NParams; k++)
            {
                diffContactWork[k] += (diffParticlesPerLevel[i][k] * numberOfCellToCheck * particlesPerCell[j] +
                                       particlesPerLevel[i] * diffNumberOfCellToCheck[k] * particlesPerCell[j] +
                                       particlesPerLevel[i] * numberOfCellToCheck * diffParticlesPerCell[j][k]);
                diffOverheadWork[k] += (diffParticlesPerLevel[i][k] * numberOfCellToCheck +
                                        particlesPerLevel[i] * diffNumberOfCellToCheck[k]);
            }
        }
    }
    if (verbosity > 0)
    {
        logger(INFO, "Contact work: %\n"
                     "Overhead work: %\n"
                     "Sum work: %\n",
               contactWork, cellCheckOverContactCheckRatio_ * overheadWork,
               contactWork + cellCheckOverContactCheckRatio_ * overheadWork, Flusher::NO_FLUSH);
        for (unsigned int j = 0; j < NParams; j++)
        {
            logger(INFO, "diff Contactwork: %", diffContactWork[j], Flusher::NO_FLUSH);
        }
        logger(INFO, "\ndiff Overheadwork: ", Flusher::NO_FLUSH);
        for (unsigned int j = 0; j < NParams; j++)
        {
            logger(INFO, " %", cellCheckOverContactCheckRatio_ * diffOverheadWork[j], Flusher::NO_FLUSH);
        }
        logger(INFO, "\ndiff sum: ", Flusher::NO_FLUSH);
        for (unsigned int j = 0; j < NParams; j++)
        {
            logger(INFO, " %", diffContactWork[j] + cellCheckOverContactCheckRatio_ * diffOverheadWork[j],
                   Flusher::NO_FLUSH);
        }
        logger(INFO, "");
    }
    dfdx.clear();
    for (unsigned int j = 0; j < NParams; j++)
    {
        dfdx.push_back(diffContactWork[j] + cellCheckOverContactCheckRatio_ * diffOverheadWork[j]);
    }
}

double HGridOptimiser::calculateWork(std::vector<double>& hGridCellSizes, HGridMethod method, int verbosity)
{
    unsigned int NLevels = hGridCellSizes.size() - 1;
    unsigned int NParams = hGridCellSizes.size();
    
    if (verbosity > 1)
    {
        logger(INFO, "Length scale=%", length_);
    }
    if (verbosity > 0)
    {
        logger(INFO, "Level sizes:\n", Flusher::NO_FLUSH);
        for (unsigned int i = 0; i < NParams; i++)
        {
            logger(INFO, "%.10 ", hGridCellSizes[i], Flusher::NO_FLUSH);
        }
        logger(INFO, "");
    }
    
    std::vector<double> particlesPerLevel;
    for (unsigned int i = 0; i < NLevels; i++)
    {
        particlesPerLevel.push_back(pdfInt(0.5 * hGridCellSizes[i], 0.5 * hGridCellSizes[i + 1], 0));
        if (particlesPerLevel.back() < 0)
        {
            particlesPerLevel.back() = 0;
        }
    }
    if (verbosity > 0)
    {
        logger(INFO, "Particles per level:\n", Flusher::NO_FLUSH);
        for (unsigned int i = 0; i < NLevels; i++)
        {
            logger(INFO, "%.10 ", particlesPerLevel[i], Flusher::NO_FLUSH);
        }
        logger(INFO, "");
    }
    
    std::vector<double> cellsPerLevel;
    for (unsigned int i = 0; i < NLevels; i++)
    {
        cellsPerLevel.push_back(pow(length_ / hGridCellSizes[i + 1], dimension_));
    }
    if (verbosity > 1)
    {
        logger(INFO, "Cells per level:\n", Flusher::NO_FLUSH);
        for (unsigned int i = 0; i < NLevels; i++)
        {
            logger(INFO, "%.10 ", cellsPerLevel[i], Flusher::NO_FLUSH);
        }
        logger(INFO, "");
    }
    
    std::vector<double> particlesPerCell;
    for (unsigned int i = 0; i < NLevels; i++)
    {
        particlesPerCell.push_back(particlesPerLevel[i] / cellsPerLevel[i]);
    }
    if (verbosity > 1)
    {
        logger(INFO, "Particles per cell:\n", Flusher::NO_FLUSH);
        for (double i : particlesPerCell)
        {
            logger(INFO, "%.10 ", i, Flusher::NO_FLUSH);
        }
        logger(INFO, "");
    }
    
    std::vector<std::vector<double> > contactWork;
    std::vector<std::vector<double> > overheadWork;
    contactWork.resize(NLevels);
    overheadWork.resize(NLevels);
    for (unsigned int i = 0; i < NLevels; i++)
    {
        contactWork[i].resize(NLevels);
        overheadWork[i].resize(NLevels);
        for (unsigned int j = 0; j < NLevels; j++)
        {
            contactWork[i][j] = 0;
            overheadWork[i][j] = 0;
        }
    }
    
    for (unsigned int i = 0; i < NLevels; i++)
    {
        contactWork[i][i] += 0.5 * pow(3, dimension_) * particlesPerCell[i] * particlesPerLevel[i];
        overheadWork[i][i] += 0.5 * (pow(3, dimension_) + 1) * particlesPerLevel[i];
        
        unsigned int startJ, endJ;
        switch (method)
        {
            case BOTTOMUP:
            {
                startJ = i + 1;
                endJ = NLevels;
                break;
            }
            case TOPDOWN:
            {
                startJ = 0;
                endJ = i;
                break;
            }
        }
        for (unsigned int j = startJ; j < endJ; j++)
        {
            double numberOfCellToCheck;
            numberOfCellToCheck = expectedCellsIntegral(0.5 * hGridCellSizes[i], 0.5 * hGridCellSizes[i + 1],
                                                        dimension_, hGridCellSizes[j + 1]);
            contactWork[i][j] += particlesPerLevel[i] * numberOfCellToCheck * particlesPerCell[j];
            overheadWork[i][j] += particlesPerLevel[i] * numberOfCellToCheck;
            
        }
    }
    
    double totalContactWork = 0, totalOverheadWork = 0;
    for (unsigned int i = 0; i < NLevels; i++)
    {
        for (unsigned int j = 0; j < NLevels; j++)
        {
            totalContactWork += contactWork[i][j];
            totalOverheadWork += overheadWork[i][j];
        }
    }
    
    if (verbosity > 0)
    {
        logger(INFO, "Contact work: %\n", Flusher::NO_FLUSH, totalContactWork);
        for (unsigned int i = 0; i < NLevels; i++)
        {
            for (unsigned int j = 0; j < NLevels; j++)
            {
                logger(INFO, "%.10 ", contactWork[i][j], Flusher::NO_FLUSH);
            }
            logger(INFO, "");
        }
        logger(INFO, "Overhead work: %\n", Flusher::NO_FLUSH, totalOverheadWork);
        for (unsigned int i = 0; i < NLevels; i++)
        {
            for (unsigned int j = 0; j < NLevels; j++)
            {
                logger(INFO, "%.10 ", overheadWork[i][j], Flusher::NO_FLUSH);
            }
            logger(INFO, "");
        }
        logger(INFO, "Total work: %", totalContactWork + totalOverheadWork);
    }
    return totalContactWork + cellCheckOverContactCheckRatio_ * totalOverheadWork;
}

void HGridOptimiser::calcDfDx(std::vector<double>& hGridCellSizes, std::vector<double>& dfdx, HGridMethod method,
                              int verbosity)
{
    if (verbosity > 0)
    {
        logger(VERBOSE, "calcDfDx\n", Flusher::NO_FLUSH);
    }
    double W = calculateWork(hGridCellSizes, method, verbosity - 1);
    double dx = 1e-7;
    double nW;
    dfdx.clear();
    for (unsigned int i = 0; i < hGridCellSizes.size(); i++)
    {
        hGridCellSizes[i] += dx;
        nW = calculateWork(hGridCellSizes, method, verbosity - 1);
        dfdx.push_back((nW - W) / dx);
        hGridCellSizes[i] -= dx;
        if (verbosity > 0)
        {
            logger(INFO, "i= Value=% Change=% dfdx=%", i, hGridCellSizes[i], nW - W, dfdx.back());
        }
    }
}

double HGridOptimiser::checkLimit(std::vector<double>& hGridCellSizes, std::vector<double>& dfdx, int verbosity)
{
    if (verbosity > 0)
    {
        logger(VERBOSE, "checkLimits part 1 remove levels");
    }
    
    dfdx[0] = 0;
    dfdx.back() = 0;
    for (unsigned int i = 1; i < hGridCellSizes.size(); i++)
    {
        if (std::abs(hGridCellSizes[i - 1] - hGridCellSizes[i]) < 1e-6)
        {
            if (dfdx[i - 1] < dfdx[i])
            {
                if (verbosity > 0)
                {
                    logger(INFO, "Remove level %", Flusher::NO_FLUSH, i);
                }
                if (i > 0.5 * hGridCellSizes.size())
                {
                    hGridCellSizes.erase(hGridCellSizes.begin() + i - 1);
                    dfdx.erase(dfdx.begin() + i - 1);
                }
                else
                {
                    hGridCellSizes.erase(hGridCellSizes.begin() + i);
                    dfdx.erase(dfdx.begin() + i);
                }
                i--;
            }
        }
    }
    
    if (verbosity > 0)
    {
        logger(VERBOSE, "checkLimits part 2 calculate limit\n", Flusher::NO_FLUSH);
        for (unsigned int i = 0; i < hGridCellSizes.size(); i++)
        {
            logger(INFO, "i=% Value=% dfdx=%\n", Flusher::NO_FLUSH, i, hGridCellSizes[i], dfdx[i]);
        }
    }
    double maxStepSize = -std::numeric_limits<double>::max();
    double nmax;
    for (unsigned int i = 1; i < hGridCellSizes.size(); i++)
    {
        nmax = (hGridCellSizes[i - 1] - hGridCellSizes[i]) / (dfdx[i] - dfdx[i - 1]);
        if (verbosity > 0)
        {
            logger(INFO, "Max f=% because otherwise %+x*% > %+x*%\n", Flusher::NO_FLUSH,
                   nmax, hGridCellSizes[i - 1], dfdx[i - 1], hGridCellSizes[i], dfdx[i]);
        }
        if (nmax < 0)
        {
            maxStepSize = std::max(nmax, maxStepSize);
        }
    }
    if (verbosity > 0)
    {
        logger(INFO, "maxStepSize=%", maxStepSize);
    }
    return maxStepSize;
}

void HGridOptimiser::applyStep(std::vector<double>& hGridCellSizes, std::vector<double>& dfdx, double stepsize,
                               int verbosity)
{
    if (verbosity > 0)
    {
        logger(VERBOSE, "Apply step:\n", Flusher::NO_FLUSH);
    }
    for (unsigned int i = 0; i < hGridCellSizes.size() - 1; i++)
    {
        hGridCellSizes[i] += stepsize * dfdx[i];
        if (verbosity > 0)
        {
            logger(INFO, "hGridCellSizes[i]+%*%=%", stepsize, dfdx[i], hGridCellSizes[i]);
        }
    }
}

double
HGridOptimiser::goldenSectionSearch(std::vector<double>& startHGridCellSizes, std::vector<double>& searchDirection,
                                    double min, double cur, double max, HGridMethod method, int verbosity)
{
    std::vector<double> curHGridCellSizes = startHGridCellSizes;
    std::vector<double> xHGridCellSizes = startHGridCellSizes;
    
    applyStep(curHGridCellSizes, searchDirection, cur, verbosity - 1);
    double curWork = calculateWork(curHGridCellSizes, method, verbosity - 1);
    
    double x;
    double resphi = 2.0 - 0.5 * (1 + std::sqrt(5));
    
    //Determine new probing point
    if (max - cur > cur - min)
    {
        //x between cur and max
        x = cur + resphi * (max - cur);
    }
    else
    {
        //x between min and cur
        x = cur - resphi * (cur - min);
    }
    
    if (std::abs(max - min) < 1e-7 * (std::abs(cur) + std::abs(x)))
    {
        return 0.5 * (min + max);
    }
    
    applyStep(xHGridCellSizes, searchDirection, x, 0);
    double xWork = calculateWork(xHGridCellSizes, method, 0);
    if (verbosity > 0)
    {
        logger(INFO, "min=% max=% cur=% curWork=% x=% xWork=%", min, max, cur, curWork, x, xWork);
    }
    if (xWork < curWork) //x is better
    {
        if (max - cur > cur - min)
        {
            return goldenSectionSearch(startHGridCellSizes, searchDirection, cur, x, max, method, verbosity);
        }
        else
        {
            return goldenSectionSearch(startHGridCellSizes, searchDirection, min, x, cur, method, verbosity);
        }
    }
    else //cur is better
    {
        if (max - cur > cur - min)
        {
            return goldenSectionSearch(startHGridCellSizes, searchDirection, min, cur, x, method, verbosity);
        }
        else
        {
            return goldenSectionSearch(startHGridCellSizes, searchDirection, x, cur, max, method, verbosity);
        }
    }
}

void HGridOptimiser::getOptimalDistribution(std::vector<double>& hGridCellSizes, unsigned int numberOfLevels,
                                            HGridMethod method, int verbosity)
{
    hGridCellSizes.clear();
    for (unsigned int i = 0; i <= numberOfLevels; i++)
    {
        hGridCellSizes.push_back(2.0 * (rMin_ + (rMax_ - rMin_) * i / numberOfLevels));
    }
    
    std::vector<double> dfdx;
    double step, max, W, oW;
    W = calculateWork(hGridCellSizes, method, verbosity - 3);
    int stepNumber = 0;
    do
    {
        oW = W;
        calculateDiffWork(hGridCellSizes, dfdx, method, verbosity - 3);
        max = checkLimit(hGridCellSizes, dfdx, verbosity - 3);
        step = goldenSectionSearch(hGridCellSizes, dfdx, max, 0.5 * max, 0, method, verbosity - 2);
        applyStep(hGridCellSizes, dfdx, step, verbosity - 3);
        W = calculateWork(hGridCellSizes, method, verbosity - 3);
        while (W > oW)
        {
            if (verbosity > 1)
            {
                logger(INFO, "% Bad step, trying closer search range\n", Flusher::NO_FLUSH, stepNumber);
            }
            applyStep(hGridCellSizes, dfdx, -step, verbosity - 3);
            max /= 2;
            step = goldenSectionSearch(hGridCellSizes, dfdx, max, 0.5 * max, 0, method, verbosity - 2);
            applyStep(hGridCellSizes, dfdx, step, verbosity - 3);
            W = calculateWork(hGridCellSizes, method, verbosity - 3);
        }
        ++stepNumber;
        if (verbosity > 1)
        {
            logger(INFO, "% Correct step: old work=% new work=% difference=& step=%\n", Flusher::NO_FLUSH,
                   stepNumber, oW, W, oW - W, step / max);
        }
    } while (oW - W > 1e-13 && stepNumber);
    calculateDiffWork(hGridCellSizes, dfdx, method, verbosity - 2);
    if (verbosity > 0)
    {
        logger(INFO, "Work=%\nOptimal cell sizes:", W, Flusher::NO_FLUSH);
        for (double hGridCellSize : hGridCellSizes)
        {
            logger(INFO, " %", hGridCellSize, Flusher::NO_FLUSH);
        }
        logger(INFO, "");
    }
}

void HGridOptimiser::histNumberParticlesPerCell(std::vector<double>& hGridCellSizes)
{
    unsigned int NLevels = hGridCellSizes.size() - 1;
    
    std::vector<double> particlesPerLevel;
    for (unsigned int i = 0; i < NLevels; i++)
    {
        particlesPerLevel.push_back(pdfInt(0.5 * hGridCellSizes[i], 0.5 * hGridCellSizes[i + 1], 0));
        if (particlesPerLevel.back() < 0)
        {
            particlesPerLevel.back() = 0;
        }
    }
    
    std::vector<double> cellsPerLevel;
    for (unsigned int i = 0; i < NLevels; i++)
    {
        cellsPerLevel.push_back(pow(length_ / hGridCellSizes[i + 1], dimension_));
    }
    
    for (unsigned int i = 0; i < NLevels; i++)
    {
        double l = particlesPerLevel[i] / cellsPerLevel[i];
        logger(INFO, "Histogram for level %:", i, Flusher::NO_FLUSH);
        for (unsigned int k = 0; k <= 10; k++)
        {
            logger(INFO, " %.8",
                   std::floor(std::pow(l, k) * std::exp(-l) / mathsFunc::factorial(k) * 1e4 * cellsPerLevel[i]),
                   Flusher::NO_FLUSH);
        }
        logger(INFO, "");
    }
}
