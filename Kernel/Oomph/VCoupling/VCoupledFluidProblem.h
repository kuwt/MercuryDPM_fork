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

#include "Mercury3D.h"
#include "Oomph/FluidProblem.h"
#include "FluidCoupling.h"

using namespace oomph;

template<class ELEMENT>
class VCoupledFluidProblem : public FluidCoupling<Mercury3D, FluidProblem<FluidCoupledElement<ELEMENT>>>
{
public:
    void computeAdditionalForces() override
    {
        unsigned long int nEl = this->fluid_mesh_pt()->nelement();
        for (unsigned int iEl = 0; iEl < nEl; iEl++)
        {
            auto elPtr = dynamic_cast<ELEMENT *>(this->fluid_mesh_pt()->element_pt(iEl));
            
            computeDragForce(elPtr);
            computePressureGradient(elPtr);
            computeShearForce(elPtr);
            computeMagnusLift(elPtr);
            computeSaffmanLift(elPtr);
    
            /*
            const unsigned int iNode = elPtr->nnode();
            oomph::Vector<double> posNode0(3, 0.0);
            oomph::Vector<double> posNodeEnd(3, 0.0);
            // Note that index of end node is iNode - 1, else segmentation fault as node_pt(iNode) does not exist
            elPtr->node_pt(0)->position(posNode0);
            elPtr->node_pt(iNode - 1)->position(posNodeEnd);
            Vec3D min = convertVecFuncs::convertToVec3D(posNode0);
            Vec3D max = convertVecFuncs::convertToVec3D(posNodeEnd);
            oomph::Vector<double> pPos;
            
            
            for (int pId = 0; pId < this->particleHandler.getNumberOfObjects(); pId++)
            {
                oomph::Vector<double> s(3, 0.0);
                GeomObject* geom_obj_pt = nullptr;
                elPtr->locate_zeta(convertVecFuncs::convertToOomphVec(this->particleHandler.getObject(pId)->getPosition()), geom_obj_pt, s);
                if (geom_obj_pt != nullptr)
                {
                    pList.push_back(this->particleHandler.getObject(pId));
                }
            }
    
            for (auto p_ : pList)
            {
                logger(INFO,"Trying to compute coupling forces for iEl % and particle %", iEl, p_->getId());
                // Compute all coupling forces
                Vec3D dummy = {0., 0., 0.};
                dummy += this->getFractionDragForce()* computeDragForce(p_, elPtr);
                
                logger(INFO, "Using Stokes, Fd(iPart = %) = (%, %, %)",p_->getId(), dummy.X,dummy.Y,dummy.Z);
                
                //dummy += this->getFractionPressureGradient()* computePressureGradient(p_);
                //dummy += this->getFractionShearForce()* computeShearForce(p_);
                //dummy += this->getFractionMagnusLift()* computeMagnusLift(p_);
                //dummy += this->getFractionSaffmanLift()* computeSaffmanLift(p_);
    
                //p->addForce(dummy);
            }
            
            //logger(INFO,"\n");
            pList.clear();
            
            // Loop to next element
             */
        }
    }
    
    void computeDragForce(ELEMENT* elPtr_)
    {
        if (this->getFractionDragForce() <= 0.0)
        {
            logger(DEBUG,"Drag force is not activated");
            return;
        }
        
        logger(DEBUG,"In computeDragForce()");
        oomph::Vector<double> s(3,0.0);
        GeomObject* geom_obj_pt = nullptr;
        oomph::Vector<double> velAtPos(3,0.0);
    
        for (int pId = 0; pId < this->particleHandler.getNumberOfObjects(); pId++)
        {
            auto p_ = this->particleHandler.getObject(pId);
            double rad = p_->getRadius();
    
            oomph::Vector<double> s(3, 0.0);
            GeomObject* geom_obj_pt = nullptr;
            elPtr_->locate_zeta(convertVecFuncs::convertToOomphVec(this->particleHandler.getObject(pId)->getPosition()), geom_obj_pt, s);
            
            if (geom_obj_pt != nullptr)
            {
                logger(DEBUG,"Trying to compute coupling forces for iEl and particle %", p_->getId());
    
                elPtr_->interpolated_u_nst(s,velAtPos);
                Vec3D relVel = p_->getVelocity() - convertVecFuncs::convertToVec3D(velAtPos);
                double relVelMagnitude = sqrt(relVel.X*relVel.X + relVel.Y*relVel.Y + relVel.Z*relVel.Z);
    
                ///FIXME For now voidInEl is set to 1;
                // double Re_p = voidInEl * this->getFluidDensity() * relVelMagnitude * 2.0 * p_->getRadius() / this->getFluidDynamicViscosity();
                double Re_p = 1.0 * this->getFluidDensity() * relVelMagnitude * 2.0 * p_->getRadius() / this->getFluidDynamicViscosity();
    
                if (Re_p != 0)
                {
                    logger(DEBUG, "Re_p = %, relVel = (% % %), rho_f = %", Re_p, relVel.X, relVel.Y, relVel.Z,
                           this->getFluidDensity());
                    
                    logger.assert_always(Re_p < 1.0e6, "Computation of Cd is not valid for Re_p = % > 1.0e6",Re_p);
                    //https://www.sciencedirect.com/science/article/pii/S1001627917302901?via%3Dihub#bib47 (Valid to Re_p < 1e6)
                    // 24./Re + (4.2798.*(0.2.*Re).^0.2374)./(1+(0.2.*Re).^0.8961) + (0.4769.*((Re/263000).^(-7.5968)))./(1+((Re./263000).^(-7.6707))) + (Re.^0.8156)./461000
                    double Cd = 24.0 / Re_p + (4.2798 * pow(Re_p / 5.0, 0.2374)) / (1.0 + pow(Re_p / 5.0, 0.8961)) +
                                (0.4769 * pow(Re_p / 263000, 0. - 7.5968)) / (1.0 + pow(Re_p / 263000, 0. - 7.6707)) +
                                (pow(Re_p, 0.8156)) / (461000);
    
                    // ------------------------------------------------- Stokes law ------------------------------------------------- //
                    // compute drag force on particle
                    Vec3D dF_stokes = {-(relVel.X),
                                       -(relVel.Y),
                                       -(relVel.Z)};
                    dF_stokes *= 1. / 2. * Cd * relVelMagnitude * this->getFluidDensity() * constants::pi * p_->getRadius() * p_->getRadius();
                    logger(DEBUG, "Using Stokes, Fd(iPart = %) = (%, %, %)",p_->getId(), dF_stokes.X,dF_stokes.Y,dF_stokes.Z);
    
                    p_->addForce(this->getFractionDragForce() * dF_stokes);
    
                    // ---------------------------------------------- Ergun + Wen-Yu ------------------------------------------------ //
                    double beta;
                    if (Re_p < 1e3)
                    {
                        Cd = 24./Re_p * (1.0+0.15*pow(Re_p,0.687));
                    }
                    else
                    {
                        Cd = 0.44;
                    }
                    
                    // FIXME voidInEl is not defined for NS-elements
                    /*
                    if (voidInEl <= 0.8) // Ergun
                    {
                        beta = 150.0 * (1.0-voidInEl)*(1.0-voidInEl)/voidInEl * getFluidDynamicViscosity()/(4.0*pRad*pRad) +
                               1.75*(1-voidInEl)*getFluidDensity()/(2.0*pRad)*relVelMagnitude;
                    }
                    else //Wen-Yu
                    {
                        beta = 0.75*Cd*(voidInEl * (1.0-voidInEl))/(2.0*pRad)*getFluidDensity()*relVelMagnitude*pow(voidInEl,-2.65);
                    }
    
                    Vec3D dF = -beta/(1.0-voidInEl) * relVel * (4./3. * constants::pi * pRad * pRad * pRad);
                    dF /= getBodyForceScaling();
    
                    if (iPart == 0)
                    {
                        logger(DEBUG, "Re_p = %", Re_p);
                        logger(DEBUG, "beta = %, voidInEl = %", beta, voidInEl);
                        logger(DEBUG, "relVelMagnitude = %, Cd = %", relVelMagnitude, Cd);
                        logger(DEBUG, "Using Ergun and Wen-Yu, Fd(iPart = %) = (%, %, %)",iPart,dF.X,dF.Y,dF.Z);
                    }
                    */
                }
            }
        }
    }
    
    
    void computePressureGradient(ELEMENT* elPtr_)
    {
        if (this->getFractionPressureGradient() <= 0.0)
        {
            logger(DEBUG,"Pressure gradient is not activated");
            return;
        }
    
        logger(DEBUG,"In computeDragForce()");
        oomph::Vector<double> s(3,0.0);
        GeomObject* geom_obj_pt = nullptr;
        oomph::Vector<double> velAtPos(3,0.0);
    
        for (int pId = 0; pId < this->particleHandler.getNumberOfObjects(); pId++)
        {
            auto p_ = this->particleHandler.getObject(pId);
            double rad = p_->getRadius();
        
            oomph::Vector<double> s(3, 0.0);
            GeomObject* geom_obj_pt = nullptr;
            elPtr_->locate_zeta(convertVecFuncs::convertToOomphVec(this->particleHandler.getObject(pId)->getPosition()), geom_obj_pt, s);
        
            if (geom_obj_pt != nullptr)
            {
                oomph::Vector<double> dp(3, 0.0);
    
                ///FIXME THIS FUNCTION DOES NOT EXIST IN GENERAL NAVIER STOKES ELEMENTS
                //dp[0] = elPtr_->interpolated_dpdx_nst(s,0);
                //dp[1] = elPtr_->interpolated_dpdx_nst(s,1);
                //dp[2] = elPtr_->interpolated_dpdx_nst(s,2);
    
                Vec3D pGrad = convertVecFuncs::convertToVec3D(dp);
                
                ///FIXME for now voidInEl is set to 1
                //Vec3D pForceOnP = -p_->getVolume()*(1.-voidInEl) * pGrad;
                Vec3D pForceOnP = -p_->getVolume()*(1.-1.0) * pGrad;
            }
        }
    }
    
    void computeShearForce(ELEMENT* elPtr_)
    {
        if (this->getFractionShearForce() <= 0.0)
        {
            logger(DEBUG,"Shear force is not activated");
            return;
        }
    }
    
    void computeMagnusLift(ELEMENT* elPtr_)
    {
        if (this->getFractionMagnusLift() <= 0.0)
        {
            logger(DEBUG,"Magnus lift is not activated");
            return;
        }
    }
    
    void computeSaffmanLift(ELEMENT* elPtr_)
    {
        if (this->getFractionSaffmanLift() <= 0.0)
        {
            logger(DEBUG,"Saffmann lift is not activated");
            return;
        }
    }
    
};

