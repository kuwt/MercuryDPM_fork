//
// Created by reza on 04-01-21.
//

#include "PIController.h"

void PIControllerBasic::set(Mdouble pGain, Mdouble iGain)
{
    this->pGain_ = pGain;
    this->iGain_ = iGain;
    //logger(INFO,"Set PI controller with pGain %, iGain %", pGain, iGain);
}

/**
 * Applies a PI-controller, strain = pGain * error + iGain * interror,
 * with the integrated stress interror += error * timeStep
 * @param error =  (stress-stressGoal)
 * @param timeStep
 * @return strain
 */
Mdouble PIControllerBasic::apply (Mdouble error, Mdouble timeStep)  {
    // Proportional Controller
    Mdouble pController = pGain_ * error;

    //Integral Controller
    iError_ += iGain_ * timeStep * error;

    // Controller Command
    return pController + iError_;
}

void PIControllerBasic::reset () {
    iError_ = 0;
}

Mdouble PIController::apply (Mdouble error, Mdouble timeStep)  {
    sumErrorSquared_ += error*error;
    nApplied_++;
    return PIControllerBasic::apply (error, timeStep);
}

void PIController::reset () {
    sumErrorSquared_ = 0.0;
    nApplied_ = 0;
    PIControllerBasic::reset();
}

double PIController::getErrorVariance() const {
    return sumErrorSquared_/nApplied_;
}

unsigned PIController::getNApplied() const {
    return nApplied_;
}