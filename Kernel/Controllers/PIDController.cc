//
// Created by reza on 04-01-21.
//
#include "PIDController.h"
Mdouble PIDController::apply (Mdouble stressError, Mdouble timeStep) {

    Mdouble filterCoefficient = 0.2;
    if (filterFlag==true) {
        filterError = stressError;
        filterFlag=false;
    }
    else {
        filterError = filterCoefficient * filterError +(1-filterCoefficient) * stressError;
    }

    // Proportional Controller
    Mdouble pController = pGain * stressError;

    //Integral Controller
    iController += iGain * timeStep * stressError;

    //Derivative Controller
    Mdouble dController=dGain*(stressError-previousError)/timeStep;
    previousError=stressError;

    // Controller Command
    return pController + iController+dController;
}

