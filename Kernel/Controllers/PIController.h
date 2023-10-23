//
// Created by reza on 04-01-21.
//

#ifndef MERCURYDPM_PICONTROLLER_H
#define MERCURYDPM_PICONTROLLER_H
#include "Math/ExtendedMath.h"

/*
 * Applies a PI-controller,
 *   control = pGain * error + iGain * intError,
 * with
 *   error = processVariable - setPoint
 * and integrated error
 *   intError += error * timeStep.
 *
 * See https://en.wikipedia.org/wiki/Proportional%E2%80%93integral%E2%80%93derivative_controller
 */
class PIControllerBasic {
    // the gain factor applied to the error
    Mdouble pGain_ = 0.0;
    // the gain factor applied the integrated error
    Mdouble iGain_ = 0.0;
    // the integrated error
    Mdouble iError_ = 0.0;

public:
    // default constructor
    PIControllerBasic () = default;
    // constructor
    PIControllerBasic(Mdouble pGain, Mdouble iGain) : pGain_(pGain), iGain_(iGain) {}
    // sets pGain and iGain
    void set(Mdouble pGain, Mdouble iGain);
    // resets iError to zero
    void reset ();
    // returns the value of the control variable
    Mdouble  apply (Mdouble error, Mdouble timeStep);
    // returns iError
    double getIError() const { return iError_; }
};

/*
 * PI controller with error estimate
 */
class PIController : public PIControllerBasic {
    // the integrated error
    Mdouble sumErrorSquared_ = 0.0;
    unsigned nApplied_ = 0;

public:
    // default constructor
    PIController () = default;
    // constructor
    PIController(Mdouble pGain, Mdouble iGain) : PIControllerBasic(pGain, iGain) {}
    // resets iError to zero
    void reset ();
    // returns the value of the control variable
    Mdouble  apply (Mdouble error, Mdouble timeStep);
    // returns variance of error
    double getErrorVariance() const;
    // returns variance of error
    unsigned getNApplied() const;
};

#endif //MERCURYDPM_PICONTROLLER_H
