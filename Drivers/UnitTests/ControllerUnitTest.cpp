//
// Created by reza on 21-01-21.
//

#include <Math/Matrix.h>
#include "Controllers/PController.h"
#include "Controllers/PIController.h"
#include "Controllers/PIDController.h"
#include "Boundaries/StressStrainControlBoundary.h"
class ControllerUnitTest{
    Mdouble x1=0.0;
    Mdouble x2=0.0;
    Mdouble m=4;
    Mdouble c=10;
    Mdouble k=100.0;
    Mdouble t=0.0;
    Mdouble dt=0.01;
    Mdouble xDesire=5;
    Mdouble F=0.0;
    Mdouble x1Error=0.0;
public:
    ControllerUnitTest () {
    }



    Mdouble f1(Mdouble x2){
        Mdouble x1dot=x2;
        return x1dot;
    }

    Mdouble f2(Mdouble x1,Mdouble x2){
        Mdouble x2dot=1/m*(F-c*x2-k*x1);
        return x2dot;
    }

    Mdouble RungeKutta(Mdouble x1,Mdouble x2,Mdouble dt){
        PIController xx(0.1,0.05);
        for (int i=0;i<200/dt;i++) {
            Mdouble k11 = dt*f1(x2);
            Mdouble k21 = dt*f2(x1,x2);
            Mdouble k12 = dt*f1(x2+0.5*k21);
            Mdouble k22 = dt*f2(x1+0.5*k11,x2+0.5*k21);
            Mdouble k13 = dt*f1(x2+0.5*k22);
            Mdouble k23 = dt*f2(x1+0.5*k12,x2+0.5*k22);
            Mdouble k14 = dt*f1(x2+k23);
            Mdouble k24 = dt*f2(x1+k13,x2+k23);
            x1 = x1+(k11+2*k12+2*k13+k14)/6;
            x2 = x2+(k21+2*k22+2*k23+k24)/6;
            t =t+ dt;
            x1Error=-(x1-xDesire);
            F = xx.apply(x1Error,dt);      // Controller Command
        }
        return x1Error;
    }

    Mdouble SpringMassControl(){
        Mdouble Err=RungeKutta(x1,x2,dt);
        return Err;
    }

};

int main()
{
    ControllerUnitTest springMass;
    //StressStrainControlBoundary ssError;
    //Mdouble steadyError=ssError.computeStressError();
    Mdouble springError=springMass.SpringMassControl();
    if (abs(springError)>0.1)
    std::cerr<<"The Controller is Removed"<<std::endl;

    /*if (steadyError>1.0)
        std::cout<<"Tune The Controller Gains"<<std::endl;*/
}