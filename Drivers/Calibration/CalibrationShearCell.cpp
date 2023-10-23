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
#include <Logger.h>
#include <Math/Helpers.h>
#include "Calibration.h"
using mathsFunc::square;
using helpers::toString;

// Here is a templated function to do the linear regression fit of data points.
template <typename T>
std::pair<T,T> getLinearFit(const std::vector<T> X, const std::vector<T> Y)
{
    size_t n = Y.size();
    logger.assert_always(X.size()==n,"X and Y data sets have to be same size");
    T xSum = 0, ySum = 0, xxSum = 0, xySum = 0;
    for (long i = 0; i < n; i++)
    {
        xSum += X[i];
        ySum += Y[i];
        xxSum += X[i] * X[i];
        xySum += X[i] * Y[i];
    }
    T slope = (n * xySum - xSum * ySum) / (n * xxSum - xSum * xSum);
    T intercept = (ySum - slope * xSum) / Y.size();
    logger(INFO,"Fit values: Slope %, intercept %",slope, intercept);
    return {slope, intercept};
}

//Here is a function to solve the quadratic equation
std::pair<double , double > solveQuad(double a, double b, double c)
{
    double disc = (b * b) - (4. * a * c); //The discriminant
    if(disc < 0) //If the discriminant is less than 0
    {
        logger(WARN,"Discriminant < 0. Roots are imaginary.");
        return {constants::NaN,constants::NaN};
    } else {
        disc = sqrt(disc);
        return {0.5 * (-b + disc) / a, 0.5 * (-b - disc) / a};
    }
}

struct ShearStageData {
    Mdouble compressiveStress; //input
    Mdouble maxShearStress; //output
};

/**
 * https://www.dietmar-schulze.de/grdle1.pdf
 */
double computeFFc (Mdouble preCompressionNormalStress, Mdouble preCompressionShearStress, 
                   std::vector<double> normalStress, std::vector<double> maxShearStress) {
    //This part does the linear fit through the shear points and get the yield locus with slope and intercept
    //Then we calculate the sigma_1 and sigma_c to compute flow function FFc
    logger(INFO, "p % % tau % %", preCompressionNormalStress, toString(normalStress), preCompressionShearStress,
           toString(maxShearStress));
    
    double sigmaC;
    {
        auto[slope, intercept] = getLinearFit(normalStress, maxShearStress);
        
        if (intercept < 0) {
            logger(WARN, "y-intercept % <0; setting FFc to 100", intercept);
            return 1000;
        }
        if (slope < 0) {
            logger(WARN, "slope % <0; this does not make sense; setting FFc to 0", slope);
            return 0;
        }
        //see ffc.nb
        double arclength = sqrt(1+slope*slope);
        sigmaC = 2*intercept*(arclength+slope);
    }
    
    double sigma1;
    {
        auto[slope, intercept] = getLinearFit<double>({preCompressionNormalStress, normalStress[0]}, {preCompressionShearStress, maxShearStress[0]});
        if (slope < 0) {
            sigma1 = preCompressionShearStress + preCompressionShearStress;
            logger(WARN, "slope % <0; setting sigma1 to %", slope, sigma1);
        } else {
            //see ffc.nb
            double arclength = sqrt(1 + slope * slope);
            sigma1 = preCompressionNormalStress * arclength * arclength
                     + intercept * slope
                     + arclength * (intercept + preCompressionNormalStress * slope);
        }
    }
    
    double ffc = sigma1/sigmaC;
    logger(INFO,"p = %, sigmaC = %, sigma1 = %, ffc = %",preCompressionNormalStress, sigmaC, sigma1, ffc);
    return ffc;
}

/**
 * https://www.dietmar-schulze.de/grdle1.pdf
 */
double computeSimpleFFc (Mdouble preCompressionNormalStress, Mdouble preCompressionShearStress,
                   std::vector<double> normalStress, std::vector<double> maxShearStress) {
    //This part does the linear fit through the shear points and get the yield locus with slope and intercept
    //Then we calculate the sigma_1 and sigma_c to compute flow function FFc
    
    logger(INFO, "p % % tau % %", preCompressionNormalStress, toString(normalStress), preCompressionShearStress,
           toString(maxShearStress));
    auto [slope, intercept] = getLinearFit(normalStress,maxShearStress);
    
    if (intercept < 0) {
        logger(WARN, "y-intercept % <0; setting FFc to 1000", intercept);
        return 1000;
    }
    if (slope < 0) {
        logger(WARN, "slope % <0; setting FFc to 1000", slope);
        return 1000;
    }
    
    //see ffc.nb
    double arclength = sqrt(1+slope*slope);
    double sigmaC = 2*intercept*(arclength+slope);
    double sigma1 = preCompressionNormalStress*arclength*arclength
                    +intercept*slope
                    +arclength*(intercept+preCompressionNormalStress*slope);
    double ffc = sigma1/sigmaC;
    logger(INFO,"p = %, sigmaC = %, sigma1 = %, ffc = %",preCompressionNormalStress, sigmaC, sigma1, ffc);
    return ffc;
}

//double computeSimpleFFc (Mdouble preCompressionNormalStress, Mdouble preCompressionShearStress,
//                   std::vector<double> normalStress, std::vector<double> maxShearStress) {
//    //This part does the linear fit through the shear points and get the yield locus with slope and intercept
//    //Then we calculate the sigma_1 and sigma_c to compute flow function FFc
//
//    auto [slope, intercept] = getLinearFit(normalStress,maxShearStress);
//    double ffc =  1./(intercept/preCompressionNormalStress+slope);
//    logger(INFO,"Fit values: Slope %, intercept %, ffc %",slope, intercept, ffc);
//    return ffc;
//}

int main(int argc, char* argv[])
{
    // read in pre-compression stress and subsequent compression stresses
    std::vector<double> normalStress =
        helpers::readVectorFromCommandLine<double>(argc,argv,"-normalStress",4, {});
    logger.assert_always(normalStress.size()>2,"You need at least one pre-compression stress and two compression stresses to compute an ffc");
    double preCompressionNormalStress = normalStress[0];
    normalStress.erase(normalStress.begin());
    
    //get the directory name
    std::string dir = argv[0];
    dir.resize(dir.size()-std::string("CalibrationShearCell").length());
    
    // run pre-compression stage
    {
        std::stringstream cmd;
        cmd << dir << "CalibrationPrecompression";
        for (int i = 1; i < argc; ++i) cmd << ' ' << argv[i];
        logger(INFO, "\nRunning %", cmd.str());
        int status = system(cmd.str().c_str());
    }

    // run shear stages
    for (const double s : normalStress) {
        std::stringstream cmd;
        cmd << dir << "CalibrationShearStage -compression " << s;
        for (int i = 1; i < argc; ++i) cmd << ' ' << argv[i];
        logger(INFO, "\nRunning %", cmd.str());
        int status = system(cmd.str().c_str());
    }
    
    //read back output
    Calibration dpm(argc, argv);
    dpm.setName("CalibrationShearCell" + dpm.getParam());
    Mdouble preCompressionShearStress = helpers::readFromFile(dpm.getName()+".out","preCompressionShearStress",0);
    Mdouble preCompressionBulkDensity = helpers::readFromFile(dpm.getName()+".out","preCompressionBulkDensity", 0);
    
    std::vector<Mdouble> maxShearStress;
    std::vector<Mdouble> bulkDensity;
    for (const int s : normalStress) {
        maxShearStress.push_back(
            helpers::readFromFile(dpm.getName() + "-" + std::to_string(s) + "Pa.out",
                              "maxShearStress", 0));
        bulkDensity.push_back(
                helpers::readFromFile(dpm.getName() + "-" + std::to_string(s) + "Pa.out",
                                      "bulkDensity", 0));
    }
    
    //Write output file
    std::stringstream out;
    if (helpers::readFromCommandLine(argc,argv,"-ffc")) {
        out << computeFFc(preCompressionNormalStress, preCompressionShearStress, normalStress, maxShearStress) << ' ';
        //out << computeSimpleFFc(preCompressionNormalStress, preCompressionShearStress, normalStress, maxShearStress) << ' ';
    } else {
        for (const auto s : maxShearStress)
            out << s << ' ';
        for (const auto b : bulkDensity)
            out << b << ' ';
    }
    helpers::writeToFile(dpm.getName()+".txt", out.str());
    logger(INFO,"Output to %: %",dpm.getName()+".txt", out.str());
    
    return 0;
}