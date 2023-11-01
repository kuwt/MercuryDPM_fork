//
// Created by mitchel on 17/Jan/19.
//

// This driver code will contain a 3D system that is used for the simulation of a die filling process.
#include "DieFilling.h"

int main(int argc, char**argv)
{
    double domainLength = 2.5e-2; // Enter in m (x-dir)
    double domainWidth = 2.5e-2; // Enter in m (y-dir)
    double domainDepth = 10e-2; // Enter in m (z-dir)
    
    const int nx = 5;
    const int ny = 5;
    const int nz = 15;
    
    DieFilling<oomph::RefineableAJQCrouzeixRaviartElement<3>> problem(0.0, domainLength, 0.0, domainWidth, 0.0, domainDepth,nx,ny,nz);

    problem.setInflowVel(0.000); // Enter in m/s
    //problem.setFluidDensity(1.225);
    //problem.setFluidDynamicViscosity(1.81e-5);
    //problem.setFluidKinematicViscosity(1.470e-5);
    
    problem.setFluidDensity(1000.0);
    problem.setFluidDynamicViscosity(8.91e-4);
    problem.setFluidKinematicViscosity(8.91e-4);
    
    double computedReynolds = domainDepth * problem.getInflowVel() * problem.getFluidDensity()/problem.getFluidDynamicViscosity();
    logger(INFO,"computed Reynolds number = %",computedReynolds);
    
    Mdouble scaleFactorPSizeMCC = 10.0;
    Mdouble dMin = 4e-3 * scaleFactorPSizeMCC;
    Mdouble g = 9.81;
    
    Mdouble tg = sqrt(dMin/g);
    Mdouble tc = tg/1000; ///\todo

    double dt_Merc = tc/50.0;
    double dt_Oomph = 1000.*dt_Merc;

    problem.setReynoldsNumber(computedReynolds);
    problem.scaleReynolds();
    problem.setReynoldsStrouhalNumber(computedReynolds);
    problem.scaleReynoldsStrouhal();
    problem.setReynolds_InverseFroudeNumber(computedReynolds);
    problem.scaleReynolds_InverseFroude();

    problem.setSettling(false);
    problem.setAdaptOn(false);
    problem.setAdaptEveryNFluidTimesteps(10);
    problem.setUpdateCouplingEveryNParticleTimesteps(100);
    problem.setInteractionForceFd(true);

    problem.setTimeStep(dt_Merc);
    problem.setTimeStepOomph(dt_Oomph);
    problem.setTimeMax(0.15);

    problem.setSaveCount(10);// Store data every so many mercury timesteps
    problem.setSaveCountOomph(10); // Store data every so many oomph timesteps
    
    problem.setGravity(Vec3D(0.0, 0.0, -g));
    problem.dataFile.setFileType(FileType::ONE_FILE);
    problem.eneFile.setFileType(FileType::NO_FILE);
    problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.setParticlesWriteVTK(false);
    problem.wallHandler.setWriteVTK(false);
    problem.setName("singleParticleTest/MercSol");
    problem.removeOldFiles();

    oomph::DocInfo doc_info;
    doc_info.set_directory("singleParticleTest");
    doc_info.number() = 0;

    problem.setupInitialConditions();
    problem.generateLists();

    //FIXME steady state newton not working due to default of voidage and body force from AJeq
    //problem.steady_newton_solve();
    //problem.doc_solution(doc_info);
    //problem.doc_voidage(doc_info);
    //problem.doc_element(doc_info);
    //doc_info.number()++;

    problem.solveSystem(doc_info);
    problem.doc_solution(doc_info);
    problem.doc_voidage(doc_info);
    problem.doc_element(doc_info);
    doc_info.number()++;
    
    return 0;
}





