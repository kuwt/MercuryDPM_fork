//
// Created by irana on 3/2/18.
//

#ifndef MERCURYDPM_PERIODICCHUTE_H
#define MERCURYDPM_PERIODICCHUTE_H

#include <fstream>
#include "Mercury3D.h"


class PeriodicChute : public Mercury3D
{
public:
    
    PeriodicChute(std::string roughBottomFile, Mdouble height, bool isMonodisperse);
    
    //Set up periodic walls, rough bottom, add flow particles
    void setupInitialConditions() override;
    
    void insertParticles();
    
    //add flow particles
    void insertParticles(unsigned int numberOfLargeToGo, unsigned int numberOfSmallToGo);
    
    //defines type of flow particles
    void insertLargeParticle(unsigned int& numberOfLargeToGo);
    
    void insertSmallParticle(unsigned int& numberOfSmallToGo);
    
    void insertOneParticle(BaseParticle& p0);
    
    virtual void setBoundaries();
    
    void setChuteProperties();
    
    void setTimeAndSaveCount();
    
    virtual void setSpeciesProperties();
    
    void printTime() const override
    {
        logger(INFO, "t = %, tmax = %", getTime(), getTimeMax());
    }
    
    void actionsAfterTimeStep() override;
    
    
private:
    unsigned int numberOfParticles;
    unsigned int numberOfSmall;
    unsigned int numberOfLarge;
    std::ofstream comFile;
    bool isMonoDisperse;
};

#endif //MERCURYDPM_PERIODICCHUTE_H
