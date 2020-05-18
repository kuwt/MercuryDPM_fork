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

#include <iostream>
#include <iomanip> 
#include <vector>
#include <map>

#include "Mercury3D.h"
#include "Mercury2D.h"
#include "Particles/BaseParticle.h"

class HGridCell
{
public:
    HGridCell(BaseParticle* p)
    {
        x_ = p->getHGridX();
        y_ = p->getHGridY();
        z_ = p->getHGridZ();
        l_ = p->getHGridLevel();
    }

    HGridCell(int x, int y, int z, unsigned int l)
    {
        x_ = x;
        y_ = y;
        z_ = z;
        l_ = l;
    }

    bool operator <(const HGridCell& rhs) const
            {
        if (x_ != rhs.x_)
        {
            return x_ < rhs.x_;
        }
        else if (y_ != rhs.y_)
        {
            return y_ < rhs.y_;
        }
        else if (z_ != rhs.z_)
        {
            return z_ < rhs.z_;
        }
        else
        {
            return l_ < rhs.l_;
        }
    }

    friend std::ostream& operator << (std::ostream& os, const HGridCell& rhs)
    {
        os<<rhs.x_<<" "<<rhs.y_<<" "<<rhs.z_<<" "<<rhs.l_;
        return os;
    }

//private:
    int x_, y_, z_;
    unsigned int l_;
};


class HGridCountCalls3D : public Mercury3D
{
public:
    void computeInternalForces(BaseParticle* P1, BaseParticle* P2)
    {
        BaseParticle *P1Per = P1->getPeriodicFromParticle();
        BaseParticle *P2Per = P2->getPeriodicFromParticle();
        if (P1Per != 0 && P2Per != 0)
        {
            return;
        }
        if (P1Per != 0 || P2Per != 0)
            cif[P1->getHGridLevel()][P2->getHGridLevel()]++;
        else
            cif[P1->getHGridLevel()][P2->getHGridLevel()] += 2;
        DPMBase::computeInternalForces(P1, P2);
    }
    
    void hGridFindContactsWithinTargetCell(int x, int y, int z, unsigned int l)
    {
        oif[l][l]++;
        Mercury3D::hGridFindContactsWithinTargetCell(x, y, z, l);
    }
    
    void hGridFindContactsWithTargetCell(int x, int y, int z, unsigned int l, BaseParticle* obj)
    {
        if (obj->getPeriodicFromParticle() == 0)
            oif[obj->getHGridLevel()][l]++;
        Mercury3D::hGridFindContactsWithTargetCell(x, y, z, l, obj);
    }

    void calculateParticlesPerLevel()
    {
        std::vector<int> ppl;
        std::vector<int> patl;
        ppl.resize(getHGrid()->getNumberOfLevels());
        for (unsigned int i = 0; i < ppl.size(); i++)
        {
            ppl[i] = 0;
        }
        for (std::vector<BaseParticle*>::iterator it = particleHandler.begin(); it != particleHandler.end(); ++it)
        {
            ppl[(*it)->getHGridLevel()]++;
        }
        std::cout << "Particles per level:" << std::endl;
        for (unsigned int i = 0; i < ppl.size(); i++)
        {
            std::cout << std::setw(10) << 1.0 * ppl[i] / particleHandler.getNumberOfObjects() << " ";
        }
        std::cout << std::endl;
    }

    void histNumberParticlesPerCell()
    {
        std::map<HGridCell, unsigned int> particlesPerCell;
        for (auto it = particleHandler.begin(); it != particleHandler.end(); ++it)
        {
            HGridCell cell = HGridCell(*it);
            auto obj = particlesPerCell.find(cell);
            if (obj == particlesPerCell.end())
            {
                particlesPerCell.insert(std::pair<HGridCell, unsigned int>(cell, 1));
            }
            else
            {
                obj->second++;
            }
        }
        int N=0;

        std::vector<std::vector<int> > histParticlesPerLevelPerCell;
        histParticlesPerLevelPerCell.resize(getHGridMaxLevels());
        for (auto it = particlesPerCell.begin(); it != particlesPerCell.end(); ++it)
        {
            if(histParticlesPerLevelPerCell[it->first.l_].size()<it->second+1)
            {
                histParticlesPerLevelPerCell[it->first.l_].resize(it->second+1,0);
            }
            histParticlesPerLevelPerCell[it->first.l_][it->second]++;
            N+=it->second;
            //std::cout << it->first << " " << it->second << std::endl;
        }

        for (auto it1 =histParticlesPerLevelPerCell.begin(); it1!=histParticlesPerLevelPerCell.end();++it1)
        {
            std::cout<<"Histogram for level "<<it1-histParticlesPerLevelPerCell.begin()<<": ";
            for (auto it2=it1->begin(); it2!=it1->end();++it2)
            {
                std::cout<<" "<<*it2;
            }
            std::cout<<std::endl;
        }
    }

    void perpareCalls()
    {
        std::cout << "Largest particle radius=" << particleHandler.getLargestParticle()->getRadius() << std::endl;

        cif.resize(getHGridMaxLevels());
        oif.resize(getHGridMaxLevels());
        for (unsigned int i = 0; i < getHGridMaxLevels(); i++)
        {
            cif[i].resize(getHGridMaxLevels());
            oif[i].resize(getHGridMaxLevels());
            for (unsigned int j = 0; j < getHGridMaxLevels(); j++)
            {
                cif[i][j] = 0;
                oif[i][j] = 0;
            }
        }
    }

    void displayCalls(int numberOfForceCalculations)
    {
        std::cout << "Level sizes:" << std::endl << std::setw(10) << 2.0 * particleHandler.getSmallestParticle()->getRadius();
        for (unsigned int i = 0; i < getHGrid()->getNumberOfLevels(); i++)
        {
            std::cout << std::setw(10) << getHGrid()->getCellSize(i) << " ";
        }
        std::cout << std::endl;

        calculateParticlesPerLevel();

        double totalContactWork = 0, totalOverheadWork = 0;
        for (unsigned int i = 0; i < getHGrid()->getNumberOfLevels(); i++)
        {
            for (unsigned int j = 0; j < getHGrid()->getNumberOfLevels(); j++)
            {
                totalContactWork += cif[i][j];
                totalOverheadWork += oif[i][j];
            }
        }

        std::cout << "Contactwork: " << 0.5 * totalContactWork / particleHandler.getNumberOfObjects() / numberOfForceCalculations << std::endl;
        for (unsigned int i = 0; i < getHGrid()->getNumberOfLevels(); i++)
        {
            for (unsigned int j = 0; j < getHGrid()->getNumberOfLevels(); j++)
            {
                std::cout << std::setw(10) << 0.5 * cif[i][j] / particleHandler.getNumberOfObjects() / numberOfForceCalculations << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "Overheadwork: " << totalOverheadWork / particleHandler.getNumberOfObjects() / numberOfForceCalculations << std::endl;
        for (unsigned int i = 0; i < getHGrid()->getNumberOfLevels(); i++)
        {
            for (unsigned int j = 0; j < getHGrid()->getNumberOfLevels(); j++)
            {
                std::cout << std::setw(10) << 1.0 * oif[i][j] / particleHandler.getNumberOfObjects() / numberOfForceCalculations << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "Total work: " << 0.5 * totalContactWork / particleHandler.getNumberOfObjects() / numberOfForceCalculations + totalOverheadWork / particleHandler.getNumberOfObjects() / numberOfForceCalculations << std::endl;
    }
    
    std::vector<std::vector<int> > cif;
    std::vector<std::vector<int> > oif;
};



class HGridCountCalls2D : public Mercury2D
{
public:
    void computeInternalForces(BaseParticle* P1, BaseParticle* P2)
    {
        BaseParticle *P1Per = P1->getPeriodicFromParticle();
        BaseParticle *P2Per = P2->getPeriodicFromParticle();
        if (P1Per != 0 && P2Per != 0)
        {
            return;
        }
        if (P1Per != 0 || P2Per != 0)
            cif[P1->getHGridLevel()][P2->getHGridLevel()]++;
        else
            cif[P1->getHGridLevel()][P2->getHGridLevel()] += 2;
        DPMBase::computeInternalForces(P1, P2);
    }

    void hGridFindContactsWithinTargetCell(int x, int y, unsigned int l)
    {
        oif[l][l]++;
        Mercury2D::hGridFindContactsWithinTargetCell(x, y, l);
    }
    
    void hGridFindContactsWithTargetCell(int x, int y, unsigned int l, BaseParticle* obj)
    {
        if (obj->getPeriodicFromParticle() == 0)
            oif[obj->getHGridLevel()][l]++;
        Mercury2D::hGridFindContactsWithTargetCell(x, y, l, obj);
    }
    
    void calculateParticlesPerLevel()
    {
        std::vector<int> ppl;
        std::vector<int> patl;
        ppl.resize(getHGrid()->getNumberOfLevels());
        for (unsigned int i = 0; i < ppl.size(); i++)
        {
            ppl[i] = 0;
        }
        for (std::vector<BaseParticle*>::iterator it = particleHandler.begin(); it != particleHandler.end(); ++it)
        {
            ppl[(*it)->getHGridLevel()]++;
        }
        std::cout << "Particles per level:" << std::endl;
        for (unsigned int i = 0; i < ppl.size(); i++)
        {
            std::cout << std::setw(10) << 1.0 * ppl[i] / particleHandler.getNumberOfObjects() << " ";
        }
        std::cout << std::endl;
    }
    
    void histNumberParticlesPerCell()
    {
        std::map<HGridCell, unsigned int> particlesPerCell;
        for (auto it = particleHandler.begin(); it != particleHandler.end(); ++it)
        {
            HGridCell cell = HGridCell(*it);
            auto obj = particlesPerCell.find(cell);
            if (obj == particlesPerCell.end())
            {
                particlesPerCell.insert(std::pair<HGridCell, unsigned int>(cell, 1));
            }
            else
            {
                obj->second++;
            }
        }
        int N=0;

        std::vector<std::vector<int> > histParticlesPerLevelPerCell;
        histParticlesPerLevelPerCell.resize(getHGridMaxLevels());
        for (auto it = particlesPerCell.begin(); it != particlesPerCell.end(); ++it)
        {
            if(histParticlesPerLevelPerCell[it->first.l_].size()<it->second+1)
            {
                histParticlesPerLevelPerCell[it->first.l_].resize(it->second+1,0);
            }
            histParticlesPerLevelPerCell[it->first.l_][it->second]++;
            N+=it->second;
            //std::cout << it->first << " " << it->second << std::endl;
        }

        for (auto it1 =histParticlesPerLevelPerCell.begin(); it1!=histParticlesPerLevelPerCell.end();++it1)
        {
            std::cout<<"Histogram for level "<<it1-histParticlesPerLevelPerCell.begin()<<": ";
            for (auto it2=it1->begin(); it2!=it1->end();++it2)
            {
                std::cout<<" "<<*it2;
            }
            std::cout<<std::endl;
        }
    }

    void perpareCalls()
    {
        cif.resize(getHGridMaxLevels());
        oif.resize(getHGridMaxLevels());
        for (unsigned int i = 0; i < getHGridMaxLevels(); i++)
        {
            cif[i].resize(getHGridMaxLevels());
            oif[i].resize(getHGridMaxLevels());
            for (unsigned int j = 0; j < getHGridMaxLevels(); j++)
            {
                cif[i][j] = 0;
                oif[i][j] = 0;
            }
        }
    }

    void displayCalls(int numberOfForceCalculations)
    {
        std::cout << "Level sizes:" << std::endl << std::setw(10) << 2.0 * particleHandler.getSmallestParticle()->getRadius();
        for (unsigned int i = 0; i < getHGrid()->getNumberOfLevels(); i++)
        {
            std::cout << std::setw(10) << getHGrid()->getCellSize(i) << " ";
        }
        std::cout << std::endl;

        calculateParticlesPerLevel();

        double totalContactWork = 0, totalOverheadWork = 0;
        for (unsigned int i = 0; i < getHGrid()->getNumberOfLevels(); i++)
        {
            for (unsigned int j = 0; j < getHGrid()->getNumberOfLevels(); j++)
            {
                totalContactWork += cif[i][j];
                totalOverheadWork += oif[i][j];
            }
        }

        std::cout << "Contactwork: " << 0.5 * totalContactWork / particleHandler.getNumberOfObjects() / numberOfForceCalculations << std::endl;
        for (unsigned int i = 0; i < getHGrid()->getNumberOfLevels(); i++)
        {
            for (unsigned int j = 0; j < getHGrid()->getNumberOfLevels(); j++)
            {
                std::cout << std::setw(10) << 0.5 * cif[i][j] / particleHandler.getNumberOfObjects() / numberOfForceCalculations << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "Overheadwork: " << totalOverheadWork / particleHandler.getNumberOfObjects() / numberOfForceCalculations << std::endl;
        for (unsigned int i = 0; i < getHGrid()->getNumberOfLevels(); i++)
        {
            for (unsigned int j = 0; j < getHGrid()->getNumberOfLevels(); j++)
            {
                std::cout << std::setw(10) << 1.0 * oif[i][j] / particleHandler.getNumberOfObjects() / numberOfForceCalculations << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "Total work: " << 0.5 * totalContactWork / particleHandler.getNumberOfObjects() / numberOfForceCalculations + totalOverheadWork / particleHandler.getNumberOfObjects() / numberOfForceCalculations << std::endl;
    }

    std::vector<std::vector<int> > cif;
    std::vector<std::vector<int> > oif;
};

