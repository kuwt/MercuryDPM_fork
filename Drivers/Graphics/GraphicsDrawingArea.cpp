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

#include <limits>
#include "GraphicsDrawingArea.h"
#include "GraphicsMainWindow.h"

GraphicsDrawingArea::GraphicsDrawingArea(MyWindow* mW_)
{
	mW=mW_;
    xDrawDir=0;
	yDrawDir=1;
    resetView();
	add_events(Gdk::SCROLL_MASK);
    add_events(Gdk::BUTTON_PRESS_MASK);
    add_events(Gdk::BUTTON_RELEASE_MASK);
}

GraphicsDrawingArea::~GraphicsDrawingArea()
{
}

bool GraphicsDrawingArea::on_draw(const Cairo::RefPtr<Cairo::Context>& cr)
{
    //std::cout<<"on_draw"<<std::endl;
    if(mW->problem)
    {
        Gtk::Allocation allocation = get_allocation();
        const int width = allocation.getCGWidthidth();
        const int height = allocation.get_height();
        double minColor,maxColor;
        switch(colorBy)
        {
        case XPos:
            minColor=mW->problem->particleHandler.getLowestPositionComponentParticle(0)->getPosition().X;
            maxColor=mW->problem->particleHandler.getHighestPositionComponentParticle(0)->getPosition().X;
            break;
        case YPos:
            minColor=mW->problem->particleHandler.getLowestPositionComponentParticle(1)->getPosition().Y;
            maxColor=mW->problem->particleHandler.getHighestPositionComponentParticle(1)->getPosition().Y;
            break;
        case ZPos:
            minColor=mW->problem->particleHandler.getLowestPositionComponentParticle(2)->getPosition().Z;
            maxColor=mW->problem->particleHandler.getHighestPositionComponentParticle(2)->getPosition().Z;
            break;
        case XVel:
            minColor=mW->problem->particleHandler.getLowestVelocityComponentParticle(0)->getVelocity().X;
            maxColor=mW->problem->particleHandler.getHighestVelocityComponentParticle(0)->getVelocity().X;
            break;
        case YVel:
            minColor=mW->problem->particleHandler.getLowestVelocityComponentParticle(1)->getVelocity().Y;
            maxColor=mW->problem->particleHandler.getHighestVelocityComponentParticle(1)->getVelocity().Y;
            break;
        case ZVel:
            minColor=mW->problem->particleHandler.getLowestVelocityComponentParticle(2)->getVelocity().Z;
            maxColor=mW->problem->particleHandler.getHighestVelocityComponentParticle(2)->getVelocity().Z;
            break;        
        case AbsVel:
            minColor=0;
            maxColor=mW->problem->particleHandler.getFastestParticle()->getVelocity().getLength();
            break;
        default:
            minColor=0;
            maxColor=1;
            break;
        }   
          
        for (std::vector<BaseParticle*>::const_iterator it=mW->problem->particleHandler.begin();it!=mW->problem->particleHandler.end();it++)
        {
            double xp=coord2pix((*it)->getPosition().getComponent(xDrawDir),domainMin.getComponent(xDrawDir),domainMax.getComponent(xDrawDir),width);
            double yp=coord2pix((*it)->getPosition().getComponent(yDrawDir),domainMin.getComponent(yDrawDir),domainMax.getComponent(yDrawDir),height);
            double rx=fabs(length2pix((*it)->getRadius(),domainMin.getComponent(xDrawDir),domainMax.getComponent(xDrawDir),width));
            double ry=fabs(length2pix((*it)->getRadius(),domainMin.getComponent(yDrawDir),domainMax.getComponent(yDrawDir),height));
            double color;
            
            switch(colorBy)
                {
                case XPos:
                    color=(*it)->getPosition().X;
                    break;
                case YPos:
                    color=(*it)->getPosition().Y;
                    break;
                case ZPos:
                    color=(*it)->getPosition().Z;
                    break;
                case XVel:
                    color=(*it)->getVelocity().X;
                    break;
                case YVel:
                    color=(*it)->getVelocity().Y;
                    break;
                case ZVel:
                    color=(*it)->getVelocity().Z;
                    break;                
                case AbsVel:
                    color=(*it)->getVelocity().getLength();
                    break;
                default:
                    color=0.5;
                    break;                
                }        

            //std::cout<<"Tyring to draw ellipse at ("<<xp<<","<<yp<<") width radii ("<<rx<<","<<ry<<")"<<std::endl;
            if(rx>0&&ry>0)
            {
                cr->save();
                cr->translate(xp, yp);
                cr->scale(rx, ry);
                cr->set_source_rgb(getColorRed(color,minColor,maxColor),getColorGreen(color,minColor,maxColor),getColorBlue(color,minColor,maxColor));
                cr->arc(0.0, 0.0, 1.0, 0.0, 2 * constants::pi);
                cr->fill();
                if(mW->isParticleSelected[(*it)->getIndex()])
                {
                    std::cout<<"In draw Selected "<<**it<<" index="<<(*it)->getIndex()<<std::endl;
                    cr->set_source_rgb(0.0,0.0,0.0);
                    cr->set_line_width(0.2);
                    cr->arc(0.0, 0.0, 1.0, 0.0, 2 * constants::pi);
                    cr->restore();
                    cr->stroke();
                }
                else
                {               
                    cr->restore();
                }
            }
        }
    }
 //   std::cout<<"Finished draw"<<std::endl;
    return true;
}

bool GraphicsDrawingArea::on_scroll_event(GdkEventScroll* event)
{
//	std::cout<<"MyArea::on_scroll_event"<<std::endl;
	if (event->direction == GDK_SCROLL_UP) 
	{ 
//		std::cout<<"Scroll up at x="<<event->x<<" y="<<event->y<<std::endl;
		Zoom(event->x,event->y,0.5);
        queue_draw();
        return true;
	} 
	else if (event->direction == GDK_SCROLL_DOWN) 
	{ 
//		std::cout<<"Scroll down at x="<<event->x<<" y="<<event->y<<std::endl;
		Zoom(event->x,event->y,2);
        queue_draw();
        return true;
	} 
	return false;
}

bool GraphicsDrawingArea::on_button_press_event(GdkEventButton* event)
{
//    std::cout<<"MyArea::on_button_press_event(GdkEventButton* event)"<<event<<std::endl;
    if(event->button == 1)
    {
        double dx =pix2coord(event->x,domainMin.getComponent(xDrawDir),domainMax.getComponent(xDrawDir),get_allocation().getCGWidthidth());
        double dy =pix2coord(event->y,domainMin.getComponent(yDrawDir),domainMax.getComponent(yDrawDir),get_allocation().get_height());
        double distance=std::numeric_limits<double>::max();
        double curdistance;
        BaseParticle *P=NULL;
        for (std::vector<BaseParticle*>::const_iterator it=mW->problem->particleHandler.begin();it!=mW->problem->particleHandler.end();it++)
        {
            curdistance=pow((*it)->getPosition().getComponent(xDrawDir)-dx,2)+pow((*it)->getPosition().getComponent(yDrawDir)-dy,2);
            if (curdistance<distance)
            {
                distance=curdistance;
                P=*it;
            }
        }
        if(P)    
        {
            if(!mW->isParticleSelected[P->getIndex()])
            {
                mW->isParticleSelected[P->getIndex()]=true;
                ParticleInfoTab* test;
                test = new ParticleInfoTab(P->getIndex(),mW->problem);
                mW->particleInfoDialog->particleInfoTabVector.push_back(test);
                std::string name="Particle: ";
                name.append(ToString(P->getIndex()));
                mW->particleInfoDialog->particleInfoTabs.append_page(*test, name);
                mW->particleInfoDialog->particleInfoTabs.show_all();  
            }
            else
            {
                mW->isParticleSelected[P->getIndex()]=false;
                for (std::vector<ParticleInfoTab*>::iterator it=mW->particleInfoDialog->particleInfoTabVector.begin();it!=mW->particleInfoDialog->particleInfoTabVector.end();it++)
                {
                    if((*it)->getParticleIndex()==P->getIndex())
                    {
                        mW->particleInfoDialog->particleInfoTabs.remove_page(std::distance(mW->particleInfoDialog->particleInfoTabVector.begin(), it));
                        mW->particleInfoDialog->particleInfoTabVector.erase(it);
                        break;
                    }
                }
                
            }
        }
        queue_draw();
        return true;
    }    
    if(event->button == 3)
    {
        rightMousePressed=true;
        xClick=event->x;
        yClick=event->y;
        return true;
    }
    return false;
}

bool GraphicsDrawingArea::on_button_release_event(GdkEventButton* event)
{
 //   std::cout<<"MyArea::on_button_release_event(GdkEventButton* event)"<<event<<std::endl;
    if(event->button == 3)
    {
        const int width = get_allocation().getCGWidthidth();
        const int height = get_allocation().get_height();
        double dx=pix2length(event->x-xClick,domainMin.getComponent(xDrawDir),domainMax.getComponent(xDrawDir),width);
        double dy=pix2length(event->y-yClick,domainMin.getComponent(yDrawDir),domainMax.getComponent(yDrawDir),height);
        domainMin.setComponent(xDrawDir,domainMin.getComponent(xDrawDir)-dx);
        domainMax.setComponent(xDrawDir,domainMax.getComponent(xDrawDir)-dx);
        domainMin.setComponent(yDrawDir,domainMin.getComponent(yDrawDir)-dy);
        domainMax.setComponent(yDrawDir,domainMax.getComponent(yDrawDir)-dy);
        rightMousePressed=false;
        queue_draw();
        return true;
    }
    return false;
}

void GraphicsDrawingArea::Zoom(double x, double y, double zoom)
{
//	cout<<"Zoom met x="<<x<<" en y="<<y<<" en zoom="<<zoom<<endl;
	Gtk::Allocation allocation = get_allocation();
	const int width = allocation.getCGWidthidth();
	const int height = allocation.get_height();
	
	double dx =pix2coord(x, domainMin.getComponent(xDrawDir), domainMax.getComponent(xDrawDir), width);
	double dy =pix2coord(y, domainMin.getComponent(yDrawDir), domainMax.getComponent(yDrawDir), height);
//	cout<<"Old domain ["<<domainMin.getComponent(xDrawDir)<<","<<domainMax.getComponent(xDrawDir)<<"]x["<<domainMin.getComponent(yDrawDir)<<","<<domainMax.getComponent(yDrawDir)<<"]"<<endl;
	domainMin.setComponent(xDrawDir,zoom*(domainMin.getComponent(xDrawDir)-dx)+dx);
	domainMax.setComponent(xDrawDir,zoom*(domainMax.getComponent(xDrawDir)-dx)+dx);
	domainMin.setComponent(yDrawDir,zoom*(domainMin.getComponent(yDrawDir)-dy)+dy);
	domainMax.setComponent(yDrawDir,zoom*(domainMax.getComponent(yDrawDir)-dy)+dy);
//	cout<<"New domain ["<<domainMin.getComponent(xDrawDir)<<","<<domainMax.getComponent(xDrawDir)<<"]x["<<domainMin.getComponent(yDrawDir)<<","<<domainMax.getComponent(yDrawDir)<<"]"<<endl;
}

void GraphicsDrawingArea::make_squared()
{
//	cout<<"make_squared"<<endl;
	Gtk::Allocation allocation = get_allocation();
	const int width = allocation.getCGWidthidth();
	const int height = allocation.get_height();
	double Xres=(domainMax.getComponent(xDrawDir)-domainMin.getComponent(xDrawDir))/width;
	double Yres=(domainMax.getComponent(yDrawDir)-domainMin.getComponent(yDrawDir))/height;
	if(Xres>Yres) //Change Ymax and Ymin
	{
		double Mean=0.5*(domainMax.getComponent(yDrawDir)+domainMin.getComponent(yDrawDir));
		domainMin.setComponent(yDrawDir,Mean-0.5*height*Xres);
		domainMax.setComponent(yDrawDir,Mean+0.5*height*Xres);
	}
	else //Change Xmax and Xmin
	{
		double Mean=0.5*(domainMax.getComponent(xDrawDir)+domainMin.getComponent(xDrawDir));
		domainMin.setComponent(xDrawDir,Mean-0.5*width*Yres);
		domainMax.setComponent(xDrawDir,Mean+0.5*width*Yres);
	}
}

double GraphicsDrawingArea::coord2pix (double x,double xmin, double xmax, int N){return (x-xmin)/(xmax-xmin)*(N-1);}
double GraphicsDrawingArea::pix2coord (double x,double xmin, double xmax, int N){return xmin+x*(xmax-xmin)/(N-1);}
double GraphicsDrawingArea::length2pix(double x,double xmin, double xmax, int N){return x/(xmax-xmin)*(N-1);}
double GraphicsDrawingArea::pix2length(double x,double xmin, double xmax, int N){return x*(xmax-xmin)/(N-1);}

void GraphicsDrawingArea::Xflip()
{
//	cout<<"Old domain ["<<domainMin.getComponent(xDrawDir)<<","<<domainMax.getComponent(xDrawDir)<<"]x["<<domainMin.getComponent(yDrawDir)<<","<<domainMax.getComponent(yDrawDir)<<"]"<<endl;
    double dummy=domainMin.getComponent(xDrawDir);
    domainMin.setComponent(xDrawDir,domainMax.getComponent(xDrawDir));
    domainMax.setComponent(xDrawDir,dummy);
//	cout<<"New domain ["<<domainMin.getComponent(xDrawDir)<<","<<domainMax.getComponent(xDrawDir)<<"]x["<<domainMin.getComponent(yDrawDir)<<","<<domainMax.getComponent(yDrawDir)<<"]"<<endl;
}
void GraphicsDrawingArea::Yflip()
{
//	cout<<"Old domain ["<<domainMin.getComponent(xDrawDir)<<","<<domainMax.getComponent(xDrawDir)<<"]x["<<domainMin.getComponent(yDrawDir)<<","<<domainMax.getComponent(yDrawDir)<<"]"<<endl;
    double dummy=domainMin.getComponent(yDrawDir);
    domainMin.setComponent(yDrawDir,domainMax.getComponent(yDrawDir));
    domainMax.setComponent(yDrawDir,dummy);
//	cout<<"New domain ["<<domainMin.getComponent(xDrawDir)<<","<<domainMax.getComponent(xDrawDir)<<"]x["<<domainMin.getComponent(yDrawDir)<<","<<domainMax.getComponent(yDrawDir)<<"]"<<endl;
}

void GraphicsDrawingArea::setColorBy(ColorBy colorBy)
{
    this->colorBy=colorBy;    
//    std::cout<<"Setting color"<<std::endl;
}

void GraphicsDrawingArea::resetView()
{
//    std::cout<<"Reset view"<<std::endl;
    if(mW->problem)
    {
        domainMin.X=mW->problem->getXMin();
        domainMin.Y=mW->problem->getYMin();
        domainMin.Z=mW->problem->getZMin();
        domainMax.X=mW->problem->getXMax();
        domainMax.Y=mW->problem->getYMax();
        domainMax.Z=mW->problem->getZMax();
    }
}

double GraphicsDrawingArea::getColorRed  (double x, double xMin, double xMax) const
{
    if(xMax<=xMin) return 0.5;
    if     (8.0*(x-xMin)<3.0*(xMax-xMin)) return  0.0;
    else if(8.0*(x-xMin)<5.0*(xMax-xMin)) return -1.5+4.0*x/xMax;
    else if(8.0*(x-xMin)<7.0*(xMax-xMin)) return  1.0;
    else                                  return  4.5-4.0*x/xMax;
}
double GraphicsDrawingArea::getColorGreen(double x, double xMin, double xMax) const
{
    if(xMax<=xMin) return 1.0;
    if     (8.0*(x-xMin)<1.0*(xMax-xMin)) return  0.0;
    else if(8.0*(x-xMin)<3.0*(xMax-xMin)) return -0.5+4.0*x/xMax;
    else if(8.0*(x-xMin)<5.0*(xMax-xMin)) return  1.0;
    else if(8.0*(x-xMin)<7.0*(xMax-xMin)) return  3.5-4.0*x/xMax;
    else                                  return  0.0;
}
double GraphicsDrawingArea::getColorBlue (double x, double xMin, double xMax) const
{
    if(xMax<=xMin) return 0.5;
    if     (8.0*(x-xMin)<1.0*(xMax-xMin)) return  0.5+4.0*x/xMax;
    else if(8.0*(x-xMin)<3.0*(xMax-xMin)) return  1.0;
    else if(8.0*(x-xMin)<5.0*(xMax-xMin)) return  2.5-4.0*x/xMax;
    else                                  return  0.0;
}
