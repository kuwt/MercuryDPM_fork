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

#ifndef GRAPHICSDRAWINGAREA_H
#define GRAPHICSDRAWINGAREA_H

#include <gtkmm/drawingarea.h>

#include "DPMBase.h"
#include "BaseFunctions.h"

class MyWindow;

class GraphicsDrawingArea : public Gtk::DrawingArea
{
public:
    GraphicsDrawingArea(MyWindow* mW_);
    virtual ~GraphicsDrawingArea();
	void Xflip();
	void Yflip();
    void setColorBy(ColorBy colorBy);
    void make_squared();
    void setXDrawDir(int xDrawDir_){xDrawDir=xDrawDir_;}
    void setYDrawDir(int yDrawDir_){yDrawDir=yDrawDir_;}
    const Vec3D& getDomainMin(){return domainMin;}
    const Vec3D& getDomainMax(){return domainMax;}
    void setDomainMin(Vec3D domainMin_){domainMin=domainMin_;}
    void setDomainMax(Vec3D domainMax_){domainMax=domainMax_;}
    void resetView();
protected:
    double getColorRed  (double x, double xMin, double xMax) const;
    double getColorGreen(double x, double xMin, double xMax) const;
    double getColorBlue (double x, double xMin, double xMax) const;
	virtual bool on_draw(const Cairo::RefPtr<Cairo::Context>& cr);
	virtual bool on_scroll_event(GdkEventScroll* event);
    virtual bool on_button_press_event(GdkEventButton* event);
    virtual bool on_button_release_event(GdkEventButton* event);
	double coord2pix(double x,double xmin, double xmax, int N);
	double pix2coord(double x,double xmin, double xmax, int N);
	double length2pix(double x,double xmin, double xmax, int N);
	double pix2length(double x,double xmin, double xmax, int N);
	void Zoom(double x, double y, double zoom);

    MyWindow *mW;
	double xClick,yClick;
    bool rightMousePressed;
	Vec3D domainMin, domainMax;
	int xDrawDir, yDrawDir;
    ColorBy colorBy;
};

#endif
