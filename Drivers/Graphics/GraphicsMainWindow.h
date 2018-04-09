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

#ifndef GRAPHICSMAINWINDOW_H
#define GRAPHICSMAINWINDOW_H

#include <gtkmm/window.h>
#include <gtkmm/uimanager.h>
#include <gtkmm/radioaction.h>
#include <gtkmm/stock.h>
#include <gtkmm/filechooserdialog.h>
#include <vector>

#include "DPMBase.h"
#include "GraphicsDrawingArea.h"
#include "ChangeDomainDialog.h"
#include "ParticleInfoDialog.h"

class MyWindow : public Gtk::Window
{
public:
    MyWindow();
    virtual ~MyWindow();
    DPMBase* problem;
    std::vector<bool> isParticleSelected;
    ParticleInfoDialog* particleInfoDialog;
    
private:
    void on_FileQuit();
    void on_ViewXdirDrawDirX();
    void on_ViewXdirDrawDirY();
    void on_ViewXdirDrawDirZ();
	void on_ViewXDirInvert();
    void on_ViewYdirDrawDirX();
    void on_ViewYdirDrawDirY();
    void on_ViewYdirDrawDirZ();    
	void on_ViewYDirInvert();
    void on_ViewColorBy(ColorBy colorBy);
	void on_ViewMakeSquared();
	void on_ViewReset();    
    void on_ViewChangeDomain();
    void on_ViewParticleInfoDialog();
    void on_FileOpen();
    bool on_key_press_event(GdkEventKey* event);
	
	GraphicsDrawingArea* area;
	Gtk::VBox m_Box;

	Glib::RefPtr<Gtk::UIManager> m_refUIManager;
	Glib::RefPtr<Gtk::ActionGroup> m_refActionGroup;	
};
#endif
