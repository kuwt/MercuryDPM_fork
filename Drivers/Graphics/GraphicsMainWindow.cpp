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

#include "GraphicsMainWindow.h"

MyWindow::MyWindow()
{
    problem = 0;
    area = new GraphicsDrawingArea(this);
    set_title("DrawingArea");
	set_default_size(800, 800);
	add_events(Gdk::KEY_PRESS_MASK);
    
	add(m_Box);	
        
	//Create actions for menus and toolbars:
	m_refActionGroup = Gtk::ActionGroup::create();

	//File menu:
	m_refActionGroup->add(Gtk::Action::create("FileMenu", "File"));
    m_refActionGroup->add(Gtk::Action::create("FileOpen", Gtk::Stock::OPEN),sigc::mem_fun(*this, &MyWindow::on_FileOpen));
	m_refActionGroup->add(Gtk::Action::create("FileQuit", Gtk::Stock::QUIT),sigc::mem_fun(*this, &MyWindow::on_FileQuit));

	//View menu:
	m_refActionGroup->add(Gtk::Action::create("ViewMenu", "View"));
    m_refActionGroup->add(Gtk::Action::create("ViewXDirMenu", "X-Dir"));
    Gtk::RadioAction::Group ViewXDirDrawDirChoice;    
    m_refActionGroup->add(Gtk::RadioAction::create(ViewXDirDrawDirChoice, "ViewXDirDrawDirX", "X"),sigc::mem_fun(*this, &MyWindow::on_ViewXdirDrawDirX));
    m_refActionGroup->add(Gtk::RadioAction::create(ViewXDirDrawDirChoice, "ViewXDirDrawDirY", "Y"),sigc::mem_fun(*this, &MyWindow::on_ViewXdirDrawDirY));
    m_refActionGroup->add(Gtk::RadioAction::create(ViewXDirDrawDirChoice, "ViewXDirDrawDirZ", "Z"),sigc::mem_fun(*this, &MyWindow::on_ViewXdirDrawDirZ));
	m_refActionGroup->add(Gtk::Action::create("ViewXDirInvert", "Invert X"),sigc::mem_fun(*this, &MyWindow::on_ViewXDirInvert));
    m_refActionGroup->add(Gtk::Action::create("ViewYDirMenu", "Y-Dir"));
    Gtk::RadioAction::Group ViewYDirDrawDirChoice;
    m_refActionGroup->add(Gtk::RadioAction::create(ViewYDirDrawDirChoice, "ViewYDirDrawDirX", "X"),sigc::mem_fun(*this, &MyWindow::on_ViewYdirDrawDirX));
    m_refActionGroup->add(Gtk::RadioAction::create(ViewYDirDrawDirChoice, "ViewYDirDrawDirY", "Y"),sigc::mem_fun(*this, &MyWindow::on_ViewYdirDrawDirY));
    m_refActionGroup->add(Gtk::RadioAction::create(ViewYDirDrawDirChoice, "ViewYDirDrawDirZ", "Z"),sigc::mem_fun(*this, &MyWindow::on_ViewYdirDrawDirZ));
    m_refActionGroup->add(Gtk::Action::create("ViewColorByMenu", "ColorBy"));
    Gtk::RadioAction::Group ViewColorByChoice;
    
    m_refActionGroup->add(Gtk::RadioAction::create(ViewColorByChoice, "ViewColorByXPos", "X-pos"),sigc::bind<ColorBy>(sigc::mem_fun(*this, &MyWindow::on_ViewColorBy),XPos));
    m_refActionGroup->add(Gtk::RadioAction::create(ViewColorByChoice, "ViewColorByYPos", "Y-pos"),sigc::bind<ColorBy>(sigc::mem_fun(*this, &MyWindow::on_ViewColorBy),YPos));
    m_refActionGroup->add(Gtk::RadioAction::create(ViewColorByChoice, "ViewColorByZPos", "Z-pos"),sigc::bind<ColorBy>(sigc::mem_fun(*this, &MyWindow::on_ViewColorBy),ZPos));
    m_refActionGroup->add(Gtk::RadioAction::create(ViewColorByChoice, "ViewColorByXVel", "X-vel"),sigc::bind<ColorBy>(sigc::mem_fun(*this, &MyWindow::on_ViewColorBy),XVel));
    m_refActionGroup->add(Gtk::RadioAction::create(ViewColorByChoice, "ViewColorByYVel", "Y-vel"),sigc::bind<ColorBy>(sigc::mem_fun(*this, &MyWindow::on_ViewColorBy),YVel));
    m_refActionGroup->add(Gtk::RadioAction::create(ViewColorByChoice, "ViewColorByZVel", "Z-vel"),sigc::bind<ColorBy>(sigc::mem_fun(*this, &MyWindow::on_ViewColorBy),ZVel));
    m_refActionGroup->add(Gtk::RadioAction::create(ViewColorByChoice, "ViewColorByAbsVel", "Velocity"),sigc::bind<ColorBy>(sigc::mem_fun(*this, &MyWindow::on_ViewColorBy),AbsVel));
	m_refActionGroup->add(Gtk::Action::create("ViewYDirInvert", "Invert Y"),sigc::mem_fun(*this, &MyWindow::on_ViewYDirInvert));
    m_refActionGroup->add(Gtk::Action::create("ViewMakeSquared", "Make squared"),sigc::mem_fun(*this, &MyWindow::on_ViewMakeSquared));
	m_refActionGroup->add(Gtk::Action::create("ViewReset", "Reset"),sigc::mem_fun(*this, &MyWindow::on_ViewReset));
    m_refActionGroup->add(Gtk::Action::create("ViewChangeDomain", "Change Domain"),sigc::mem_fun(*this, &MyWindow::on_ViewChangeDomain));
    m_refActionGroup->add(Gtk::Action::create("ViewParticleInfoDialog", "Particle Information"),sigc::mem_fun(*this, &MyWindow::on_ViewParticleInfoDialog));


	m_refUIManager = Gtk::UIManager::create();
	m_refUIManager->insert_action_group(m_refActionGroup);

	add_accel_group(m_refUIManager->get_accel_group());

    //Layout the actions in a menubar and toolbar: 
    Glib::ustring ui_info = 
        "<ui>"
        "  <menubar name='MenuBar'>"
        "    <menu action='FileMenu'>"
        "      <menuitem action='FileOpen'/>"
        "      <menuitem action='FileQuit'/>"
        "    </menu>"
        "    <menu action='ViewMenu'>"
        "      <menu action='ViewXDirMenu'>"
        "        <menuitem action='ViewXDirDrawDirX'/>"
        "        <menuitem action='ViewXDirDrawDirY'/>"
        "        <menuitem action='ViewXDirDrawDirZ'/>"        
        "        <separator/>"
        "        <menuitem action='ViewXDirInvert'/>"
        "      </menu>"
        "      <menu action='ViewYDirMenu'>"
        "        <menuitem action='ViewYDirDrawDirX'/>"
        "        <menuitem action='ViewYDirDrawDirY'/>"
        "        <menuitem action='ViewYDirDrawDirZ'/>"        
        "        <separator/>"        
        "        <menuitem action='ViewYDirInvert'/>"
        "      </menu>"
        "      <menu action='ViewColorByMenu'>"
        "        <menuitem action='ViewColorByXPos'/>"
        "        <menuitem action='ViewColorByYPos'/>"
        "        <menuitem action='ViewColorByZPos'/>"        
        "        <menuitem action='ViewColorByXVel'/>"
        "        <menuitem action='ViewColorByYVel'/>"
        "        <menuitem action='ViewColorByZVel'/>"        
        "        <menuitem action='ViewColorByAbsVel'/>"        
        "      </menu>" 
        "      <menuitem action='ViewMakeSquared'/>"
        "      <menuitem action='ViewReset'/>"
        "      <menuitem action='ViewChangeDomain'/>"
        "      <menuitem action='ViewParticleInfoDialog'/>"        
        "    </menu>"
        "  </menubar>"
        "</ui>";

    try
    {
        m_refUIManager->add_ui_from_string(ui_info);
    }
    catch(const Glib::Error& ex)
    {
        std::cerr << "building menus failed: " <<  ex.what();
    }
  
    /*ParticleInfoTab* test;
    test = new ParticleInfoTab(0,problem);
    particleInfoTabVector.push_back(test);
    particleInfoTabs.append_page(*particleInfoTabVector.back(), "page 1");
    test = new ParticleInfoTab(1,problem);
    particleInfoTabVector.push_back(test);
    particleInfoTabs.append_page(*particleInfoTabVector.back(), "page 2");*/
    
    
    Gtk::Widget* pMenubar = m_refUIManager->getCGWidthidget("/MenuBar");
    m_Box.pack_start(*pMenubar, Gtk::PACK_SHRINK);
    m_Box.pack_start(*area);
    show_all_children();
    
    particleInfoDialog=new ParticleInfoDialog;
}

MyWindow::~MyWindow()
{
}

void MyWindow::on_FileQuit()
{
  delete problem;
  hide(); //Closes the main window to stop the Gtk::Main::run().
}

void MyWindow::on_ViewXdirDrawDirX(){area->setXDrawDir(0);queue_draw();}
void MyWindow::on_ViewXdirDrawDirY(){area->setXDrawDir(1);queue_draw();}
void MyWindow::on_ViewXdirDrawDirZ(){area->setXDrawDir(2);queue_draw();}
void MyWindow::on_ViewYdirDrawDirX(){area->setYDrawDir(0);queue_draw();}
void MyWindow::on_ViewYdirDrawDirY(){std::cout<<"on_ViewYdirDrawDirY"<<std::endl;area->setYDrawDir(1);queue_draw();}
void MyWindow::on_ViewYdirDrawDirZ(){area->setYDrawDir(2);queue_draw();}
void MyWindow::on_ViewXDirInvert(){area->Xflip();queue_draw();}
void MyWindow::on_ViewYDirInvert(){area->Yflip();queue_draw();}
void MyWindow::on_ViewColorBy(ColorBy colorBy)
{
    {
        std::cout<<"on_ViewColorBy"<<colorBy<<std::endl;
        area->setColorBy(colorBy);
        queue_draw();
    }
}
void MyWindow::on_ViewMakeSquared(){area->make_squared();queue_draw();}
void MyWindow::on_ViewReset(){area->resetView();queue_draw();}
void MyWindow::on_ViewChangeDomain()
{
    ChangeDomainDialog dialog(area->getDomainMin(),area->getDomainMax());
    int result = dialog.run();

    switch(result)
    {
        case(Gtk::RESPONSE_OK):
        {
            Vec3D domainMin,domainMax;
            dialog.getDomainMinMax(domainMin,domainMax);
            area->setDomainMin(domainMin);
            area->setDomainMax(domainMax);
            break;
        }
        case(Gtk::RESPONSE_CANCEL):
        {
            break;
        }
        default:
        {
            std::cout<<"Unexpected button clicked. ("<<result<<")"<<std::endl;
            break;
        }
    }        
}
void MyWindow::on_ViewParticleInfoDialog()
{
    particleInfoDialog->show(); 
}
void MyWindow::on_FileOpen()
{
    Gtk::FileChooserDialog dialog("Choose File",Gtk::FILE_CHOOSER_ACTION_OPEN);
    dialog.add_button("Cancel", Gtk::RESPONSE_CANCEL);
    dialog.add_button("Select", Gtk::RESPONSE_OK);
    dialog.set_current_folder("/home/dinant/HgMD/NewLocation/DRIVERS/Graphics");
    int result = dialog.run();
    switch(result)
    {
        case(Gtk::RESPONSE_OK):
        {
            std::cout << "Open clicked." << std::endl;
            std::string filename = dialog.get_filename();
            std::cout << "File selected: " <<  filename << std::endl;
            delete problem;
            problem = new DPMBase;
            problem->set_data_filename(filename);
            problem->open_data_file(std::fstream::in);
            problem->readNextDataFile();
            isParticleSelected.resize(problem->particleHandler.getNumberOfObjects());
            area->resetView();
            break;
        }
        case(Gtk::RESPONSE_CANCEL):
        {
            std::cout << "Cancel clicked." << std::endl;
            break;
        }
        default:
        {
            std::cout << "Unexpected button clicked. ("<<result<<")" << std::endl;
            break;
        }
    }    
}

bool MyWindow::on_key_press_event(GdkEventKey* event)
{
    switch(event->keyval)
    {
        case GDK_KEY_Escape:
        {
            delete problem;
            hide();
            return true;
        }
        case GDK_KEY_n:
        {
            if(problem)
            {
                problem->readNextDataFile();
                isParticleSelected.resize(problem->particleHandler.getNumberOfObjects());
                queue_draw();
            }
            return true;
        }
        default:
        {
            std::cout<<"Unknown key: "<<event->keyval<<std::endl;
            std::cout<<"Current working keys are:"<<std::endl;
            std::cout<<"n           go to next time step"<<std::endl;
            std::cout<<"Escpae      exits program"<<std::endl;
            return true;
        }
    }
}
