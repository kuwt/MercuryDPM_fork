//
// Created by bert on 5/15/17.
//
#include <iostream>
#include <fstream>
void WriteVTK(int numTimeSteps, std::string name);


int main()
{
    int Tsteps=841;
    WriteVTK(Tsteps,"TriolietDietFeederProcessor_0_Particle");
    WriteVTK(Tsteps,"TriolietDietFeederProcessor_1_Particle");
    WriteVTK(Tsteps,"TriolietDietFeederProcessor_2_Particle");
    WriteVTK(Tsteps,"TriolietDietFeederProcessor_3_Particle");
}

void WriteVTK(int numTimeSteps, std::string name)
{
    std::ofstream VTUfile;
    //write headers
    VTUfile.open (name+".pvd");
    VTUfile << "<VTKFile type=\"Collection\">\n"
            "   <Collection>\n";
    //write all filenames with time steps upto numTimeSteps
    for ( int i = 0; i < numTimeSteps; i++ ) {
        VTUfile <<  "      <DataSet timestep=\""<<i<<"\" file=\""<<name<<"_"<<i<<".vtu\"/>\n";
    }

    VTUfile <<"   </Collection> \n</VTKFile>\n";

    //close file
    VTUfile.close();
}
