/*
#include "UnderResolvedMesh.h"

// TODO make all these doc_funcs part of UnderResolvedCoupling as it needs mesh_pt() to point at this class, add Mesh::output_voidage_byEl from mesh.cc from desktop here

//==============start_doc===========================================
/// Doc the solution
//==================================================================
template <class ELEMENT>
void UnderResolvedMesh<ELEMENT>::doc_solution(oomph::DocInfo& doc_info)
{
*/
/*    std::ofstream some_file;
    char filename[100];

    // Number of plot points
    unsigned npts;
    npts=5;

    // Output solution
    sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),doc_info.number());
    some_file.open(filename);
    mesh_pt()->output(some_file,npts);
    some_file.close();*//*

} //end doc_solution


//==============start_doc_void======================================
/// Doc the voidage
//==================================================================
template <class ELEMENT>
void UnderResolvedMesh<ELEMENT>::doc_voidage(oomph::DocInfo& doc_info)
{
*/
/*    //std::cout << "In doc_voidage" << std::endl;
    std::ofstream some_file;
    char filename[100];

    // number of plot points
    unsigned npts = 1;

    // Output solution if using get_voidage_byEl
    sprintf(filename,"%s/voidagen%i.dat",doc_info.directory().c_str(),doc_info.number());
    some_file.open(filename);
    mesh_pt()->output_voidage_byEl(some_file,npts);
    some_file.close();*//*

} //end doc_voidage


//==============start_doc_element======================================
/// Doc the element data
//==================================================================
template <class ELEMENT>
void UnderResolvedMesh<ELEMENT>::doc_element(oomph::DocInfo& doc_info)
{
*/
/*    std::ofstream some_file;
    char filename[100];

    // number of plot points
    unsigned npts = 1;

    // Output solution
    sprintf(filename,"%s/elements%i.dat",doc_info.directory().c_str(),doc_info.number());
    some_file.open(filename);
    mesh_pt()->output(some_file,npts);
    some_file.close();*//*

} //end doc_voidage
*/
