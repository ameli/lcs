/*
 * =====================================================================================
 *
 *       Filename:  TestIdentifyStructures.cxx
 *
 *    Description:  Test for IdentifyStructures class
 *
 *        Version:  1.0
 *        Created:  04/10/2014 04:57:13 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University Of California, Berkeley
 *
 * =====================================================================================
 */

// =======
// Headers
// =======

#include "IdentifyStructures.h"
#include <vtkSmartPointer.h>
#include <vtkStructuredPointsReader.h>
#include <vtkPolyDataWriter.h>

#ifndef HERE
#define HERE std::cout << "DEBUG: " << __FILE__ << ",  line  " << __LINE__ << std::endl;
#endif

// ====
// Main
// ====

int main(int argc,char *argv[])
{
    // Reader
    vtkSmartPointer<vtkStructuredPointsReader> Reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
    Reader->SetFileName("../data/input-sample.vtk");
    Reader->Update();

    // Identify Structures
    vtkSmartPointer<IdentifyStructures> Filter = vtkSmartPointer<IdentifyStructures>::New();
    Filter->SetInputConnection(Reader->GetOutputPort());
    Filter->Update();

    // Writer
    vtkSmartPointer<vtkPolyDataWriter> Writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    Writer->SetInputConnection(Filter->GetOutputPort());
    Writer->SetFileName("../data/output-sample.vtk");
    Writer->Update();
    
    return EXIT_SUCCESS;
}
