/*
 * =====================================================================================
 *
 *       Filename:  main.cxx
 *
 *    Description:  Lagrangian Coherence Structures
 *
 *        Version:  1.0
 *        Created:  02/12/2014 05:05:20 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University Of California, Berkeley
 *
 * =====================================================================================
 */

// ===========
// Definitions
// ===========

// Disable MSDN security features (WINDOWS only)
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

// =======
// Headers
// =======

#include <vtkSmartPointer.h>
#include <vtkStructuredPointsReader.h>
#include "Deformation.h"
#include "IdentifyStructures.h"
#include <vtkPolyDataWriter.h>
#include <string.h>

// ====
// Main
// ====

int main(int argc, char *argv[])
{
    // Inuput filenames
    char InputFilename[256];
    if(argc > 1)
    {
        strcpy(InputFilename,argv[1]);
    }
    else
    {
        strcpy(InputFilename,"../data/input-101.vtk");
    }

    // Output filename
    char OutputFilename[256];
    if(argc > 2)
    {
        strcpy(OutputFilename,argv[2]);
    }
    else
    {
        strcpy(OutputFilename,"../data/output-101.vtk");
    }

    // Read Data
    vtkSmartPointer<vtkStructuredPointsReader> Reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
    Reader->SetFileName(InputFilename);
    Reader->Update();

    // Deformation
    vtkSmartPointer<Deformation> DeformationFilter = vtkSmartPointer<Deformation>::New();
    DeformationFilter->SetInputConnection(Reader->GetOutputPort());
    DeformationFilter->SmoothingOn();
    DeformationFilter->SmoothingOff();
    // DeformationFilter->SetSmoothingKernelSize(31);
    DeformationFilter->ProgressOn();
    DeformationFilter->Update();

    // Identify Structures
    vtkSmartPointer<IdentifyStructures> IdentifyStructuresFilter = vtkSmartPointer<IdentifyStructures>::New();
    IdentifyStructuresFilter->SetInputConnection(DeformationFilter->GetOutputPort());
    IdentifyStructuresFilter->SetStructureToMaxStrain();
    // IdentifyStructuresFilter->SetStructureToMinStrain();
    // IdentifyStructuresFilter->SmoothingOn();
    IdentifyStructuresFilter->SmoothingOff();
    // IdentifyStructuresFilter->SetSmoothingKernelSize(13);
    IdentifyStructuresFilter->DebugOff();
    IdentifyStructuresFilter->ProgressOn();
    IdentifyStructuresFilter->Update();

    // Write Data
    vtkSmartPointer<vtkPolyDataWriter> Writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    Writer->SetInputConnection(IdentifyStructuresFilter->GetOutputPort());
    Writer->SetFileName(OutputFilename);
    Writer->Update();

    return EXIT_SUCCESS;
}
