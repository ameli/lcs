/*
 * =====================================================================================
 *
 *       Filename:  TestSmoothStructuredPoints.cxx
 *
 *    Description:  Test for Smooth Structured Points
 *
 *        Version:  1.0
 *        Created:  04/27/2014 02:00:35 PM
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

#include <vtkSmartPointer.h>
#include "SmoothStructuredPoints.h"
#include "Deformation.h"
#include <vtkDoubleArray.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPoints.h>
#include <vtkPointData.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkImageGradient.h>


// ====
// Main
// ====

int main(int arhc,char *argv[])
{
    // Reader
    vtkSmartPointer<vtkStructuredPointsReader> Reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
    Reader->SetFileName("../data/input-101-10.vtk");
    Reader->Update();

    // Deformaton
    vtkSmartPointer<Deformation> DeformationFilter = vtkSmartPointer<Deformation>::New();
    DeformationFilter->SetInputConnection(Reader->GetOutputPort());
    DeformationFilter->SmoothingOff();
    DeformationFilter->ProgressOff();
    DeformationFilter->Update();

    // Get Structured Points Data
    vtkSmartPointer<vtkStructuredPoints> InputImageData = vtkStructuredPoints::SafeDownCast(DeformationFilter->GetOutput());

    // Get Input Scalar Field
    vtkSmartPointer<vtkDoubleArray> InputScalarField = 
        vtkDoubleArray::SafeDownCast(InputImageData->GetPointData()->GetArray("MaxStrainValues"));

    // Get Input Vector Field
    vtkSmartPointer<vtkDoubleArray> InputVectorField =
        vtkDoubleArray::SafeDownCast(InputImageData->GetPointData()->GetArray("MaxStrainVectors"));

    // Settings
    int *GridResolution = InputImageData->GetDimensions();
    unsigned int KernelSize = 7;

    // 1- Smoothing Input Scalar Field
    vtkSmartPointer<vtkDoubleArray> OutputScalarField = vtkSmartPointer<vtkDoubleArray>::New();
    SmoothStructuredPoints::SmoothScalarField(
            KernelSize,
            GridResolution,
            InputScalarField,
            OutputScalarField);

    // 2- Smoothing Input Scalar Field using VTK
    vtkSmartPointer<vtkDoubleArray> OutputScalarFieldUseVTK = vtkSmartPointer<vtkDoubleArray>::New();
    SmoothStructuredPoints::SmoothScalarFieldUseVTK(
            KernelSize,
            GridResolution,
            InputScalarField,
            OutputScalarFieldUseVTK);

    // 3- Smooth Unit Vector Field
    vtkSmartPointer<vtkDoubleArray> OutputVectorField = vtkSmartPointer<vtkDoubleArray>::New();
    SmoothStructuredPoints::SmoothUnitVectorField(
            KernelSize,
            GridResolution,
            InputVectorField,
            OutputVectorField);

    // 4- Smooth Gradient Field
    vtkSmartPointer<vtkDoubleArray> OutputGradientField = vtkSmartPointer<vtkDoubleArray>::New();
    SmoothStructuredPoints::SmoothGradientOfScalarField(
            KernelSize,
            GridResolution,
            InputScalarField,
            OutputGradientField);

    // 5- Non-Smooth Gradient Field
    vtkSmartPointer<vtkDataArray> NonSmoothedGradientField;
    vtkSmartPointer<vtkStructuredPoints> InputStructuredPoints = vtkSmartPointer<vtkStructuredPoints>::New();
    InputStructuredPoints->SetOrigin(0,0,0);
    InputStructuredPoints->SetSpacing(1,1,1);
    InputStructuredPoints->SetDimensions(GridResolution);
    InputStructuredPoints->GetPointData()->SetScalars(InputScalarField);
    vtkSmartPointer<vtkImageGradient> GradientFilter = vtkSmartPointer<vtkImageGradient>::New();
    #if VTK_MAJOR_VERSION <= 5
        GradentFilter->SetInput(InputStructuredPoints);
    #else
        GradientFilter->SetInputData(InputStructuredPoints);
    #endif
    GradientFilter->SetDimensionality(3);
    GradientFilter->HandleBoundariesOn();
    GradientFilter->Update();
    NonSmoothedGradientField = GradientFilter->GetOutput()->GetPointData()->GetScalars();
    NonSmoothedGradientField->SetName("NonSMoothedGradientField");

    // 6- Smooth Hessian
    double **HessianAtPoint = new double*[3];
    for(int i = 0; i< 3; i++)
    {
        HessianAtPoint[i] = new double[3];
    }

    int PointId = 0;  // a test point
    SmoothStructuredPoints::SmoothHessianAtPoint(
            KernelSize,
            GridResolution,
            PointId,
            InputScalarField,
            HessianAtPoint);

    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            std::cout << HessianAtPoint[i][j] << ",\t";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    for(int i = 0; i < 3; i++)
    {
        delete [] HessianAtPoint[i];
    }
    delete [] HessianAtPoint;

    // Create Output ImageData
    vtkSmartPointer<vtkStructuredPoints> OutputImageData = vtkSmartPointer<vtkStructuredPoints>::New();
    OutputImageData->CopyStructure(InputImageData);
    OutputImageData->GetPointData()->AddArray(InputScalarField);
    OutputImageData->GetPointData()->AddArray(OutputScalarField);
    OutputImageData->GetPointData()->AddArray(OutputScalarFieldUseVTK);
    OutputImageData->GetPointData()->AddArray(InputVectorField);
    OutputImageData->GetPointData()->AddArray(OutputVectorField);
    OutputImageData->GetPointData()->AddArray(OutputGradientField);
    OutputImageData->GetPointData()->AddArray(NonSmoothedGradientField);

    // Write Output
    vtkSmartPointer<vtkStructuredPointsWriter> Writer = vtkSmartPointer<vtkStructuredPointsWriter>::New();
    Writer->SetInputData(OutputImageData);
    Writer->SetFileName("../data/OutputImageData.vtk");
    Writer->Update();

    std::cout << "Output has been written to ../data/OutputImageData.vtk" << std::endl;

    return EXIT_SUCCESS;
}
