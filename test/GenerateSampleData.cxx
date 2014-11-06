/*
 * =====================================================================================
 *
 *       Filename:  GenerateSampleData.cxx
 *
 *    Description:  Generate sample data for testing Marching cubes (without PCA)
 *
 *        Version:  1.0
 *        Created:  04/10/2014 09:39:59 AM
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
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <cstdlib>   // for std::div
#include <cmath>  // for GeneratePointData

// ======
// Macros
// ======

#define DIMENSION 3
// #define HERE std::cout << "DEBUG: " << __FILE__ << ",  line  " << __LINE__ << std::endl;

// ==========
// Prototypes
// ==========

void Id2Index(int *GridResolution,int Id,int *Index);
int Index2Id(int *GridResolution,int Index);
void GeneratePointData(int *GridResolution,int *Index,double & Scalar,double Vector[DIMENSION]);

// ====
// Main
// ====

int main(int argc, char *argv[])
{
    // Grid Info
    int GridResolution[DIMENSION] = {31,31,31};
    double Origin[DIMENSION] = {0.0,0.0,0.0};
    double Corner[DIMENSION] = {1.0,1.0,1.0};
    double Spacing[DIMENSION];
    for(int i = 0; i < DIMENSION; i++)
    {
        Spacing[i] = (Corner[i] - Origin[i])/(GridResolution[i]-1);;
    }

    // Generate Grid
    vtkSmartPointer<vtkStructuredPoints> StructuredPoints = vtkSmartPointer<vtkStructuredPoints>::New();
    StructuredPoints->SetDimensions(GridResolution);
    StructuredPoints->SetOrigin(Origin);
    StructuredPoints->SetSpacing(Spacing);

    // Generate Scalars
    int NumPoints = StructuredPoints->GetNumberOfPoints();
    vtkSmartPointer<vtkDoubleArray> Scalars = vtkSmartPointer<vtkDoubleArray>::New();
    Scalars->SetNumberOfTuples(NumPoints);
    Scalars->SetNumberOfComponents(1);
    Scalars->SetName("MaxStrainValues");

    vtkSmartPointer<vtkDoubleArray> Vectors = vtkSmartPointer<vtkDoubleArray>::New();
    Vectors->SetNumberOfComponents(DIMENSION);
    Vectors->SetNumberOfTuples(NumPoints);
    Vectors->SetName("MaxStrainVectors");

    for(int Id = 0; Id < NumPoints; Id++)
    {
        int Index[DIMENSION];
        Id2Index(GridResolution,Id,Index);

        // std::cout << "Id: " << Id << ", \tIndex: " << Index[0] << ", " << Index[1] << ", " << Index[2] << std::endl;

        // Generate Scalars
        double Scalar;
        double Vector[DIMENSION];
        GeneratePointData(GridResolution,Index,Scalar,Vector);

        Scalars->SetTuple1(Id,Scalar);
        Vectors->SetTuple(Id,Vector);
        // Vectors->SetTuple3(Id,Vector[0],Vector[1],Vector[2]);
    }

    // Add Salars to Grid
    StructuredPoints->GetPointData()->AddArray(Scalars);
    StructuredPoints->GetPointData()->AddArray(Vectors);

    // Write
    vtkSmartPointer<vtkStructuredPointsWriter> Writer = vtkSmartPointer<vtkStructuredPointsWriter>::New();
    Writer->SetInputData(StructuredPoints);
    Writer->SetFileName("../data/input-sample.vtk");
    Writer->Update();

    return EXIT_SUCCESS;
}

// ==========
// Id 2 Index
// ==========

void Id2Index(int *GridResolution,int Id,int *Index)
{
    for(int i = 0; i < DIMENSION; i++)
    {
        div_t div;
        div = std::div(Id,GridResolution[i]);
        Index[i] = div.rem;
        Id = div.quot;
    }
}

// ==========
// Index 2 Id
// ==========

int Index2Id(int *GridResolution,int *Index)
{
    int Id = 0;
    int Product = 1;

    for(int i = 0; i < DIMENSION; i++)
    {
        Id += Index[i] * Product;
        Product *= GridResolution[i];
    }

    return Id;
}

// ================
// Generate Scalars
// ================

void GeneratePointData(int *GridResolution,int *Index,double & Scalar,double Vector[DIMENSION])
{
    // Point coordinates
    double x = double(Index[0]) / GridResolution[0];
    double y = double(Index[1]) / GridResolution[1];
    double z = double(Index[2]) / GridResolution[2];

    // Output
    // Scalar = x+y+z-0.5;
    // Scalar = sin(4*x)*sin(1.0/(x+1))-cos(2*z)*sin(3*y);
    Scalar = sin(5*x)*x + sin(5*y)*y + z*z-1;
    Vector[0] = x;
    Vector[1] = y;
    Vector[2] = z;
}
