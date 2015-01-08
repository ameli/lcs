/*
 * =====================================================================================
 *
 *       Filename:  SmoothStructuredPoints.cxx
 *
 *    Description:  Smooth Structured Points
 *
 *        Version:  1.0
 *        Created:  04/17/2014 10:32:21 AM
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

#include "SmoothStructuredPoints.h"

// Pipeline
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDataObject.h>
#include <vtkObjectFactory.h>
#include <vtkDemandDrivenPipeline.h>

// General
#include <vtkSmartPointer.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkMath.h>

// Algorithms
#include <vtkImageGaussianSmooth.h>

// Arrays
#include <vtkDoubleArray.h>

// Data
#include <vtkStructuredData.h>
#include <vtkImageData.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>

// STL
#include <cstdlib>  // for exit
#include <cmath>    // for cos, sin, sqrt, floor

// ======
// Macros
// ======

vtkStandardNewMacro(SmoothStructuredPoints);
vtkCxxRevisionMacro(SmoothStructuredPoints,"$Revision 1.0$");
#define DIMENSION 3

// ===========
// Constructor
// ===========

SmoothStructuredPoints::SmoothStructuredPoints()
{
}

// ==========
// Destructor
// ==========

SmoothStructuredPoints::~SmoothStructuredPoints()
{
}

// ==========
// Print Self
// ==========

void SmoothStructuredPoints::PrintSelf(ostream &os,vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}

// ============
// Request Data
// ============

int SmoothStructuredPoints::RequestData(
        vtkInformation *vtkNotUsed(Request),
        vtkInformationVector **InputVector,
        vtkInformationVector *OutputVector)
{

    return 1;
}

// ===============
// Get Point Index
// ===============

// Decsription:
// Converts the list style ID if a point into coordinate-wise indices if the point
// in a structured grid. Grid Resolution is number of points in each direction.
// Note: PointId is 0-offset. So the first point is 0. Therefore, the PointIndex are
// also offset from 0. So the first point is (0,0,0) in 3D.

void StructuredPointsHelper::GetPointIndex(
        int *GridResolution,
        int PointId,
        int *PointIndex)  // Output
{
    for(unsigned int DimensionIterator = 0; DimensionIterator < DIMENSION; DimensionIterator++)
    {
        div_t DivisionResult;
        DivisionResult = std::div(static_cast<int>(PointId),GridResolution[DimensionIterator]);
        PointIndex[DimensionIterator] = DivisionResult.rem;
        PointId = DivisionResult.quot;
    }
}

// ============
// Get Point Id
// ============

// Description:
// Given coordinate index for a point, it returns the point Id in structured grid.
// All components of PointIndex are offset form 0. So the output PointId will be 
// offset from 0. GridResolution is number of points at each direction.

inline int StructuredPointsHelper::GetPointId(
        int *GridResolution,
        int *PointIndex)
{
    int PointId = 0;
    unsigned int DimensionProduct = 1;

    for(unsigned int DimensionIterator = 0; DimensionIterator < DIMENSION; DimensionIterator++)
    {
        PointId += DimensionProduct * PointIndex[DimensionIterator];
        DimensionProduct *= GridResolution[DimensionIterator];
    }

    return PointId;
}

// ===============
// Get Stencil Ids
// ===============

// Description:
// Finds the point Ids of stencil neighbors  of a point on a structured grid.
// Note: If point is on the boundary, it considers the point itself instead of a stencil outdside
// the domain. Thus the central point duplicates.
// Note: The stencil is only two points in a specific direction. Direction = 0 returns the two 
// stencilpoints in x-axis, Direction = 1 returns the two stencils in y-axis, etc.
// Output is an array with two elements as following
// StencilIds[0] is the back point Id.
// StencilIds[1] is the front point Id.

void StructuredPointsHelper::GetStencilIds(
        int *GridResolution,
        unsigned int PointId,
        unsigned int TargetDirection,
        unsigned int *StencilIds)   // Output
{
    // Get point index
    int PointIndex[DIMENSION];
    StructuredPointsHelper::GetPointIndex(GridResolution,PointId,PointIndex);

    // Stencils Index
    int BackStencilIndex[DIMENSION];
    int FrontStencilIndex[DIMENSION];

    for(unsigned int DimensionIterator = 0; DimensionIterator < DIMENSION; DimensionIterator++)
    {
        // Initialize Stencil indices the same as point index
        BackStencilIndex[DimensionIterator] = PointIndex[DimensionIterator];
        FrontStencilIndex[DimensionIterator] = PointIndex[DimensionIterator];

        // Set stencil indices in the requested direction
        if(DimensionIterator == TargetDirection)
        {
            // Check right boundaries for front stencil
            if(PointIndex[TargetDirection] < GridResolution[TargetDirection]-1)
            {
                FrontStencilIndex[TargetDirection]++;
            }

            // Check left boundaries for back stencil
            if(PointIndex[TargetDirection] > 0)
            {
                BackStencilIndex[TargetDirection]--;
            }
        }
    }

    // Convert Stencil index to list Id
    StencilIds[0] = vtkStructuredData::ComputePointId(GridResolution,BackStencilIndex);  // Back stencil
    StencilIds[1] = vtkStructuredData::ComputePointId(GridResolution,FrontStencilIndex); // Front stencil
}

// ====================
// Clip Index To Bounds
// ====================

// Desciption:
// Indices are 0-offset. So the left bound is 0 and the most right index is RightBound-1.
// The index is clipped between 0 and RightBound-1. This function produces Newmann boundary
// condition.

void StructuredPointsHelper::ClipIndexToBounds(
        int RightBound,
        int &Index)   // Input / Output
{
    int LeftBound = 0;
    if(Index < LeftBound)
    {
        Index = 0;
    }
    else if(Index > RightBound -1)
    {
        Index = RightBound - 1;
    }
}

// ====================
// Tile Index To Bounds
// ====================

// Description:
// This will produce a peridoc boundary condition.

void StructuredPointsHelper::TileIndexToBounds(
        int RightBound,
        int &Index)
{
    int LeftBound = 0;
    if(Index < LeftBound)
    {
        Index = RightBound + Index;
    }
    else if(Index > RightBound - 1)
    {
        Index = Index - RightBound;
    }
}

// ==============
// Gauss Function
// ==============

// Description:
// By setting dimension, the integral under Guass function in that dimension will be 1.
// For separable kernels, we use dimension = 1;

double SmoothingHelper::GaussFunction(
        double Distance,
        double StandardDeviation,
        unsigned int Dimension)
{
    double Coefficient = 1.0 / pow(StandardDeviation * sqrt(2.0 * vtkMath::Pi()),int(Dimension));
    double Exponent = exp(-static_cast<double>(Distance * Distance) / 
            (2.0 * StandardDeviation * StandardDeviation));

    return Coefficient * Exponent;
}

// ===============================
// Gauss Function First Derivative
// ===============================

double SmoothingHelper::GaussFunctionFirstDerivative(
        double Distance,
        double StandardDeviation,
        unsigned int Dimension)
{
    double FirstDerivativeCoefficient = -Distance / (StandardDeviation * StandardDeviation);
    return FirstDerivativeCoefficient * 
        SmoothingHelper::GaussFunction(Distance,StandardDeviation,Dimension);
}

// ================================
// Gauss Function Second Derivative
// ================================

double SmoothingHelper::GaussFunctionSecondDerivative(
        double Distance,
        double StandardDeviation,
        unsigned int Dimension)
{
    double SecondDerivativeCoefficient = (Distance*Distance - StandardDeviation * StandardDeviation) / 
        pow(StandardDeviation,4);
    return SecondDerivativeCoefficient * 
        SmoothingHelper::GaussFunction(Distance,StandardDeviation,Dimension);
}

// ===========================
// Initialize Smoothing Kernel
// ===========================

// Description:
// This method assumes that:
// 1- Kernel is to used in 1-D convolution, such as in Separated kernels.
// 2- Standard Deviation is set such that the Kernel size spans the 99% confidence bound.
// To do so, Standard Deviation is set to  one-third of Kernel radius.

void SmoothingHelper::InitializeSmoothingKernel(
        unsigned int KernelSize,
        double *SmoothingKernel)
{
}

// ====================================
// Initialize Smoothing Gradient Kernel
// ====================================

// Description:
// This method assumes that:
// 1- Kernel is to used in 1-D convolution, such as in Separated kernels.
// 2- Standard Deviation is set such that the Kernel size spans the 99% confidence bound.
// To do so, Standard Deviation is set to  one-third of Kernel radius.

void SmoothingHelper::InitializeSmoothingGradientKernel(
        unsigned int KernelSize,
        double *SmoothingGradientKernel)
{
    // Check Kernel size is an odd number
    if(KernelSize < 3)
    {
        std::cerr << "KernelSize should be greater or equal to 3." << std::endl;
        exit(0);
    }
    else if((floor(static_cast<double>(KernelSize-1)) / 2.0) != floor(static_cast<double>(KernelSize) / 2.0))
    {
        std::cerr << "Kernel size should be an odd number." << std::endl;
        exit(0);
    }

    int KernelCenter = floor(static_cast<double>(KernelSize)/2.0);
    double StandardDeviation = static_cast<double>(KernelCenter) / 3.0; // 99% bound
    unsigned int Dimension = 1; // Separated kernel

    // Initialize Kernel
    for(int Iterator = 0; Iterator < static_cast<int>(KernelSize); Iterator++)
    {
        double Distance = static_cast<double>(Iterator - KernelCenter);
        SmoothingGradientKernel[Iterator] = SmoothingHelper::GaussFunctionFirstDerivative(
                Distance,StandardDeviation,Dimension);
    }
}

// ===================================
// Initialize Smoothing Hessian Kernel
// ===================================

// Description:
// This method assumes that:
// 1- Kernel is to used in 1-D convolution, such as in Separated kernels.
// 2- Standard Deviation is set such that the Kernel size spans the 99% confidence bound.
// To do so, Standard Deviation is set to  one-third of Kernel radius.

void SmoothingHelper::InitializeSmoothingHessianKernel(
        unsigned int KernelSize,
        double *SmoothingHessianKernel)
{
    // Check Kernel size is an odd number
    if(KernelSize < 3)
    {
        std::cerr << "KernelSize should be greater or equal to 3." << std::endl;
        exit(0);
    }
    else if((floor(static_cast<double>(KernelSize-1)) / 2.0) != floor(static_cast<double>(KernelSize) / 2.0))
    {
        std::cerr << "Kernel size should be an odd number." << std::endl;
        exit(0);
    }

    int KernelCenter = floor(static_cast<double>(KernelSize) / 2.0);
    double StandardDeviation = static_cast<double>(KernelCenter) / 3.0; // 99% bound
    unsigned int Dimension = 1; // Separated kernel

    // Initialize Kernel
    for(int Iterator = 0; Iterator < static_cast<int>(KernelSize); Iterator++)
    {
        double Distance = static_cast<double>(Iterator - KernelCenter);
        SmoothingHessianKernel[Iterator] = SmoothingHelper::GaussFunctionSecondDerivative(
                Distance,StandardDeviation,Dimension);
    }
}

// =================================
// Convolute Kernel In One Direction
// =================================

// Description:
// This method will set the number of components and tuples of the output. It is not
// necessary to set these values apriori.

void SmoothingHelper::ConvoluteKernelInOneDirection(
        double *Kernel,
        unsigned int KernelSize,
        unsigned int Direction,
        int *GridResolution,
        vtkDoubleArray *OriginalScalarField,
        vtkDoubleArray *ConvolutedScalarField)
{
    unsigned int NumberOfPoints = OriginalScalarField->GetNumberOfTuples();
    int KernelRadius = static_cast<int>((KernelSize-1) / 2);

    // Prepare Output Scalar Field
    ConvolutedScalarField->SetNumberOfComponents(1);
    ConvolutedScalarField->SetNumberOfTuples(NumberOfPoints);

    // Traverse over points
    for(unsigned int CenterPointId = 0; CenterPointId < NumberOfPoints; CenterPointId++)
    {
        // Get Point Index which is at the center of stencil
        int CenterPointIndex[DIMENSION];
        StructuredPointsHelper::GetPointIndex(GridResolution,CenterPointId,CenterPointIndex);
       
        // Convolute data with Kernel in the specified direction
        double Sum = 0;
        for(int ConvolutionIterator = -KernelRadius;
            ConvolutionIterator <= KernelRadius;
            ConvolutionIterator++)
        {
            unsigned int KernelIndex = ConvolutionIterator + KernelRadius;
            ConvolutionPointIndex[Direction] = static_cast<int>(CenterPointIndex[Direction] - ConvolutionIterator);

            // Handling Boundary
            // StructuredPointsHelper::TileIndexToBounds(GridResolution[Direction],ConvolutionPointIndex[Direction]);
            StructuredPointsHelper::ClipIndexToBounds(GridResolution[Direction],ConvolutionPointIndex[Direction]);

            // Get Convolution Point Id
            int ConvolutionPointId = StructuredPointsHelper::GetPointId(GridResolution,ConvolutionPointIndex);

            #ifdef _OPENMP  // Thread safe region

            // Get Convolution point Scalar data
            double ScalarData[1];
            OriginalScalarField->GetTuple(ConvolutionPointId,ScalarData);

            // Convolute
            Sum += Kernel[KernelIndex] * ScalarData[0];

            #else  // not thread safe region

            // Get Convolution point Scalar data
            double ScalarData = OriginalScalarField->GetTuple1(ConvolutionPointId);

            // Convolute
            Sum += Kernel[KernelIndex] * ScalarData;
            #endif

        }

        // Store the convolution sum in new array
        ConvolutedScalarField->SetTuple1(CenterPointId,Sum);
    }
}

// ===================
// Smooth Scalar Field
// ===================

void SmoothStructuredPoints::SmoothScalarField(
        unsigned int KernelSize,
        int *GridResolution,
        vtkDoubleArray *InputScalarField,
        vtkDoubleArray *OutputScalarField)
{
    // Check Input
    if(InputScalarField == NULL)
    {
        std::cerr << "InputScalarField is NULL." << std::endl;
        exit(0);
    }
    else if(InputScalarField->GetNumberOfTuples() < 1)
    {
        std::cerr << "InputScalarField has no tuples." << std::endl;
        exit(0);
    }
    else if(InputScalarField->GetNumberOfComponents() != 1)
    {
        std::cerr << "InputScalarField is not a scalar array." << std::endl;
        exit(0);
    }

    // Initialize Smoothing Kernel
    double *SmoothingKernel = new double[KernelSize];
    SmoothingHelper::InitializeSmoothingKernel(
            KernelSize,
            SmoothingKernel);   // Output

    // Declare intermediate fields for convolution steps
    vtkDoubleArray *ScalarField[DIMENSION+1];
    ScalarField[0] = InputScalarField;
    ScalarField[DIMENSION] = OutputScalarField;

    // Intermediate Scalar Fields
    for(unsigned int FieldIterator = 1; FieldIterator < DIMENSION; FieldIterator++)
    {
        ScalarField[FieldIterator] = vtkDoubleArray::New();
    }

    // Smooth in each dimension direction
    for(unsigned int DimensionIterator = 0; DimensionIterator < DIMENSION; DimensionIterator++)
    {
        // Assign Input and output Scalar fields for convolution
        vtkDoubleArray *OriginalScalarField = ScalarField[DimensionIterator];
        vtkDoubleArray *ConvolutedScalarField = ScalarField[DimensionIterator+1];

        // Convolute Kernel
        SmoothingHelper::ConvoluteKernelInOneDirection(
                SmoothingKernel,
                KernelSize,
                DimensionIterator,
                GridResolution,
                OriginalScalarField,
                ConvolutedScalarField);
    }

    // Name output
    OutputScalarField->SetName("SmoothedScalarField");

    // Free memory
    delete [] SmoothingKernel;
    for(int FieldIterator = 1; FieldIterator < DIMENSION; FieldIterator++)
    {
        ScalarField[FieldIterator]->Delete();
    }
}

// ========================
// Self Smooth Scalar Field
// ========================

void SmoothStructuredPoints::SelfSmoothScalarField(
        unsigned int KernelSize,
        int *GridResolution,
        vtkDoubleArray *ScalarField)   // Input / Output
{
    SmoothStructuredPoints::SmoothScalarField(
            KernelSize,
            GridResolution,
            ScalarField,
            ScalarField);
}

// ===========================
// Smooth Scalar Field Use VTK
// ===========================

// Description:
// This method uses vtkImageGaussianSmooth. The purpose of this method is to
// compare SmoothScalarField mathid with it. See /test/TestSmoothStructuredPoints.cxx

void SmoothStructuredPoints::SmoothScalarFieldUseVTK(
        unsigned int KernelSize,
        int *GridResolution,
        vtkDoubleArray *InputScalarField,
        vtkDoubleArray *OutputScalarField)
{
    // Create Image Data
    vtkSmartPointer<vtkImageData> ImageData = vtkSmartPointer<vtkImageData>::New();
    ImageData->SetOrigin(0,0,0);
    ImageData->SetSpacing(1,1,1);
    ImageData->SetDimensions(GridResolution);
    ImageData->GetPointData()->SetScalars(InputScalarField);

    // Filter Settings
    // Note: In this code we use KernelRadius, but in the vtk filter they use sigma*radius.
    // So we devide radius by sigma to make the both radia consistent for fair comparison.
    int KernelRadius = floor(static_cast<double>(KernelSize) / 2.0);
    double StandardDeviation = static_cast<double>(KernelRadius) / 3.0;
    // double RadiusFactor = static_cast<double>(KernelRadius) / StandardDeviation;
    double RadiusFactor = 3.01;

    // vtk's Gaussian Filter
    vtkSmartPointer<vtkImageGaussianSmooth> SmoothingFilter = vtkSmartPointer<vtkImageGaussianSmooth>::New();
    #if VTK_MAJOR_VERSION <= 5
        SmoothingFilter->SetInput(ImageData);
    #else
        SmoothingFilter->SetInputData(ImageData);
    #endif
    SmoothingFilter->SetRadiusFactor(RadiusFactor);
    SmoothingFilter->SetStandardDeviation(StandardDeviation);
    SmoothingFilter->SetDimensionality(DIMENSION);
    SmoothingFilter->Update();

    // Output
    OutputScalarField->DeepCopy(SmoothingFilter->GetOutput()->GetPointData()->GetScalars());
    OutputScalarField->SetName("SmoothedScalarFieldUseVTK");
}

// ========================
// Smooth Unit Vector Field
// ========================

// Description:
// This method converts the carteesian representartion of unit vectors into the 
// spherical coordinates (Theta, Phi). Next it smoothes the angle fields. It has the
// advantages of smoothing 2 scalar fields instead of 3. Moreover it preserves the unity
// of smoothed vectors.

void SmoothStructuredPoints::SmoothUnitVectorField(
        unsigned int KernelSize,
        int *GridResolution,
        vtkDoubleArray *InputVectorField,
        vtkDoubleArray *OutputVectorField)   // Output
{
    // Check Input
    if(InputVectorField == NULL)
    {
        std::cerr << "InputVectorField is NULL." << std::endl;
        exit(0);
    }
    else if(InputVectorField->GetNumberOfTuples() < 1)
    {
        std::cerr << "InputVectorField has no tuples." << std::endl;
        exit(0);
    }

    // Set output
    unsigned int NumberOfTuples = InputVectorField->GetNumberOfTuples();
    OutputVectorField->SetNumberOfComponents(DIMENSION);
    OutputVectorField->SetNumberOfTuples(NumberOfTuples);
    OutputVectorField->SetName("SmoothedVectorField");
    

    // Declare Spherical Coordinates Angle arrays
    vtkSmartPointer<vtkDoubleArray> PolarAngleField = vtkSmartPointer<vtkDoubleArray>::New();
    PolarAngleField->SetNumberOfComponents(1);
    PolarAngleField->SetNumberOfTuples(NumberOfTuples);

    vtkSmartPointer<vtkDoubleArray> AzimuthalAngleField = vtkSmartPointer<vtkDoubleArray>::New();
    AzimuthalAngleField->SetNumberOfComponents(1);
    AzimuthalAngleField->SetNumberOfTuples(NumberOfTuples);

    // Convert Cartesian vectors into spherical representation
    for(unsigned int TupleIterator = 0; TupleIterator < NumberOfTuples; TupleIterator++)
    {
        // Get Vector in cartesian coordinates
        #ifdef _OPENMP  // Thread safe region
        double CartesianVector[DIMENSION];
        InputVectorField->GetTuple(TupleIterator,CartesianVector);

        #else  // Not thread safe region
        double *CartesianVector = InputVectorField->GetTuple(TupleIterator);
        #endif

        // Convert Cartesian to Spherical Coordinates
        double SphericalVector[2];
        SmoothingHelper::ConvertCartesianToSphericalCoordinates(CartesianVector,SphericalVector);

        // Store Spherical Fields
        PolarAngleField->SetTuple1(TupleIterator,SphericalVector[0]);
        AzimuthalAngleField->SetTuple1(TupleIterator,SphericalVector[1]);
    }

    // Smooth Scalar fields
    SmoothStructuredPoints::SelfSmoothScalarField(
            KernelSize,
            GridResolution,
            PolarAngleField);

    SmoothStructuredPoints::SelfSmoothScalarField(
            KernelSize,
            GridResolution,
            AzimuthalAngleField);

    // Return data back to Cartesian Coordinates
    for(unsigned int TupleIterator = 0; TupleIterator < NumberOfTuples; TupleIterator++)
    {
        #ifdef _OPENMP  // Thread safe region
        // Get smoothed values
        double Theta[1];
        double Phi[1];
        PolarAngleField->GetTuple(TupleIterator,Theta);
        AzimuthalAngleField->GetTuple(TupleIterator,Phi);

        // Convert back to Cartesian coordinates
        double SmoothedSphericalVector[2] = {Theta[0],Phi[0]};

        #else  // Not thread safe region
        // Get smoothed values
        double Theta = PolarAngleField->GetTuple1(TupleIterator);
        double Phi = AzimuthalAngleField->GetTuple1(TupleIterator);

        // Convert back to Cartesian coordinates
        double SmoothedSphericalVector[2] = {Theta,Phi};
        #endif

        double SmoothedCartesianVector[3];
        SmoothingHelper::ConvertSphericalToCartesianCoordinates(
                SmoothedSphericalVector,
                SmoothedCartesianVector);

        // Store smoothed cartesian vector
        OutputVectorField->SetTuple(TupleIterator,SmoothedCartesianVector);
    }
}

// =============================
// Self Smooth Unit Vector Field
// =============================

void SmoothStructuredPoints::SelfSmoothUnitVectorField(
        unsigned int KernelSize,
        int *GridResolution,
        vtkDoubleArray *VectorField)   // Input / Output
{
    SmoothStructuredPoints::SmoothUnitVectorField(
            KernelSize,
            GridResolution,
            VectorField,
            VectorField);
}

// ===========================
// Smooth Vector Field Use VTK
// ===========================

// Description:
// This filter over-writes the anchoral vector field.
// Also the length of vectors are not preserved.
// The method uses VTK's smoothing class. Both StandardDeviation and Radius Factor
// should be specified.

void SmoothStructuredPoints::SmoothVectorFieldUseVTK(
        vtkDataSet *InputDataSet,
        vtkDoubleArray *VectorField) // Input and Output
{
    // Copy only the structure of dataset to a newer dataset
    vtkSmartPointer<vtkDataSet> CloneInputDataSet = InputDataSet->NewInstance();
    CloneInputDataSet->CopyStructure(InputDataSet);

    // Decompose Vector field to a single component scalar fields
    unsigned int NumberOfComponents = VectorField->GetNumberOfComponents();
    if(NumberOfComponents == 0)
    {
        std::cerr << "Vector field has no component." << std::endl;
        exit(0);
    }

    unsigned int NumberOfTuples = VectorField->GetNumberOfTuples();
    if(NumberOfTuples < 1)
    {
        std::cerr << "Vector field does not have any tuples." << std::endl;
        exit(0);
    }

    for(unsigned int ComponentIterator = 0; ComponentIterator < NumberOfComponents; ComponentIterator++)
    {
        // Declare a single component array
        vtkSmartPointer<vtkDoubleArray> SingleComponentArray = vtkSmartPointer<vtkDoubleArray>::New();
        SingleComponentArray->SetNumberOfComponents(1);
        SingleComponentArray->SetNumberOfTuples(NumberOfTuples);
        SingleComponentArray->SetName("SingleComponentArray");

        // Use Smoothing Filter
        vtkSmartPointer<vtkImageGaussianSmooth> SmoothFilter = vtkSmartPointer<vtkImageGaussianSmooth>::New();
        #if VTK_MAJOR_VERSION <= 5
            SmoothFilter->SetInput(CloneInputDataSet);
        #else
            SmoothFilter->SetInputData(CloneInputDataSet);
        #endif
        SmoothFilter->SetStandardDeviation(0);
        SmoothFilter->SetRadiusFactor(0);
        SmoothFilter->Update();

        // Smoothed Array
        vtkSmartPointer<vtkDoubleArray> SmoothedSingleComponentArray = 
            vtkDoubleArray::SafeDownCast(CloneInputDataSet->GetPointData()->GetArray("SingleComponentArray"));

        // Test
        for(unsigned int TupleIterator = 0; TupleIterator < NumberOfTuples; TupleIterator++)
        {
            if(VectorField->GetComponent(TupleIterator,ComponentIterator) - 
                    SmoothedSingleComponentArray->GetComponent(TupleIterator,ComponentIterator))
            {
                std::cout << VectorField->GetComponent(TupleIterator,ComponentIterator) << ", " 
                    << SmoothedSingleComponentArray->GetComponent(TupleIterator,ComponentIterator) << std::endl;
            }
        }
    }
}


// =============================
// Smooth Symmetric Tensor Field
// =============================

void SmoothStructuredPoints::SmoothSymmetricTensorField(
        unsigned int KernelSize,
        int *GridResolution,
        vtkDoubleArray *InputTensorField,
        vtkDoubleArray *OutputTensorField)   // Output
{
    // Check Input
    if(InputTensorField == NULL)
    {
        std::cerr << "InputTensorField is NULL." << std::endl;
        exit(0);
    }
    else if(InputTensorField->GetNumberOfTuples() < 1)
    {
        std::cerr << "InputTensorField has no tuples." << std::endl;
        exit(0);
    }
    else if(InputTensorField->GetNumberOfComponents() < 1)
    {
        std::cerr << "InputTensorField has no components." << std::endl;
        exit(0);
    }

    unsigned int NumberOfComponents = InputTensorField->GetNumberOfComponents();
    unsigned int NumberOfTuples = InputTensorField->GetNumberOfTuples();
    unsigned int TensorSize = static_cast<unsigned int>(
            floor(0.5+sqrt(static_cast<double>(NumberOfComponents))));

    // Iterate over each element of tensor
    for(unsigned int ColumnIterator = 0; ColumnIterator < TensorSize; ColumnIterator++)
    {
        for(unsigned int RowIterator = ColumnIterator; RowIterator < TensorSize; RowIterator++)
        {
            // Get component index
            unsigned int ComponentIndex = ColumnIterator + RowIterator * TensorSize;
            unsigned int TransposedComponentIndex = RowIterator + ColumnIterator * TensorSize;

            // Copy Tensor element into a separate array
            vtkSmartPointer<vtkDoubleArray> ScalarField = vtkSmartPointer<vtkDoubleArray>::New();
            ScalarField->SetNumberOfComponents(1);
            ScalarField->SetNumberOfTuples(NumberOfTuples);

            // Copy Tensor element to Scalar array
            for(unsigned int TupleIterator = 0; TupleIterator < NumberOfTuples; TupleIterator++)
            {
                ScalarField->SetTuple1(TupleIterator,InputTensorField->GetComponent(TupleIterator,ComponentIndex));
            }

            // Smooth Scalar Array
            SmoothStructuredPoints::SelfSmoothScalarField(
                    KernelSize,
                    GridResolution,
                    ScalarField);

            // Store smoothed field in new Tensor's component
            for(unsigned int TupleIterator = 0; TupleIterator < NumberOfTuples; TupleIterator++)
            {
                #ifdef _OPENMP  // Thread safe region
                double Value[1];
                ScalarField->GetTuple(TupleIterator,Value);
                OutputTensorField->SetComponent(TupleIterator,ComponentIndex,Value[0]);

                // Set value of transposed part for symmetric tensor
                if(RowIterator != ColumnIterator)
                {
                    OutputTensorField->SetComponent(TupleIterator,TransposedComponentIndex,Value[0]);
                }

                #else  // not thread safe region
                double Value = ScalarField->GetTuple1(TupleIterator);
                OutputTensorField->SetComponent(TupleIterator,ComponentIndex,Value);

                // Set value of transposed part for symmetric tensor
                if(RowIterator != ColumnIterator)
                {
                    OutputTensorField->SetComponent(TupleIterator,TransposedComponentIndex,Value);
                }

                #endif
            }
        }
    }
}

// ==================================
// Self Smooth Symmetric Tensor Field
// ==================================

void SmoothStructuredPoints::SelfSmoothSymmetricTensorField(
        unsigned int KernelSize,
        int *GridResolution,
        vtkDoubleArray *TensorField)   // Input / Output
{
    SmoothStructuredPoints::SmoothSymmetricTensorField(
            KernelSize,
            GridResolution,
            TensorField,
            TensorField);
}

// ==========================================
// Convert Cartesian To Spherical Coordinates
// ==========================================

// Description:
// This method is used for smoothing unit vector field.

void SmoothingHelper::ConvertCartesianToSphericalCoordinates(
        double *CartesianVector,
        double SphericalVector[2])
{
    double x = CartesianVector[0];
    double y = CartesianVector[1];
    double z = CartesianVector[2];

    double Phi = acos(z);
    double R = sqrt(x*x + y*y);
    double Theta = acos(x/R);

    if(y < 0)
    {
        Theta = 2 * vtkMath::Pi() - Theta;
    }

    SphericalVector[0] = Theta;
    SphericalVector[1] = Phi;
}

// ==========================================
// Convert Spherical To Cartesian Coordinates
// ==========================================

// Description:
// This method is used for smoothing unit vector field.

void SmoothingHelper::ConvertSphericalToCartesianCoordinates(
        double *SphericalVector,
        double CartesianVector[3])
{
    double Theta = SphericalVector[0];
    double Phi = SphericalVector[1];

    double x = sin(Phi) * cos(Theta);
    double y = sin(Phi) * sin(Theta);
    double z = cos(Phi);

    CartesianVector[0] = x;
    CartesianVector[1] = y;
    CartesianVector[2] = z;
}

// ===============================
// Smooth Gradient Of Scalar Field
// ===============================

void SmoothStructuredPoints::SmoothGradientOfScalarField(
        unsigned int KernelSize,
        int *GridResolution,
        vtkDoubleArray *InputScalarField,
        vtkDoubleArray *OutputGradientField)
{
    // Check Input
    if(InputScalarField == NULL)
    {
        std::cerr << "InputScalarField is NULL." << std::endl;
        exit(0);
    }
    else if(InputScalarField->GetNumberOfTuples() < 1)
    {
        std::cerr << "InputScalalrField has no tuples." << std::endl;
        exit(0);
    }
    else if(InputScalarField->GetNumberOfComponents() != 1)
    {
        std::cerr << "InputScalarField is not scalar." << std::endl;
        exit(0);
    }

    // Set Output
    unsigned int NumberOfTuples = InputScalarField->GetNumberOfTuples();
    OutputGradientField->SetNumberOfComponents(DIMENSION);
    OutputGradientField->SetNumberOfTuples(NumberOfTuples);
    OutputGradientField->SetName("SmoothedGradientField");

    // Initialize Smoothing Kernel
    double *SmoothingKernel = new double[KernelSize];
    SmoothingHelper::InitializeSmoothingKernel(
            KernelSize,
            SmoothingKernel);   // Output

    // Initialize Smoothing Derivative Kernel
    double *SmoothingGradientKernel = new double[KernelSize];
    SmoothingHelper::InitializeSmoothingGradientKernel(
            KernelSize,
            SmoothingGradientKernel);   // Output

    // Iterate over each dimension for taking gradient
    for(unsigned int GradientIterator = 0; GradientIterator < DIMENSION; GradientIterator++)
    {
        // Declare Output Scalar Field as gradient component
        vtkSmartPointer<vtkDoubleArray> OutputScalarField = vtkSmartPointer<vtkDoubleArray>::New();
        OutputScalarField->SetNumberOfComponents(1);
        OutputScalarField->SetNumberOfTuples(NumberOfTuples);

        // Declare intermediate scalar fields
        vtkDoubleArray *ScalarField[DIMENSION+1];
        ScalarField[0] = InputScalarField;
        ScalarField[DIMENSION] = OutputScalarField;

        // Intermediate scalar fields
        for(unsigned int FieldIterator = 1; FieldIterator < DIMENSION; FieldIterator++)
        {
            ScalarField[FieldIterator] = vtkDoubleArray::New();
        }

        // Smoothing over each direction
        for(unsigned int DimensionIterator = 0; DimensionIterator < DIMENSION; DimensionIterator++)
        {
            // Assign Input and Output Scalar field for convolution
            vtkDoubleArray *OriginalScalarField = ScalarField[DimensionIterator];
            vtkDoubleArray *ConvolutedScalarField = ScalarField[DimensionIterator+1];

            // Treat the gradient direction differently
            double *Kernel = NULL;
            if(DimensionIterator == GradientIterator)
            {
                // Use Smoothing Derivative Kernel
                Kernel = SmoothingGradientKernel;
            }
            else
            {
                // Use smoothing Kernel
                Kernel = SmoothingKernel;
            }

            // Convolute Kernel
            SmoothingHelper::ConvoluteKernelInOneDirection(
                    Kernel,
                    KernelSize,
                    DimensionIterator,
                    GridResolution,
                    OriginalScalarField,
                    ConvolutedScalarField);
        }

        // Copy the Output Scalar Field into the Gradient component
        for(unsigned int TupleIterator = 0; TupleIterator < NumberOfTuples; TupleIterator++)
        {
            #ifdef _OPENMP  // Thread safe region
            double Value[1];
            OutputScalarField->GetTuple(TupleIterator,Value);
            OutputGradientField->SetComponent(TupleIterator,GradientIterator,Value[0]);

            #else  // Not thread safe region
            double Value = OutputScalarField->GetTuple1(TupleIterator);
            OutputGradientField->SetComponent(TupleIterator,GradientIterator,Value);
            #endif
        }

        // Delete Scalar Field
        for(unsigned int FieldIterator = 1; FieldIterator < DIMENSION; FieldIterator++)
        {
            ScalarField[FieldIterator]->Delete();
        }
    }

    // Free memory
    delete [] SmoothingKernel;
    delete [] SmoothingGradientKernel;
}

// =======================
// Smooth Hessian At Point
// =======================

// Description:
// This method does not allocate the output. The HessianAtPoint array should be
// allocated before calling this function. The hessian is calculated only at one point.

void SmoothStructuredPoints::SmoothHessianAtPoint(
        unsigned int KernelSize,
        int *GridResolution,
        int PointId,
        vtkDoubleArray *InputScalarField,
        double **HessianAtPoint)   // Output
{
    // Initialize Smoothing Kernel
    double *SmoothingKernel = new double[KernelSize];
    SmoothingHelper::InitializeSmoothingKernel(
            KernelSize,
            SmoothingKernel);   // Output

    // Initialize Smoothing Gradient Kernel
    double *SmoothingGradientKernel = new double[KernelSize];
    SmoothingHelper::InitializeSmoothingGradientKernel(
            KernelSize,
            SmoothingGradientKernel);   // Output

    // Initializ Smoothing Hessian Kernel
    double *SmoothingHessianKernel = new double[KernelSize];
    SmoothingHelper::InitializeSmoothingHessianKernel(
            KernelSize,
            SmoothingHessianKernel);

    // Create three kernels, each 1 dimensional, for hessian components
    double **Kernels = new double*[DIMENSION];

    // Iterate over Rows of Tensor
    for(unsigned int RowIterator = 0; RowIterator < DIMENSION; RowIterator++)
    {
        // Iterate over columns of tensor
        for(unsigned int ColumnIterator = 0; ColumnIterator < DIMENSION; ColumnIterator++)
        {
            // Pre-assign all kernels to be smoothing non-derivative kernel
            for(int Direction = 0; Direction < DIMENSION; Direction++)
            {
                Kernels[Direction] = SmoothingKernel;
            }

            // If two derivatives in one direction, use Hessian Kernel
            if(RowIterator == ColumnIterator)
            {
                Kernels[RowIterator] = SmoothingHessianKernel;
            }
            else
            {
                // Each direction has a first order derivative, use gradient kernel for each
                Kernels[RowIterator] = SmoothingGradientKernel;
                Kernels[ColumnIterator] = SmoothingGradientKernel;
            }

            // Compute convolution
            HessianAtPoint[RowIterator][ColumnIterator] = SmoothingHelper::Convolute3DKernelAtPoint(
                    KernelSize,
                    Kernels,
                    GridResolution,
                    PointId,
                    InputScalarField);
        }
    }

    // Free memory
    delete [] SmoothingKernel;
    delete [] SmoothingGradientKernel;
    delete [] SmoothingHessianKernel;
    delete [] Kernels;
}

// ============================
// Convolute 3D Kernel At Point
// ============================

// Description:
// Performs 3-dimensional comvolution at once for a specific point in the grid.
// The difference with other functions is that the kernel is not separated, so the 
// whole 3-dimensional kernel operates on data. The result is also not convolved over 
// all grid points, but only for a specific point.
// Kernels is an array that i-th row is a kernel in the i-th direction. Each kernel is 
// 1-dimensional kernel. 

double SmoothingHelper::Convolute3DKernelAtPoint(
        unsigned int KernelSize,
        double **Kernels,
        int *GridResolution,
        int PointId,
        vtkDoubleArray *ScalarField)
{
    // Main Point Index
    int MainPointIndex[DIMENSION];
    StructuredPointsHelper::GetPointIndex(GridResolution,PointId,MainPointIndex);

    // Kernel Radius
    int KernelRadius = static_cast<int>((KernelSize-1) / 2);
    
    // Convolution Point Index
    int ConvolutionPointIndex[DIMENSION];

    // Sum of all convoluted points form kernel
    double Sum = 0;

    // Kernel Index
    int KernelIndex[DIMENSION];

    // Iterate over all 3D kernel elements
    for(int XIterator = -KernelRadius; XIterator <= KernelRadius; XIterator++)
    {
        KernelIndex[0] = XIterator + KernelRadius;
        ConvolutionPointIndex[0] = MainPointIndex[0] + XIterator;
        StructuredPointsHelper::ClipIndexToBounds(GridResolution[0],ConvolutionPointIndex[0]);

        // Iterate over y
        for(int YIterator = -KernelRadius; YIterator <= KernelRadius; YIterator++)
        {
            KernelIndex[1] = YIterator + KernelRadius;
            ConvolutionPointIndex[1] = MainPointIndex[1] + YIterator;
            StructuredPointsHelper::ClipIndexToBounds(GridResolution[1],ConvolutionPointIndex[1]);

            // Iterate over z
            for(int ZIterator = -KernelRadius; ZIterator <= KernelRadius; ZIterator++)
            {
                KernelIndex[2] = ZIterator + KernelRadius;
                ConvolutionPointIndex[2] = MainPointIndex[2] + ZIterator;
                StructuredPointsHelper::ClipIndexToBounds(GridResolution[2],ConvolutionPointIndex[2]);

                // Compute 3D Kernel Value from 1D kernels  at convolution point
                double KernelValue = Kernels[0][KernelIndex[0]] * 
                                     Kernels[1][KernelIndex[1]] *
                                     Kernels[2][KernelIndex[2]];

                // Get Convolution PointId
                int ConvolutionPointId = StructuredPointsHelper::GetPointId(GridResolution,ConvolutionPointIndex);

                #ifdef _OPENMP  // Thread safe region

                // Get convolution point scalar data
                double ScalarData[1];
                ScalarField->GetTuple(ConvolutionPointId,ScalarData);
                
                // Convolute
                Sum += KernelValue * ScalarData[0];

                #else  // Not thread safe region

                // Get convolution point scalar data
                double ScalarData = ScalarField->GetTuple1(ConvolutionPointId);
                
                // Convolute
                Sum += KernelValue * ScalarData;
                #endif
            }
        }
    }
   
    return Sum;
}
