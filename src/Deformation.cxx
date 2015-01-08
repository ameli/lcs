/*
 * =====================================================================================
 *
 *       Filename:  Deformation.cxx
 *
 *    Description:  Computes Deformation of flowmap
 *
 *        Version:  1.0
 *        Created:  02/12/2014 05:19:54 PM
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

#include "Deformation.h"
#include "SmoothStructuredPoints.h"

// Pipeline
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDataObject.h>
#include <vtkObjectFactory.h>
#include <vtkDemandDrivenPipeline.h>

// Data
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkImageData.h>

// Arrays
#include <vtkDataArray.h>
#include <vtkDataArrayCollection.h>
#include <vtkDoubleArray.h>

// General
#include <vtkSmartPointer.h>
#include <vtkCallbackCommand.h>
#include <vector>
#include <algorithm>  // for std::sort
#include <cmath>      // for sqrt, log
#include <sstream>    // for stringstream
#ifdef _OPENMP
#include <omp.h>
#endif

// Algorithms
#include <vtkMath.h>

// Debug
#include <time.h>

// ======
// Macros
// ======

vtkStandardNewMacro(Deformation);
vtkCxxRevisionMacro(Deformation,"$Revision 1.0$");
#define DIMENSION 3

// ===========
// Constructor
// ===========

Deformation::Deformation()
{
    // Member Data
    this->SmoothingStatus = true;
    this->SmoothingKernelSize = 7;
    this->ProgressStatus = false;

    // Callbacks
    vtkSmartPointer<vtkCallbackCommand> ProgressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
    ProgressCallback->SetCallback(this->ProgressFunction);
    this->AddObserver(vtkCommand::ProgressEvent,ProgressCallback);
}

// ==========
// Destructor
// ==========

Deformation::~Deformation()
{
}

// ==========
// Print Self
// ==========

void Deformation::PrintSelf(ostream &os,vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}

// ============================
// Fill Output Port Information
// ============================

int Deformation::FillOutputPortInformation(int port,vtkInformation *info)
{
    if(port == 0)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(),"vtkStructuredPoints");
        return 1;
    }

    return 0;
}

// =================
// Progress Function
// =================

void Deformation::ProgressFunction(
        vtkObject *Caller,
        unsigned long int vtkNotUsed(EventId),
        void *vtkNotUsed(ClientData),
        void *vtkNotUsed(CallData))
{
    Deformation *Filter = Deformation::SafeDownCast(Caller);
    int ProgressValuePercent = 10*floor(0.5+10*Filter->GetProgress());

    if(Filter->GetProgressStatus() == true && ProgressValuePercent >= 10)
    {
        std::cout << "Deformation Filter Progress: "
                  << std::setw(3) << std::right 
                  << ProgressValuePercent << " %"<< std::endl;
    }
}

// ======================
// Filter Update Progress
// ======================

void Deformation::FilterUpdateProgress(
        unsigned int Step,
        unsigned int NumberOfSteps)
{
    unsigned int IntervalPercent = 10;

    // Update Progress
    if(this->ProgressStatus == true)
    {
        if(Step % static_cast<int>(NumberOfSteps/IntervalPercent) == 0 ||
           Step == NumberOfSteps-1)
        {
            double ProgressValue = static_cast<double>(Step)/(NumberOfSteps-1);
            this->UpdateProgress(ProgressValue);
        }
    }
}

// ============
// Request Data
// ============

int Deformation::RequestData(
        vtkInformation *vtkNotUsed(Request),
        vtkInformationVector **InputVector,
        vtkInformationVector *OutputVector)
{
    // Input
    vtkInformation *InputInfo = InputVector[0]->GetInformationObject(0);
    vtkDataObject *InputDataObject = InputInfo->Get(vtkDataObject::DATA_OBJECT());
    vtkDataSet *InputDataSet = vtkDataSet::SafeDownCast(InputDataObject);

    // Output
    vtkInformation *OutputInfo = OutputVector->GetInformationObject(0);
    vtkDataObject *OutputDataObject = OutputInfo->Get(vtkDataObject::DATA_OBJECT());
    vtkDataSet *OutputDataSet = vtkDataSet::SafeDownCast(OutputDataObject);

    // Compute EigenValues, EigenVectors
    vtkSmartPointer<vtkDoubleArray> EigenValues = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDataArrayCollection> EigenVectorsCollection = vtkSmartPointer<vtkDataArrayCollection>::New();
    vtkSmartPointer<vtkDoubleArray> EigenVector[DIMENSION];

    for(unsigned int EigenVectorsIterator = 0; EigenVectorsIterator < DIMENSION; EigenVectorsIterator++)
    {
        EigenVector[EigenVectorsIterator] = vtkDoubleArray::New();
        EigenVectorsCollection->AddItem(EigenVector[EigenVectorsIterator]);
    }

    this->ComputeEigenValuesAndVectors(
            InputDataSet,
            EigenValues,               // Output
            EigenVectorsCollection);   // Output

    // Declare Deformation variables
    vtkSmartPointer<vtkDoubleArray> MaxStrainValues  = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> MaxStrainVectors = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> MinStrainValues  = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> MinStrainVectors = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> ShearValues      = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> ShearVectors1    = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> ShearVectors2    = vtkSmartPointer<vtkDoubleArray>::New();

    // Compute Deformation Values and Vectors
    this->ComputeDeformationValuesAndVectors(
            EigenValues,
            EigenVectorsCollection,
            MaxStrainValues,        // Output
            MaxStrainVectors,       // Output
            MinStrainValues,        // Output
            MinStrainVectors,       // Output
            ShearValues,            // Output
            ShearVectors1,          // Output
            ShearVectors2);         // Output

    // Clear Memory
    for(unsigned int EigenVectorsIterator = 0;
        EigenVectorsIterator < DIMENSION;
        EigenVectorsIterator++)
    {
        EigenVector[EigenVectorsIterator]->Delete();
    }

    // Generate Output
    this->GenerateOutputDataSet(
            InputDataSet,
            OutputDataSet,    // Output
            MaxStrainValues,
            MaxStrainVectors,
            MinStrainValues,
            MinStrainVectors,
            ShearValues,
            ShearVectors1,
            ShearVectors2);

    // Set Output
    OutputInfo->Set(vtkDataObject::DATA_OBJECT(),OutputDataSet);

    return 1;
}

// ===============================
// Compute EigenValues And Vectors
// ===============================

void Deformation::ComputeEigenValuesAndVectors(
        vtkDataSet *InputDataSet,
        vtkDoubleArray *EigenValues,  // Output
        vtkDataArrayCollection *EigenVectorsCollection)  // Output
{
    // Check Input DataSet
    if(InputDataSet == NULL)
    {
        vtkErrorMacro("InputDataSet is NULL.");
    }
    else if(InputDataSet->GetNumberOfPoints() < 1)
    {
        vtkErrorMacro("InputDataSet has no points.");
    }

    // Get Tensor attributes
    vtkSmartPointer<vtkDoubleArray> TensorsArray = vtkDoubleArray::SafeDownCast(InputDataSet->GetPointData()->GetTensors("CG"));

    // Check Tensors
    if(TensorsArray == NULL)
    {
        vtkErrorMacro("Tensors array is NULL.");
    }
    else if(TensorsArray->GetNumberOfTuples() < 2)
    {
        vtkErrorMacro("Tensors does not have any tuples.");
    }
    else if(TensorsArray->GetNumberOfComponents() != DIMENSION*DIMENSION)
    {
        vtkErrorMacro("TensorArray does not have enough components.");
    }

    // Smooth Tensor Field
    if(this->SmoothingStatus == true)
    {
        // Get Grid Resolution
        vtkImageData *InputImageData = vtkImageData::SafeDownCast(InputDataSet);
        if(InputImageData == NULL)
        {
            vtkErrorMacro("InputImageData is NULL.");
        }
        else if(InputImageData->GetNumberOfPoints() < 1)
        {
            vtkErrorMacro("InputImageData has no points.");
        }

        // int *GridResolution = InputImageData->GetDimensions();
        int GridResolution[DIMENSION];
        InputImageData->GetDimensions(GridResolution);

        // Smooth Tensor
        // SmoothStructuredPoints::SelfSmoothSymmetricTensorField(
        //     this->SmoothingKernelSize,
        //     GridResolution,
        //     TensorsArray);
    }

    // Allocate Memory for Outputs //
    unsigned int NumberOfPoints = TensorsArray->GetNumberOfTuples();

    // EigenValues
    EigenValues->SetNumberOfComponents(DIMENSION);
    EigenValues->SetNumberOfTuples(NumberOfPoints);
    EigenValues->SetName("EigenValues");

    // EigenVectors
    vtkSmartPointer<vtkDoubleArray> EigenVector[DIMENSION];

    // Get pointer to each of EigenVectors collections
    for(unsigned int EigenVectorsIterator = 0;
        EigenVectorsIterator < DIMENSION;
        EigenVectorsIterator++)
    {
        EigenVector[EigenVectorsIterator] = vtkDoubleArray::SafeDownCast(
                EigenVectorsCollection->GetItem(EigenVectorsIterator));
        EigenVector[EigenVectorsIterator]->SetNumberOfComponents(DIMENSION);
        EigenVector[EigenVectorsIterator]->SetNumberOfTuples(NumberOfPoints);

        std::stringstream ArrayName;
        ArrayName << "EigenVector" << EigenVectorsIterator;
        EigenVector[EigenVectorsIterator]->SetName(ArrayName.str().c_str());
    }

    // Loop over points
    double Tensor[DIMENSION * DIMENSION];
    unsigned int ProgressCounter = 0;

    #ifdef _OPENMP
    #pragma omp parallel for \
    default (none) \
    shared(NumberOfPoints,TensorsArray,EigenValues,EigenVector,ProgressCounter) \
    private(Tensor)
    #endif
    for(unsigned int PointIterator = 0; PointIterator < NumberOfPoints; PointIterator++)
    {
        // Get tensor for a point
        TensorsArray->GetTuple(PointIterator,Tensor);

        // Convert 1D array to a square matrix
        double **TensorMatrix = DeformationHelper::CreateMatrix(DIMENSION,DIMENSION);
        for(unsigned int Row = 0; Row < DIMENSION; Row++)
        {
            for(unsigned int Column = 0; Column < DIMENSION; Column++)
            {
                TensorMatrix[Row][Column] = double(Tensor[Column + Row*DIMENSION]);
            }
        }

        // Diagonalize Tensor Matrix
        double TensorEigenValues[DIMENSION];
        double **TensorEigenVectors = DeformationHelper::CreateMatrix(DIMENSION,DIMENSION);

        // Note: Originl matrix will be over written to lower triangular matrix);
        vtkMath::Jacobi(TensorMatrix,TensorEigenValues,TensorEigenVectors); 

        // Transpose EigenVectors array (for easier operation on array)
        // NOTE: From now, in all code the EigenVector array is transposed.
        DeformationHelper::TransposeSquareMatrix(TensorEigenVectors,DIMENSION);
 
        // Sort EigenValues and Vectors
        this->SortEigenValuesAndVectors(
                &TensorEigenValues[0],   // Input and Output
                TensorEigenVectors,  // Input and Output
                DIMENSION);

        // Check Negative EigenValues
        for(int DimensionIterator = 0; DimensionIterator < DIMENSION; DimensionIterator++)
        {
            if(TensorEigenValues[DimensionIterator] < 0)
            {
                // Allow tolerance
                if(TensorEigenValues[DimensionIterator] > -1e-5)
                {
                    TensorEigenValues[DimensionIterator] = 0.0;
                }
                else
                {
                    vtkErrorMacro("Negative eigenvalue detected.");
                }
            }
        }
        // Write to Output
        EigenValues->SetTuple(PointIterator,TensorEigenValues);
        for(unsigned int EigenVectorsIterator = 0;
            EigenVectorsIterator < DIMENSION;
            EigenVectorsIterator++)
        {
            EigenVector[EigenVectorsIterator]->SetTuple(PointIterator,TensorEigenVectors[EigenVectorsIterator]);
        }

        // Free Matrix
        DeformationHelper::FreeMatrix(TensorMatrix,DIMENSION);
        DeformationHelper::FreeMatrix(TensorEigenVectors,DIMENSION);

        // Update Progress
        #ifdef _OPENMP
        #pragma omp atomic
        #endif
        ProgressCounter++;
        this->FilterUpdateProgress(ProgressCounter,NumberOfPoints);
    }
}

// ============================
// Sort EigenValues and Vectors
// ============================

// Description:
// Sort eigenvalues in accesnding order. Therefore EigenValues[0] become smallest
// and EigenValues[DIMENSION-1] become the largest. EigenVectors are sorted according
// to their associated eigenvalues.
// Note: This function over-write the input argument arrays. Thus inputs can not be
// constant arrays.

void Deformation::SortEigenValuesAndVectors(
        double *TensorEigenValues,    // Input and Output
        double **TensorEigenVectors,  // Input and Output
        unsigned int TensorDimension)
{
}

// ======================================
// Compute Deformation Values And Vectors
// ======================================

// Description:
// Assumed that EigenValues are sorted in accending manner.

void Deformation::ComputeDeformationValuesAndVectors(
        vtkDoubleArray *EigenValues,
        vtkDataArrayCollection *EigenVectorsCollection,
        vtkDoubleArray *MaxStrainValues,   // Output
        vtkDoubleArray *MaxStrainVectors,  // Output
        vtkDoubleArray *MinStrainValues,   // Output
        vtkDoubleArray *MinStrainVectors,  // Output
        vtkDoubleArray *ShearValues,       // Output
        vtkDoubleArray *ShearVectors1,     // Output
        vtkDoubleArray *ShearVectors2)     // Output
{
    // Number of Points
    unsigned int NumberOfPoints = EigenValues->GetNumberOfTuples();

    // Check EigenValues
    if(EigenValues == NULL)
    {
        vtkErrorMacro("EigenValues is NULL.");
    }
    else if(EigenValues->GetNumberOfTuples() < 1)
    {
        vtkErrorMacro("EigenValues array have no tuple.");
    }

    // Check EigenVectors Collection
    if(EigenVectorsCollection == NULL)
    {
        vtkErrorMacro("EigenVectors collection is NULL.")
    }
    else if(EigenVectorsCollection->GetNumberOfItems() != DIMENSION)
    {
        vtkErrorMacro("EigenVectors collection does not have enough items.");
    }

    // Get EigenVector pointers
    vtkSmartPointer<vtkDoubleArray> EigenVector[DIMENSION];

    for(unsigned int EigenVectorsIterator = 0; EigenVectorsIterator < DIMENSION; EigenVectorsIterator++)
    {
        EigenVector[EigenVectorsIterator] = vtkDoubleArray::SafeDownCast(
                EigenVectorsCollection->GetItem(EigenVectorsIterator));

        // Check EigenVector
        if(EigenVector[EigenVectorsIterator] == NULL)
        {
            vtkErrorMacro("EigenVector is NULL.");
        }
        else if(EigenVector[EigenVectorsIterator]->GetNumberOfTuples() != 
                static_cast<int>(NumberOfPoints))
        {
            vtkErrorMacro("EigenVector does not have corrent number of points.");
        }
    }

    // Allocate Deformation Variables //

    // Min Strain Values
    MinStrainValues->SetNumberOfComponents(1);
    MinStrainValues->SetNumberOfTuples(NumberOfPoints);
    MinStrainValues->SetName("MinStrainValues");

    // Min Strain Vectors
    MinStrainVectors->SetNumberOfComponents(DIMENSION);
    MinStrainVectors->SetNumberOfTuples(NumberOfPoints);
    MinStrainVectors->SetName("MinStrainVectors");

    // Max Strain Values
    MaxStrainValues->SetNumberOfComponents(1);
    MaxStrainValues->SetNumberOfTuples(NumberOfPoints);
    MaxStrainValues->SetName("MaxStrainValues");

    // Max Strain Vectors
    MaxStrainVectors->SetNumberOfComponents(DIMENSION);
    MaxStrainVectors->SetNumberOfTuples(NumberOfPoints);
    MaxStrainVectors->SetName("MaxStrainVectors");

    // Shear Values
    ShearValues->SetNumberOfComponents(1);
    ShearValues->SetNumberOfTuples(NumberOfPoints);
    ShearValues->SetName("ShearValues");

    // Shear Vectors 1
    ShearVectors1->SetNumberOfComponents(DIMENSION);
    ShearVectors1->SetNumberOfTuples(NumberOfPoints);
    ShearVectors1->SetName("ShearVectors1");

    // Shear Vectors 2
    ShearVectors2->SetNumberOfComponents(DIMENSION);
    ShearVectors2->SetNumberOfTuples(NumberOfPoints);
    ShearVectors2->SetName("ShearVectors2");

    // Loop over all points
    for(unsigned int PointIterator = 0; PointIterator < NumberOfPoints; PointIterator++)
    {
        // Get data at the point
        double MinEigenValueAtPoint   = EigenValues->GetComponent(PointIterator,0);
        double MaxEigenValueAtPoint   = EigenValues->GetComponent(PointIterator,DIMENSION-1);
        double *MinEigenVectorAtPoint = EigenVector[0]->GetTuple(PointIterator); // RISKY
        double *MaxEigenVectorAtPoint = EigenVector[DIMENSION-1]->GetTuple(PointIterator); // RISKY

        // Min Strain Values
        double MinStrainValueAtPoint = sqrt(MinEigenValueAtPoint);
        MinStrainValues->SetTuple1(PointIterator,MinStrainValueAtPoint);

        // Min Strain Vectors
        MinStrainVectors->SetTuple(PointIterator,MinEigenVectorAtPoint);

        // Max Strain Values
        double MaxStrainValueAtPoint = sqrt(MaxEigenValueAtPoint);
        MaxStrainValues->SetTuple1(PointIterator,MaxStrainValueAtPoint);

        // Max Strain Vectors
        MaxStrainVectors->SetTuple(PointIterator,MaxEigenVectorAtPoint);

        // Shear Values
        double ShearValueAtPoint = MaxStrainValueAtPoint - MinStrainValueAtPoint;
        ShearValues->SetTuple1(PointIterator,ShearValueAtPoint);

        // Shear Vectors
        double MaxCoefficient = sqrt(MaxStrainValueAtPoint / (MaxStrainValueAtPoint + MinStrainValueAtPoint));
        double MinCoefficient = sqrt(MinStrainValueAtPoint / (MaxStrainValueAtPoint + MinStrainValueAtPoint));

        // Note: MultiplyScalar overwrites data
        vtkMath::MultiplyScalar(MaxEigenVectorAtPoint,MaxCoefficient);
        vtkMath::MultiplyScalar(MinEigenVectorAtPoint,MinCoefficient);

        double ShearVector1AtPoint[DIMENSION];
        double ShearVector2AtPoint[DIMENSION];

        vtkMath::Add(MaxEigenVectorAtPoint,MinEigenVectorAtPoint,ShearVector1AtPoint);
        vtkMath::Subtract(MaxEigenVectorAtPoint,MinEigenVectorAtPoint,ShearVector2AtPoint);

        ShearVectors1->SetTuple(PointIterator,ShearVector1AtPoint);
        ShearVectors2->SetTuple(PointIterator,ShearVector2AtPoint);

        // Update Progress
        this->FilterUpdateProgress(PointIterator,NumberOfPoints);
    }
}

// =======================
// Generate Output DataSet
// =======================

void Deformation::GenerateOutputDataSet(
        vtkDataSet *InputDataSet,
        vtkDataSet *OutputDataSet,    // Output
        vtkDoubleArray *MaxStrainValues,
        vtkDoubleArray *MaxStrainVectors,
        vtkDoubleArray *MinStrainValues,
        vtkDoubleArray *MinStrainVectors,
        vtkDoubleArray *ShearValues,
        vtkDoubleArray *ShearVectors1,
        vtkDoubleArray *ShearVectors2)
{
    // Copy structure from input to output
    OutputDataSet->CopyStructure(InputDataSet);

    // Add Deformation Arrays
    OutputDataSet->GetPointData()->AddArray(MaxStrainValues);
    OutputDataSet->GetPointData()->AddArray(MaxStrainVectors);
    OutputDataSet->GetPointData()->AddArray(MinStrainValues);
    OutputDataSet->GetPointData()->AddArray(MinStrainVectors);
    OutputDataSet->GetPointData()->AddArray(ShearValues);
    OutputDataSet->GetPointData()->AddArray(ShearVectors1);
    OutputDataSet->GetPointData()->AddArray(ShearVectors2);
}

// =============
// Create Matrix
// =============

double ** DeformationHelper::CreateMatrix(unsigned int NumberOfRows,unsigned int NumberOfColumns)
{
    // Declaration
    double **Matrix = new double *[NumberOfColumns];

    // Append row arrays to each column
    for(unsigned int ColumnIterator = 0; ColumnIterator < NumberOfColumns; ColumnIterator++)
    {
        Matrix[ColumnIterator] = new double[NumberOfRows];
    }

    return Matrix;
}

// ===========
// Free Matrix
// ===========

void DeformationHelper::FreeMatrix(double **Matrix,unsigned int NumberOfColumns)
{
    // Deleting Columns
    for(unsigned int ColumnIterator = 0; ColumnIterator < NumberOfColumns; ColumnIterator++)
    {
        delete [] Matrix[ColumnIterator];
    }

    // Delete first row of columns
    delete [] Matrix;
}

// =======================
// Transpose Square Matrix
// =======================

void DeformationHelper::TransposeSquareMatrix(double **SquareMatrix,unsigned int MatrixSize)
{
    // Loop over upper half of matrix
    for(unsigned int RowIterator = 0; RowIterator < MatrixSize; RowIterator++)
    {
        for(unsigned int ColumnIterator = RowIterator +1; ColumnIterator < MatrixSize; ColumnIterator++)
        {
            double temp = SquareMatrix[RowIterator][ColumnIterator];
            SquareMatrix[RowIterator][ColumnIterator] = SquareMatrix[ColumnIterator][RowIterator];
            SquareMatrix[ColumnIterator][RowIterator] = temp;
        }
    }
}
