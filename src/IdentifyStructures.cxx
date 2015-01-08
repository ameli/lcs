/*
 * =====================================================================================
 *
 *       Filename:  IdentifyStructures.cxx
 *
 *    Description:  Finds Coherent Structures
 *
 *        Version:  1.0
 *        Created:  02/14/2014 12:21:46 AM
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
#include "SmoothStructuredPoints.h"

// Pipeline
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDataObject.h>
#include <vtkObjectFactory.h>
#include <vtkDemandDrivenPipeline.h>

// Data
#include <vtkDataSet.h>
#include <vtkStructuredData.h>  // for vtkStructuredData::ComputePointId
#include <vtkStructuredPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include "PointDataType.h"

// Elements
// #include <vtkCellIterator.h>
#include <vtkCell.h>
#include <vtkGenericCell.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>

// Arrays
#include <vtkIdList.h>
#include <vtkDataArray.h>
#include <vtkDataArrayCollection.h>
#include <vtkDoubleArray.h>

// General
#include <vtkSmartPointer.h>
#include <vtkCallbackCommand.h>
#include <vector>
#include <algorithm>  // for std::sort
#include <cmath>      // for sqrt
#include <stdlib.h>   // for std::div
#include <sstream>    // for stringstream
#ifdef _OPENMP
#include <omp.h>
#endif

// Algorithms
#include <vtkGradientFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkMath.h>
#include "CubeCell.h"

// Debug
#include <time.h>

// ======
// Macros
// ======

vtkStandardNewMacro(IdentifyStructures);
vtkCxxRevisionMacro(IdentifyStructures,"$Revision 1.0$");

#define DIMENSION 3
#define CUBESIZE 8

// ===========
// Constructor
// ===========

IdentifyStructures::IdentifyStructures()
{
    // Default member data
    this->StructureMode = StructureType::MAX_STRAIN;
    this->SmoothingStatus = true;
    this->SmoothingKernelSize = 7;
    this->DebugStatus = false;
    this->ProgressStatus = false;

    // Callbacks
    vtkSmartPointer<vtkCallbackCommand> ProgressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
    ProgressCallback->SetCallback(this->ProgressFunction);
    this->AddObserver(vtkCommand::ProgressEvent,ProgressCallback);
}

// ==========
// Destructor
// ==========

IdentifyStructures::~IdentifyStructures()
{
}

// ====================
// Accessors / Mutators
// ====================

// Set Structure Mode
void IdentifyStructures::SetStructureMode(int InputStructureMode)
{
    switch(InputStructureMode)
    {
        // Max Strain
        case static_cast<int>(StructureType::MAX_STRAIN):
        {
            this->StructureMode = StructureType::MAX_STRAIN;
            break;
        }

        // Min Strain
        case static_cast<int>(StructureType::MIN_STRAIN):
        {
            this->StructureMode = StructureType::MIN_STRAIN;
            break;
        }

        // Shear
        case static_cast<int>(StructureType::SHEAR):
        {
            this->StructureMode = StructureType::SHEAR;
            break;
        }

        // Undefined mode
        default:
        {
            vtkErrorMacro("Structure mode is undefined.");
        }
    }
}

// Debug On
void IdentifyStructures::DebugOn()
{
    this->DebugStatus = true;
}

// Debug Off
void IdentifyStructures::DebugOff()
{
    this->DebugStatus = false;
}

// ==========
// Print Self
// ==========

void IdentifyStructures::PrintSelf(ostream &os,vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}

// ============================
// Fill Output Port Information
// ============================

int IdentifyStructures::FillOutputPortInformation(int port,vtkInformation *info)
{
    if(port == 0)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(),"vtkPolyData");
        return 1;
    }

    return 0;
}

// ===================
// Request Data Object
// ===================

int IdentifyStructures::RequestDataObject(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **vtkNotUsed(InputInfo),
        vtkInformationVector *OutputVector)
{
    vtkInformation *OutputInfo = OutputVector->GetInformationObject(0);
    vtkPolyData *OutputPolyData = vtkPolyData::SafeDownCast(OutputInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Create new instance
    if(OutputPolyData == NULL)
    {
        OutputPolyData = vtkPolyData::New();
        OutputInfo->Set(vtkDataObject::DATA_OBJECT(),OutputPolyData);
        OutputPolyData->FastDelete();
        this->GetOutputPortInformation(0)->Set(vtkDataObject::DATA_EXTENT_TYPE(),OutputPolyData->GetExtentType());
    }

    return 1;
}

// =================
// Progress Function
// =================

void IdentifyStructures::ProgressFunction(
        vtkObject *Caller,
        unsigned long int vtkNotUsed(EventId),
        void *vtkNotUsed(ClientData),
        void *vtkNotUsed(CallData))
{
    IdentifyStructures *Filter = IdentifyStructures::SafeDownCast(Caller);
    int ProgressValuePercent = 10*floor(0.5+10*Filter->GetProgress());

    if(Filter->GetProgressStatus() == true && ProgressValuePercent >= 10)
    {
        std::cout << "IdentifyStructures Filter Progress: "
                  << std::setw(3) << std::right 
                  << ProgressValuePercent << " %"<< std::endl;
    }
}

// ======================
// Filter Update Progress
// ======================

void IdentifyStructures::FilterUpdateProgress(
        unsigned int Step,
        unsigned int NumberOfSteps)
{
    unsigned int ProgressPercentIncrement = 10;

    // Update Progress
    if(this->ProgressStatus == true)
    {
        if(Step % static_cast<int>(NumberOfSteps/ProgressPercentIncrement) == 0 ||
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

int IdentifyStructures::RequestData(
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
    vtkPolyData *OutputPolyData = vtkPolyData::SafeDownCast(OutputDataObject);

    // Define Deformation Value and Vector
    vtkDoubleArray *DeformationValue = NULL;
    vtkDoubleArray *DeformationVector = NULL;

    // Get Deformation Value and Vector from InputDataSet
    this->GetDeformationValueAndVector(InputDataSet,DeformationValue,DeformationVector);

    // Input Grid Resolution
    vtkStructuredPoints *InputStructuredPoints = vtkStructuredPoints::SafeDownCast(InputDataSet);
    int *GridResolution = InputStructuredPoints->GetDimensions();

    // Smooth Deformation Vectors
    // if(this->SmoothingStatus == true)
    // {
    //     SmoothStructuredPoints::SelfSmoothUnitVectorField(
    //             this->SmoothingKernelSize,
    //             GridResolution,
    //             DeformationVector);
    // }

    // Compute Gradient of Deformation Value
    vtkSmartPointer<vtkDoubleArray> DeformationValueGradient = vtkSmartPointer<vtkDoubleArray>::New();

    if(this->SmoothingStatus == true)
    {
        SmoothStructuredPoints::SmoothGradientOfScalarField(
                this->SmoothingKernelSize,
                GridResolution,
                DeformationValue,
                DeformationValueGradient);
    }
    else
    {
        // this->ComputeDeformationValueGradient(
        //         InputDataSet,
        //         DeformationValue,
        //         DeformationValueGradient); // Output
        this->ComputeDeformationValueGradient2(
                InputDataSet,
                DeformationValueGradient); // Output
    }

    // Find Manifold Structures
    this->FindManifoldStructures(
            InputDataSet,
            DeformationValue,
            DeformationValueGradient,
            DeformationVector,
            OutputPolyData);   // Output

    // Add normals to the PolyData
    this->AddNormalsArray(OutputPolyData);

    // Compute and Add DeformationVector Deviation
    this->AddDeformationVectorDeviation(OutputPolyData);

    return 1;
}

// ================================
// Get Deformation Value And Vector
// ================================

// Description:
// Pointer to DeformationValue and DeformationVectors are declared outside the function. But
// they will be changed in this function to point to arrays in InputDataSet. So we need to
// "pass pointers by reference". Otherwise they will be null again outside this function.

void IdentifyStructures::GetDeformationValueAndVector(
        vtkDataSet *InputDataSet,
        vtkDoubleArray *& DeformationValue,    // Output
        vtkDoubleArray *& DeformationVector)   // Output
{
    // Get Point Data
    vtkPointData *InputPointData = InputDataSet->GetPointData();

    // Assign Deformation variables based on Structure mode to be identified
    switch(this->StructureMode)
    {
        // Max Strain
        case StructureType::MAX_STRAIN:
        {
            DeformationValue = vtkDoubleArray::SafeDownCast(InputPointData->GetArray("MaxStrainValues"));
            DeformationVector = vtkDoubleArray::SafeDownCast(InputPointData->GetArray("MaxStrainVectors"));
            break;
        }

        // Min Strain
        case StructureType::MIN_STRAIN:
        {
            DeformationValue = vtkDoubleArray::SafeDownCast(InputPointData->GetArray("MinStrainValues"));
            DeformationVector = vtkDoubleArray::SafeDownCast(InputPointData->GetArray("MinStrainVectors"));
            break;
        }

        // Undefined Structure mode
        default:
        {
            vtkErrorMacro("The Structure mode is undefined.");
        }
    }

    // Check DeformationValue Array
    if(DeformationValue == NULL)
    {
        vtkErrorMacro("Deformation Value array is NULL.");
    }
    else if(DeformationValue->GetNumberOfTuples() < 2)
    {
        vtkErrorMacro("DeformationValue array has no tuples.");
    }

    // Check DeformationVector Array
    if(DeformationVector == NULL)
    {
        vtkErrorMacro("Deformation Vector Array is NULL.");
    }
    else if(DeformationVector->GetNumberOfTuples() < 2)
    {
        vtkErrorMacro("Deformation Vector array has no tuples.");
    }
}

// ====================================
// Compute Deformation Value Gradient 2
// ====================================

void IdentifyStructures::ComputeDeformationValueGradient2(
        vtkDataSet *InputDataSet,
        vtkDoubleArray *DeformationValueGradient)   // Output
{
    // Declare Deformation Value array name
    char DeformationValueArrayName[256];

    // Get Deformation Value
    switch(this->StructureMode)
    {
        // Max Strain
        case StructureType::MAX_STRAIN:
        {
            strcpy(DeformationValueArrayName,"MaxStrainValues");
            break;
        }

        // Min Strain
        case StructureType::MIN_STRAIN:
        {
            strcpy(DeformationValueArrayName,"MinStrainValues");
            break;
        }

        // Shear
        case StructureType::SHEAR:
        {
            strcpy(DeformationValueArrayName,"ShearValues");
            break;
        }

        // Undefined Structure type
        default:
        {
            vtkErrorMacro("Structure type is undefined.");
        }
    }

    // use vtk's Gradient Filter
    vtkSmartPointer<vtkGradientFilter> GradientFilter = vtkSmartPointer<vtkGradientFilter>::New();
    #if VTK_MAJOR_VERSION <= 5
        GradientFilter->SetInput(InputDataSet);
    #else
        GradientFilter->SetInputData(InputDataSet);
    #endif
    GradientFilter->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,DeformationValueArrayName);
    GradientFilter->Update();
    
    // Get the output dataset of GradientFilter
    vtkSmartPointer<vtkDataSet> OutputDataSet = GradientFilter->GetOutput();
    
    // Set output of the function
    DeformationValueGradient->DeepCopy(vtkDoubleArray::SafeDownCast(OutputDataSet->GetPointData()->GetVectors("Gradients")));

    // Check Deformation Value Gradient
    if(DeformationValueGradient == NULL)
    {
        vtkErrorMacro("Deformation value gradient is NULL.");
    }
    else if(DeformationValueGradient->GetNumberOfTuples() < 1)
    {
        vtkErrorMacro("DeformationValueGradient does not have any tuples.");
    }
}

// ==================================
// Compute Deformation Value Gradient
// ==================================

void IdentifyStructures::ComputeDeformationValueGradient(
        vtkDataSet *InputDataSet,
        vtkDoubleArray *DeformationValue,
        vtkDoubleArray *DeformationValueGradient)   // Output
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

    // Check DeformationValue
    if(DeformationValue == NULL)
    {
        vtkErrorMacro("DeformationValue is NULL.");
    }
    else if(DeformationValue->GetNumberOfTuples() < 2 || 
            DeformationValue->GetNumberOfComponents() < 1)
    {
        vtkErrorMacro("DeformationValues does not have data.");
    }

    // Cast input to structured points
    vtkSmartPointer<vtkStructuredPoints> InputStructuredPoints = vtkStructuredPoints::SafeDownCast(InputDataSet);

    // Check casting
    if(InputStructuredPoints == NULL)
    {
        vtkErrorMacro("InputStructuredGrid is NULL.");
    }
    else if(InputStructuredPoints->GetNumberOfPoints() < 1)
    {
        vtkErrorMacro("InputStructuredGrid has no points.");
    }

    // Structured Points info
    unsigned int NumberOfPoints = InputDataSet->GetNumberOfPoints();
    int *GridResolution = InputStructuredPoints->GetDimensions();

    // Set the DeformationValueGradient array
    DeformationValueGradient->SetNumberOfComponents(DIMENSION);
    DeformationValueGradient->SetNumberOfTuples(NumberOfPoints);
    DeformationValueGradient->SetName("DeformationValueGradient");

    // Loop over points
    for(unsigned int PointId = 0; PointId < NumberOfPoints; PointId++)
    {
        // declare Gradient at point
        double Gradient[DIMENSION];

        // Loop over dimension
        for(unsigned int DimensionIterator = 0; DimensionIterator < DIMENSION; DimensionIterator++)
        {
            // Get Stencil Ids in the current dimension
            vtkIdType StencilIds[2];
            this->GetStencilIds(GridResolution,PointId,DimensionIterator,StencilIds);

            // Get DeformationValue on 1D stencil
            double BackValue = DeformationValue->GetTuple1(StencilIds[0]);
            double FrontValue = DeformationValue->GetTuple1(StencilIds[1]);

            // Get Positions on 1D stencil
            double BackPoint[DIMENSION];
            double FrontPoint[DIMENSION];
            InputDataSet->GetPoint(StencilIds[0],BackPoint);
            InputDataSet->GetPoint(StencilIds[1],FrontPoint);

            // Compute Gradient in specific direction
            Gradient[DimensionIterator] = (FrontValue - BackValue) / 
                (FrontPoint[DimensionIterator] - BackPoint[DimensionIterator]);
        }

        // Set gradient to output array
        DeformationValueGradient->SetTuple(PointId,Gradient);
    }
}

// ===============
// Get Point Index
// ===============

// Decsription:
// Converts the list style ID if a point into coordinate-wise indices if the point
// in a structured grid. Grid Resolution is number of points in each direction.
// Note: PointId is 0-offset. So the first point is 0. Therefore, the PointIndex are
// also offset from 0. So the first point is (0,0,0) in 3D.

void IdentifyStructures::GetPointIndex(
        int *GridResolution,
        vtkIdType PointId,
        unsigned int *PointIndex)  // Output
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

inline vtkIdType IdentifyStructures::GetPointId(
        int *GridResolution,
        vtkIdType *PointIndex)
{
    vtkIdType PointId = 0;
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

void IdentifyStructures::GetStencilIds(
        int *GridResolution,
        vtkIdType PointId,
        unsigned int TargetDirection,
        vtkIdType *StencilIds)   // Output
{
    // Get point index
    unsigned int PointIndex[DIMENSION];
    this->GetPointIndex(GridResolution,PointId,PointIndex);

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
            if(PointIndex[TargetDirection] < static_cast<unsigned int>(GridResolution[TargetDirection]-1))
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

// ======================
// Normalize Vector Array
// ======================

// Description:
// This filter over-writes the anchoral array.
// This operates on tuple vectors of vtkDoubleArray.

void IdentifyStructures::NormalizeVectorArray(vtkDoubleArray *VectorArray)
{
    // Check Input
    if(VectorArray == NULL)
    {
        vtkErrorMacro("VectorArray is NULL.");
        exit(0);
    }
    else if(VectorArray->GetNumberOfTuples() < 1)
    {
        vtkErrorMacro("VectroArray does not have any tuples.");
        exit(0);
    }

    // Number of Tuples
    unsigned int NumberOfTuples = VectorArray->GetNumberOfTuples();

    // Loop over tuples
    for(unsigned int TupleIterator = 0; TupleIterator < NumberOfTuples; TupleIterator++)
    {
        // Get a pointer to the tuple
        double *Tuple = VectorArray->GetTuple(TupleIterator);

        // Normalize tuple
        vtkMath::Normalize(Tuple);

        // Set normalized tuples again to the array (this is necessary)
        VectorArray->SetTuple(TupleIterator,Tuple);
    }
}

// ========================
// Find Manifold Structures
// ========================

// Description:
// This is the main method of this class.
// Note: This method works on Structured Points dataset.

void IdentifyStructures::FindManifoldStructures(
        vtkDataSet *InputDataSet,
        vtkDoubleArray *DeformationValue,
        vtkDoubleArray *DeformationValueGradient,
        vtkDoubleArray *DeformationVector,
        vtkPolyData *OutputPolyData)
{
    // Check Input DataSet
    if(InputDataSet == NULL)
    {
        vtkErrorMacro("Input DataSet is NULL.");
    }
    else if(InputDataSet->GetNumberOfCells() < 1)
    {
        vtkErrorMacro("Input DataSet does not have cell.");
    }

    // Declare outputs of this function
    PointsListType IsoSurfacePointsList;
    TrianglesListType IsoSurfaceTrianglesList;

    // Get Grid Resolution
    vtkSmartPointer<vtkStructuredPoints> InputStructuredPoints = vtkStructuredPoints::SafeDownCast(InputDataSet);
    int *GridResolution = InputStructuredPoints->GetDimensions();

    // Define an iterator object for cells
    // vtkSmartPointer<vtkCellIterator> CellIterator = InputDataSet->NewCellIterator();
    unsigned int NumberOfCells = InputDataSet->GetNumberOfCells();   // CHANGED
    unsigned int ProgressCounter = 0;

    // Iterate over cubes
    // for(CellIterator->InitTraversal(); CellIterator->IsDoneWithTraversal(); CellIterator->GoToNextCell())
    #ifdef _OPENMP
    #pragma omp parallel for \
    default (none) \
    shared(NumberOfCells,InputDataSet,GridResolution,DeformationValue,DeformationValueGradient, \
            DeformationVector,IsoSurfacePointsList,IsoSurfaceTrianglesList,ProgressCounter)
    #endif
    for(unsigned int CellIterator = 0; CellIterator < NumberOfCells; CellIterator++)  // CHANGED
    {
        // Check cell type to be a Voxel
        #ifdef _OPENMP
        vtkSmartPointer<vtkGenericCell> Cell = vtkSmartPointer<vtkGenericCell>::New();
        InputDataSet->GetCell(CellIterator,Cell);
        #else
        vtkSmartPointer<vtkCell> Cell = InputDataSet->GetCell(CellIterator);
        #endif

        // if(CellIterator->GetCellType() != VTK_VOXEL)
        if(Cell->GetCellType() != VTK_VOXEL)  // CHANGED
        {
            vtkErrorMacro("Cell type is not supported for Marching Cubes algorithm.");
        }

        // Get a cube (as 8 Ids)
        // vtkSmartPointer<vtkIdList> CubeIds = CellIterator->GetPointIds();
        vtkSmartPointer<vtkIdList> CubeIds = Cell->GetPointIds();  // CHANGED

        // Find Anchor vertex of Cube
        int AnchorVertexId = FindMinimumOfIds(CubeIds);

        // Find Cube Ids in order
        vtkSmartPointer<vtkIdList> OrderedCubeIds = vtkSmartPointer<vtkIdList>::New();
        this->FindOrderedCubePointIds(GridResolution,AnchorVertexId,OrderedCubeIds);

        // Get Anchor Index
        unsigned int AnchorVertexIndex[DIMENSION];
        this->GetPointIndex(GridResolution,AnchorVertexId,AnchorVertexIndex);

        // Get Points Coordinates on the cube
        double PointsCoordinatesOnCube[CUBESIZE][DIMENSION];
        this->GetPointsCoordinatesOnCube(
                InputDataSet,
                OrderedCubeIds,
                PointsCoordinatesOnCube);   // Output

        // Get Deformation values on cube vertices
        double DeformationValuesOnCube[CUBESIZE];
        this->GetScalarDataOnCube(
                DeformationValue,
                OrderedCubeIds,
                DeformationValuesOnCube);

        // Get Deformation Value Gradients on the Cube vertices
        double DeformationValueGradientsOnCube[CUBESIZE][DIMENSION];
        this->GetVectorDataOnCube(
                DeformationValueGradient,
                OrderedCubeIds,
                DeformationValueGradientsOnCube);  // Output

        // Get Deformation Vectors on the Cube vertices
        double DeformationVectorsOnCube[CUBESIZE][DIMENSION];
        // double ** DeformationVectorsOnCube;   // CHANGED
        this->GetVectorDataOnCube(
                DeformationVector,
                OrderedCubeIds,
                DeformationVectorsOnCube);  // Output

        // Define a Cube object
        CubeCell Cube;
        Cube.SetVertexPointsCoordinates(PointsCoordinatesOnCube);
        Cube.SetVertexPointsOriginalIds(OrderedCubeIds);
        Cube.SetDeformationValues(DeformationValue);
        Cube.SetDeformationValuesOnCube(DeformationValuesOnCube);
        Cube.SetDeformationValueGradientsOnCube(DeformationValueGradientsOnCube);
        Cube.SetDeformationVectorsOnCube(DeformationVectorsOnCube);
        Cube.SetIsoSurfaceLevel(0.0);
        Cube.SetGridResolution(GridResolution,DIMENSION);
        Cube.SetAnchorVertexIndex(AnchorVertexIndex,DIMENSION); 

        // Apply the Cube algorithm
        bool IsoSurfaceStatus = Cube.FindIsoSurface();

        // Get output of Cube algorithm
        if(IsoSurfaceStatus == true)
        {
            PointsListType IntersectionPointsListOnCube = Cube.GetIntersectionPointsList();
            TrianglesListType IntersectionTrianglesListOnCube = Cube.GetIntersectionTrianglesList();

            // Store Intersection point
            this->StoreIntersectionPoints(
                    GridResolution,
                    AnchorVertexId,
                    IntersectionPointsListOnCube,
                    IsoSurfacePointsList);  // Output

            // Store Intersection Triangles
            this->StoreIntersectionTriangles(
                    GridResolution,
                    AnchorVertexId,
                    IntersectionTrianglesListOnCube,
                    IsoSurfaceTrianglesList);   // Output
        }

        // Clear Cell
        Cube.ClearCell();

        // Update Progress
        #ifdef _OPENMP
        #pragma omp atomic
        #endif
        ProgressCounter++;
        this->FilterUpdateProgress(ProgressCounter,NumberOfCells);
    }

    // Create Output PolyData
    this->CreateOutputPolyData(
            IsoSurfacePointsList,
            IsoSurfaceTrianglesList,
            OutputPolyData);

    // Delete Cell Iterator
    // CellIterator->Delete();  // CHANGED
}

// ===================
// Find Minimum Of Ids
// ===================

// Description:
// Find minim Id among list of Ids. Ids are vertex of cube Ids. The minimum id indicates
// the anchor (or bae vertex) of cube, which is used for re-ordering of vertices.

inline vtkIdType IdentifyStructures::FindMinimumOfIds(vtkIdList *CubeIds)
{
    vtkIdType MinimumId = CubeIds->GetId(0);

    // Loop over Ids
    for(unsigned int IdIterator = 1;
        IdIterator < static_cast<unsigned int>(CubeIds->GetNumberOfIds());
        IdIterator++)
    {
        if(MinimumId > CubeIds->GetId(IdIterator))
        {
            MinimumId = CubeIds->GetId(IdIterator);
        }
    }

    return MinimumId;
}

// ===========================
// Find Ordered Cube Point Ids
// ===========================

// Description:
// Adds vertices of cube in a special order.
// Starting from anchor vertex as vertex 0, a counter clockwise counting of verices on buttom face
// of cube adds points 1,2 and 3. The same for points in the upper face of cube adds vertices 4,5,6 
// and 7.

void IdentifyStructures::FindOrderedCubePointIds(
        int *GridResolution,
        vtkIdType AnchorVertexId,
        vtkIdList *OrderedCubeIds)  // Output
{
}

// ==============================
// Get Points Coordinates On Cube
// ==============================

// Note: This method only performs shallow copy. It does not allocate memory to copy point
// coordinates to output. It only copeis the pointers to point coordinates and set them into
// a 2D array.

void IdentifyStructures::GetPointsCoordinatesOnCube(
        vtkDataSet *InputDataSet,
        vtkIdList *OrderedCubeIds,
        double PointsCoordinatesOnCube[CUBESIZE][DIMENSION])   // Output
{
    // Loop over cube points
    for(unsigned int PointIterator = 0; PointIterator < CUBESIZE; PointIterator++)
    {
        #ifdef _OPENMP
        double APointOnCube[DIMENSION];
        InputDataSet->GetPoint(OrderedCubeIds->GetId(PointIterator),APointOnCube);
        #else
        double *APointOnCube = InputDataSet->GetPoint(OrderedCubeIds->GetId(PointIterator));
        #endif

        // Deep Copy
        for(unsigned int DimensionIterator = 0; DimensionIterator < DIMENSION; DimensionIterator++)
        {
            PointsCoordinatesOnCube[PointIterator][DimensionIterator] = APointOnCube[DimensionIterator];
        }
    }
}

// =======================
// Get Scalar Data On Cube
// =======================

void IdentifyStructures::GetScalarDataOnCube(
        vtkDoubleArray *InputDoubleArray,
        vtkIdList *OrderedCubeIds,
        double ScalarDataOnCube[CUBESIZE])
{
    // Loop over cube points
    for(unsigned int PointIterator = 0; PointIterator < CUBESIZE; PointIterator++)
    {
        #ifdef _OPENMP
        double TempScalar[1];
        InputDoubleArray->GetTuple(OrderedCubeIds->GetId(PointIterator),TempScalar);
        ScalarDataOnCube[PointIterator] = TempScalar[0];
        #else
        ScalarDataOnCube[PointIterator] = InputDoubleArray->GetTuple1(OrderedCubeIds->GetId(PointIterator));
        #endif
    } 
}

// =======================
// Get Vector Data On Cube
// =======================

// Gets vtkDoubleArray as input. It can be a Point Data of the actual data set.
// Then it gets data from the Ids List that is provided in OrderedCubeIds.
// The output is an array of Vectors class. Each Vector holds a tuple from DoubleArray.
// Note: This perfomrs a shallow copy; just puts the pointers of actual data tuples
// into a 2D array, but does not allocate any memory.

void IdentifyStructures::GetVectorDataOnCube(
        vtkDoubleArray *InputDoubleArray,
        vtkIdList *OrderedCubeIds,
        double VectorDataOnCube[CUBESIZE][DIMENSION])   // Output
{
    // Loop over cube points
    for(unsigned int PointIterator = 0; PointIterator < CUBESIZE; PointIterator++)
    {
        #ifdef _OPENMP
        double TempTuple[DIMENSION];
        InputDoubleArray->GetTuple(OrderedCubeIds->GetId(PointIterator),TempTuple);
        for(unsigned int DimensionIterator = 0; DimensionIterator < DIMENSION; DimensionIterator++)
        {
            VectorDataOnCube[PointIterator][DimensionIterator] = TempTuple[DimensionIterator];
        }
        #else
        for(unsigned int DimensionIterator = 0; DimensionIterator < DIMENSION; DimensionIterator++)
        {
                InputDoubleArray->GetComponent(OrderedCubeIds->GetId(PointIterator),DimensionIterator);
        }
        #endif
    }
}

// ==================
// Get Global Edge Id
// ==================

// Description:
// Gives a unique Id for an edge of any cube in the whole Structured point data set. Although cubes 
// share edges, but a unique Id can be given for an edge, irregardless of which cube(s) it belongs to.
// Input to this method is "Local" edge Id on a cube (from 0 to 11). 

unsigned int IdentifyStructures::GetGlobalEdgeId(
        int *GridResolution,
        vtkIdType AnchorVertexId,
        vtkIdType LocalEdgeId)   // Output
{
    // Get Index of Anchor vertex of the cube w.r.t integer coordinates of structured point
    unsigned int AnchorVertexIndex[DIMENSION];
    this->GetPointIndex(GridResolution,AnchorVertexId,AnchorVertexIndex);

    // Use simpler variables for index of Anchor vertex of cube
    vtkIdType Idx = AnchorVertexIndex[0];
    vtkIdType Idy = AnchorVertexIndex[1];
    vtkIdType Idz = AnchorVertexIndex[2];

    // Output of function
    vtkIdType GlobalEdgeId;

    // Edge Reference Vertex is a point that edge gets its global id from.
    // Edge Global Id = 3 * reference Index Id + 0, 1 or 2 depending on direction.
    vtkIdType EdgeReferenceVertexId;

    return 1;
}

// =========================
// Store Intersection Points
// =========================

// Description:
// Note: If we use the NON-efficient method for Finding Intersection Points
// in the CubeCell class, we need to check and remove the repetitive points 
// in this method. Curently it has been enabled.

void IdentifyStructures::StoreIntersectionPoints(
        int *GridResolution,
        unsigned int AnchorVertexId,
        PointsListType IntersectionPointsListOnCube,
        PointsListType & IsoSurfacePointsList)   // Output
{
    // Check if no point is available
    if(IntersectionPointsListOnCube.size() > 0)
    {
        // Loop over intersection points
        PointsListType::iterator PointsListOnCubeIterator;
        for(PointsListOnCubeIterator = IntersectionPointsListOnCube.begin();
            PointsListOnCubeIterator != IntersectionPointsListOnCube.end();
            PointsListOnCubeIterator++)
        {
            vtkIdType LocalEdgeId = PointsListOnCubeIterator->first;
           
            // Convert Local Edge Id to global Edge Id
            unsigned int GlobalEdgeId = this->GetGlobalEdgeId(
                    GridResolution,
                    AnchorVertexId,
                    LocalEdgeId);

            // Check if the GlobalEdgeId already included in the list
            PointsListType::iterator PointsListIterator;
            PointsListIterator = IsoSurfacePointsList.find(GlobalEdgeId);
            
            // Add the new point only if it's not being included before
            if(PointsListIterator == IsoSurfacePointsList.end())
            {
                // Define Point and Id Pair object
                PointDataType PointData(PointsListOnCubeIterator->second);
                
                // Insert GlobalEdgeId and PointData object to the output map
                #ifdef _OPENMP
                #pragma omp critical
                #endif
                IsoSurfacePointsList.insert(std::make_pair(GlobalEdgeId,PointData));
            }
        }
    } 
}

// ============================
// Store Intersection Triangles
// ============================

// Description:
// Converts local edge Ids of triangle points to the global edge Ids.
// Then stores all triangles in a large list of IsoSurfaceTrianglesList.

void IdentifyStructures::StoreIntersectionTriangles(
        int *GridResolution,
        unsigned int AnchorVertexId,
        TrianglesListType IntersectionTrianglesListOnCube,
        TrianglesListType & IsoSurfaceTrianglesList)   // Output
{
    // Check if no triangle is available
    if(IntersectionTrianglesListOnCube.size() > 0)
    {
        // Loop over intersection triangles
        TrianglesListType::iterator TriangleIterator;
        for(TriangleIterator = IntersectionTrianglesListOnCube.begin();
            TriangleIterator != IntersectionTrianglesListOnCube.end();
            TriangleIterator++)
        {
            // Iterator over Triangle points
            TriangleType::iterator TrianglePointIterator;
            for(TrianglePointIterator = TriangleIterator->begin();
                TrianglePointIterator != TriangleIterator->end();
                TrianglePointIterator++)
            {
                // Get the Local Edge Id that triangle point resides on
                vtkIdType LocalEdgeId = *TrianglePointIterator;

                // Convert local edge Id to Global edge Id
                vtkIdType GlobalEdgeId = this->GetGlobalEdgeId(
                        GridResolution,
                        AnchorVertexId,
                        LocalEdgeId);

                // Overwrite the Triangle with GlobalEdgeId
                *TrianglePointIterator = GlobalEdgeId;
            }

            // Insert Triangle to the output list
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            IsoSurfaceTrianglesList.push_back(*TriangleIterator);
        }
    }
}

// ======================
// Create Output PolyData
// ======================

void IdentifyStructures::CreateOutputPolyData(
        PointsListType & IsoSurfacePointsList,
        TrianglesListType & IsoSurfaceTrianglesList,
        vtkPolyData *OutputPolyData)
{
    if(this->DebugStatus == true)
    {
        std::cout << "Number of IsoSurface Points: "  << IsoSurfacePointsList.size() << std::endl;
        std::cout << "Number of IsoSurface Triangles: " << IsoSurfaceTrianglesList.size() << std::endl;
        std::cout << "Rearrange points ..." << std::endl;
    }

    // Organize Point Ids
    this->OrganizeIsoSurfacePointIds(IsoSurfacePointsList);

    // Organize Triangle Point Ids
    this->OrganizeIsoSurfaceTrianglePointIds(IsoSurfacePointsList,IsoSurfaceTrianglesList);

    // Create PolyData //
    
    // 1- Points and PointData

    // Points
    vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();
    Points->SetNumberOfPoints(IsoSurfacePointsList.size());

    // Deformation Values
    vtkSmartPointer<vtkDoubleArray> DeformationValues = vtkSmartPointer<vtkDoubleArray>::New();
    DeformationValues->SetNumberOfComponents(1);
    DeformationValues->SetNumberOfTuples(IsoSurfacePointsList.size());
    DeformationValues->SetName("DeformationValues");

    // Deformation Vectors
    vtkSmartPointer<vtkDoubleArray> DeformationVectors = vtkSmartPointer<vtkDoubleArray>::New();
    DeformationVectors->SetNumberOfComponents(DIMENSION);
    DeformationVectors->SetNumberOfTuples(IsoSurfacePointsList.size());
    DeformationVectors->SetName("DeformationVectors");

    // Deformation Value Curvature
    vtkSmartPointer<vtkDoubleArray> DeformationValueCurvatures = vtkSmartPointer<vtkDoubleArray>::New();
    DeformationValueCurvatures->SetNumberOfComponents(1);
    DeformationValueCurvatures->SetNumberOfTuples(IsoSurfacePointsList.size());
    DeformationValueCurvatures->SetName("DeformationValueCurvatures");

    // Iterate over points
    PointsListType::iterator PointsListIterator;
    for(PointsListIterator = IsoSurfacePointsList.begin();
        PointsListIterator != IsoSurfacePointsList.end();
        PointsListIterator++)
    {
        // Get PointData object
        PointDataType PointData = PointsListIterator->second;

        // Get Organized Id of point
        unsigned int OrganizedPointId = PointData.GetOrganizedPointId();

        // Insert coordinate of point
        Vector PointCoordinatesAsVector = PointData.GetCoordinates();
        double *PointCoordinates = PointCoordinatesAsVector.GetComponents();
        Points->SetPoint(OrganizedPointId,PointCoordinates);

        // Insert Deformation Value
        double PointDeformationValue = PointData.GetDeformationValue();
        DeformationValues->SetTuple1(OrganizedPointId,PointDeformationValue);

        // Insert Deformation Vector
        Vector PointDeformationVectorAsVector = PointData.GetDeformationVector();
        for(unsigned int DimensionIterator = 0; DimensionIterator < DIMENSION; DimensionIterator++)
        {
            DeformationVectors->SetComponent(
                OrganizedPointId,
                DimensionIterator,
                PointDeformationVectorAsVector[DimensionIterator]);
        }

        // Insert Deformation Value Curvature
        double PointDeformationValueCurvature = PointData.GetDeformationValueCurvature();
        DeformationValueCurvatures->SetTuple1(OrganizedPointId,PointDeformationValueCurvature);
    }

    // 2- Triangles

    // Create list of all triangle cells
    vtkSmartPointer<vtkCellArray> TriangleCells = vtkSmartPointer<vtkCellArray>::New();

    // Iterator over each triangle
    TrianglesListType::iterator TriangleIterator;
    for(TriangleIterator = IsoSurfaceTrianglesList.begin();
        TriangleIterator != IsoSurfaceTrianglesList.end();
        TriangleIterator++)
    {
        // Create a single Triangle Cell
        vtkSmartPointer<vtkTriangle> TriangleCell = vtkSmartPointer<vtkTriangle>::New();

        // Iterator over Triangle Points
        TriangleType::iterator TrianglePointIterator;
        unsigned int PointCounter = 0;
        for(TrianglePointIterator = TriangleIterator->begin();
            TrianglePointIterator != TriangleIterator->end();
            TrianglePointIterator++)
        {
            // Add Point Id t0 the triangle cell
            TriangleCell->GetPointIds()->SetId(PointCounter,*TrianglePointIterator);
            PointCounter++;
        }
        
        // Add Triangle cell to the list of Triangle Cells
        TriangleCells->InsertNextCell(TriangleCell);     
    }

    // Add Data to PolyData
    OutputPolyData->SetPoints(Points);
    OutputPolyData->SetPolys(TriangleCells);
    OutputPolyData->GetPointData()->AddArray(DeformationValues);
    OutputPolyData->GetPointData()->AddArray(DeformationVectors);
    OutputPolyData->GetPointData()->AddArray(DeformationValueCurvatures);

    if(this->DebugStatus == true)
    {
        std::cout << "Done!" << std::endl;
    }
}

// ==============================
// Organize Iso-Surface Point Ids
// ==============================

// Description:
// Gives an ordered Id for each point starting from 1.

void IdentifyStructures::OrganizeIsoSurfacePointIds(
        PointsListType & IsoSurfacePointsList)
{
    // Id Counter starting from zero
    unsigned int IdCounter = 0;

    // Loop over Points
    PointsListType::iterator PointsListIterator;
    for(PointsListIterator = IsoSurfacePointsList.begin();
        PointsListIterator != IsoSurfacePointsList.end();
        PointsListIterator++)
    {
        // Set Organized Id to the counter Id
        PointsListIterator->second.SetOrganizedPointId(IdCounter);

        // Update Counter Id
        IdCounter++;
    }
}

// =======================================
// Organize Iso-Surface Triangle Point Ids
// =======================================

void IdentifyStructures::OrganizeIsoSurfaceTrianglePointIds(
        PointsListType & IsoSurfacePointsList,
        TrianglesListType & IsoSurfaceTrianglesList)
{
    // Iterate over triangles
    TrianglesListType::iterator TriangleIterator;

    for(TriangleIterator = IsoSurfaceTrianglesList.begin();
        TriangleIterator != IsoSurfaceTrianglesList.end();
        TriangleIterator++)
    {
        // Iterate over points of each triangle
        TriangleType::iterator TrianglePointIterator;
        for(TrianglePointIterator = TriangleIterator->begin();
            TrianglePointIterator != TriangleIterator->end();
            TrianglePointIterator++)
        {
            // Get GlobalEdgeId of triangle point
            unsigned int GlobalEdgeId = *TrianglePointIterator;

            // Find the The PointData object from GlobalEdgeId
            PointsListType::iterator PointsListIterator;
            PointsListIterator = IsoSurfacePointsList.find(GlobalEdgeId);

            // Check if the key found
            if(PointsListIterator == IsoSurfacePointsList.end())
            {
                vtkErrorMacro("No GlobalEdgeId found in the Points List.");

                if(this->DebugStatus == true)
                {
                    int ProblemAnchorVertexId = floor(static_cast<double>(GlobalEdgeId / 3));
                    int ProblemLocalEdgeId = GlobalEdgeId - 3*ProblemAnchorVertexId;
                    unsigned int ProblemAnchorVertexIndex[3];
                    int GridResolution[3] = {101,101,101};       // BIG TODO Get GridResolution from arguments!!!
                    this->GetPointIndex(GridResolution,ProblemAnchorVertexId,ProblemAnchorVertexIndex);

                    std::cout << "GlobalEdgeID: "   << GlobalEdgeId << std::endl;
                    std::cout << "LocalEdgeId: "    << ProblemLocalEdgeId << std::endl;
                    std::cout << "AnchorVertexId: " << ProblemAnchorVertexId << std::endl;
                    std::cout << "AnchorVertexIndex: ";
                    std::cout << ProblemAnchorVertexIndex[0] << ", " ;
                    std::cout << ProblemAnchorVertexIndex[1] << ", ";
                    std::cout << ProblemAnchorVertexIndex[2] << std::endl;
                }
            }

            // Get the PointAndId object
            PointDataType PointData;
            PointData = PointsListIterator->second;

            // Set triangle point to Organized PointId
            *TrianglePointIterator = PointData.GetOrganizedPointId();
        }
    }
}

// =================
// Add Normals Array
// =================

// Description:
// Computes an adds normals to the output polydata. Moreover, it re-orientates triangles
// consistently so that all normals are lie in consistent directions. Finally it replaces
// the original triangles on output polydata with the re-oriented triangles.

void IdentifyStructures::AddNormalsArray(vtkPolyData *OutputPolyData)
{
    if(this->DebugStatus == true)
    {
        std::cout << "Compute normals ..." << std::endl;
    }

    if(OutputPolyData->GetNumberOfPoints() > 0)
    {
        vtkSmartPointer<vtkPolyDataNormals> PolyDataNormalsFilter = vtkSmartPointer<vtkPolyDataNormals>::New();
        #if VTK_MAJOR_VERSION <= 5
            PolyDataNormalsFilter->SetInput(OutputPolyData);
        #else
            PolyDataNormalsFilter->SetInputData(OutputPolyData);
        #endif
        PolyDataNormalsFilter->ComputePointNormalsOn();
        PolyDataNormalsFilter->ComputeCellNormalsOff();
        PolyDataNormalsFilter->ConsistencyOn();
        PolyDataNormalsFilter->SplittingOff();
        PolyDataNormalsFilter->SetFeatureAngle(10);
        PolyDataNormalsFilter->AutoOrientNormalsOn();
        PolyDataNormalsFilter->FlipNormalsOff();
        PolyDataNormalsFilter->NonManifoldTraversalOn();
        PolyDataNormalsFilter->Update();

        // Get output of the filter
        vtkSmartPointer<vtkPolyData> PolyDataWithNormals;
        PolyDataWithNormals = PolyDataNormalsFilter->GetOutput();

        // Check output
        if(PolyDataWithNormals == NULL)
        {
            vtkErrorMacro("PolyDataWithNormals is NULL.");
        }
        else if(PolyDataWithNormals->GetNumberOfPoints() < 1)
        {
            vtkErrorMacro("PolyDataWithNormals is empty.");
        }

        // Get Normals array
        vtkSmartPointer<vtkDataArray> Normals;
        Normals = PolyDataWithNormals->GetPointData()->GetArray("Normals");

        // Check Normals
        if(Normals == NULL)
        {
            vtkErrorMacro("Normals is NULL.");
        }
        else if(Normals->GetNumberOfTuples() < 1)
        {
            vtkErrorMacro("Normals has no tuples.");
        }

        // Add normals array to the OutputPolyData
        OutputPolyData->GetPointData()->SetNormals(Normals);

        // Reset Cells with Oriented-Consistent cells
        OutputPolyData->SetPolys(PolyDataWithNormals->GetPolys());
    }
}

// ================================
// Add Deformation Vector Deviation
// ================================

// Description:
// The Deviation angle is the angle between the Iso-surface normal and the Deformation
// vector for each point. It's favorable that these two vectors align. However, for
// non-integrable manifolds the angle is not necessarily zero.
// This function adds a new array to the output PolyData called DeviationAngle.

void IdentifyStructures::AddDeformationVectorDeviation(vtkPolyData *OutputPolyData)
{
    unsigned int NumberOfPoints = OutputPolyData->GetNumberOfPoints();

    // Get Normals
    vtkSmartPointer<vtkDataArray> Normals = OutputPolyData->GetPointData()->GetArray("Normals");

    // Check Normals
    if(Normals == NULL)
    {
        vtkErrorMacro("No Normals found in PolyData.");
    }
    else if(Normals->GetNumberOfTuples() < 1)
    {
        vtkErrorMacro("Normals array is empty.");
    }

    // Get Deformation Vectors
    vtkSmartPointer<vtkDoubleArray> DeformationVectors = vtkDoubleArray::SafeDownCast(
            OutputPolyData->GetPointData()->GetArray("DeformationVectors"));

    // Check Deformation Vectors
    if(DeformationVectors == NULL)
    {
        vtkErrorMacro("No DeformationVectors found in PolyData.");
    }
    else if(DeformationVectors->GetNumberOfTuples() < 1)
    {
        vtkErrorMacro("DeformationVectors is empty.");
    }

    // Define Deformation Vectors Deviation array
    vtkSmartPointer<vtkDoubleArray> DeformationVectorsDeviation = vtkSmartPointer<vtkDoubleArray>::New();
    DeformationVectorsDeviation->SetNumberOfTuples(1);
    DeformationVectorsDeviation->SetNumberOfTuples(NumberOfPoints);
    DeformationVectorsDeviation->SetName("DeformationVectorsDeviation");

    // Iterate over points
    for(unsigned int PointId = 0; PointId < NumberOfPoints; PointId++)
    {
        // Get Deformation Vector
        double *DeformationVector = DeformationVectors->GetTuple(PointId);

        // Get Normal
        double *Normal = Normals->GetTuple(PointId);

        // Deviation Angle
        double DeviationAngle = fabs(vtkMath::Dot(DeformationVector,Normal));

        // Add Deviation angle to the array
        DeformationVectorsDeviation->SetTuple1(PointId,DeviationAngle);
    }

    // Add Deformayion Vectors Deviation array to the PolyData
    OutputPolyData->GetPointData()->AddArray(DeformationVectorsDeviation);
}
