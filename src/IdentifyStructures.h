/*
 * =====================================================================================
 *
 *       Filename:  IdentifyStructures.h
 *
 *    Description:  Finds Coherent Structures
 *
 *        Version:  1.0
 *        Created:  02/14/2014 12:20:10 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University Of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __IdentifyStructures_h
#define __IdentifyStructures_h

// ======
// Macros
// ======
//
// Disable MSDN security features (WINDOWS only)
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

// Debugging
#ifndef HERE
#define HERE std::cout << "\033[0;44m" << "DEBUG: " << __FILE__ << " at " << __LINE__ << "\033[0m" << std::endl;
#endif

// Dimension
#ifndef DIMENSION
#define DIMENSION 3
#endif

#ifndef CUBESIZE
#define CUBESIZE 8
#endif

// =====================
// Foreward Declarations
// =====================

// Complete declarations
#include <vtkDataSetAlgorithm.h>
#include <vtkSmartPointer.h>
#include "Vector.h"
#include <vector>
// #include <array>  // CHANGED
#include <map>

// Incomplete declarations
class vtkDoubleArray;
class vtkIdList;
class PointDataType;

// =====
// Types
// =====

// Structure Type
class StructureType
{
    public:
        enum StructureModeType
        {
            MAX_STRAIN = 0,
            MIN_STRAIN,
            SHEAR,
            NumberOfStructureTypes 
        };
};

// Type for Intersection Points List
// The format is: <LocalEdgeId,Point>
// (avoid conflict of defintion in CubeCell.h)
#ifndef __POINTS_LIST_TYPE
#define __POINTS_LIST_TYPE
typedef std::map<unsigned int,PointDataType> PointsListType;
#endif

// Triangle Type
#ifndef __TRIANGLE_TYPE
#define __TRIANGLE_TYPE
typedef std::vector<unsigned int> TriangleType;
#endif
// typedef std::array<unsigned int,3> TriangleType; // CHANGED

// Type for Iso Surface Triangles List
// (avoid conflict of defintion in CubeCell.h)
#ifndef __TRIANGLES_LIST_TYPE
#define __TRIANGLES_LIST_TYPE
typedef std::vector<TriangleType> TrianglesListType;
#endif

// =========================
// Identify Structures Class
// =========================

class IdentifyStructures : public vtkDataSetAlgorithm
{
    public:
        static IdentifyStructures *New();
        vtkTypeRevisionMacro(IdentifyStructures,vtkDataSetAlgorithm);
        virtual void PrintSelf(ostream &os,vtkIndent Indent);

        // Accessros / Mutators
        vtkGetMacro(StructureMode,StructureType::StructureModeType);
        void SetStructureMode(int InputStructureMode);
        void SetStructureToMaxStrain() { this->StructureMode = StructureType::MAX_STRAIN; };
        void SetStructureToMinStrain() { this->StructureMode = StructureType::MIN_STRAIN; };
        void SetStructureToShear() { this->StructureMode = StructureType::SHEAR; };

        vtkGetMacro(SmoothingStatus,bool);
        vtkSetMacro(SmoothingStatus,bool);
        void SmoothingOn() { this->SmoothingStatus = true; };
        void SmoothingOff() { this->SmoothingStatus = false; };

        vtkGetMacro(SmoothingKernelSize,unsigned int);
        vtkSetMacro(SmoothingKernelSize,unsigned int);

        vtkGetMacro(DebugStatus,bool);
        vtkSetMacro(DebugStatus,bool);
        void DebugOn();
        void DebugOff();

        vtkGetMacro(ProgressStatus,bool);
        vtkSetMacro(ProgressStatus,bool);
        void ProgressOn() { this->ProgressStatus = true; };
        void ProgressOff() { this->ProgressStatus = false; };

    protected:
        IdentifyStructures();
        virtual ~IdentifyStructures();

        // Pipeline Executives
        virtual int FillOutputPortInformation(
                int port,
                vtkInformation *info);

        virtual int RequestDataObject(
                vtkInformation *request,
                vtkInformationVector **InputInfo,
                vtkInformationVector *OutputInfo);

        virtual int RequestData(
                vtkInformation *Request,
                vtkInformationVector **InputVector,
                vtkInformationVector *OutputVector);

        static void ProgressFunction(
                vtkObject *Caller,
                unsigned long int EventId,
                void *ClientData,
                void *CallData);
        
        void FilterUpdateProgress(
                unsigned int Step,
                unsigned int NumberOfSteps);

        // Methods
        void GetDeformationValueAndVector(
                vtkDataSet *InputDataSet,
                vtkDoubleArray *& DeformationValue,     // Output
                vtkDoubleArray *& DeformationVector);   // Output

        void ComputeDeformationValueGradient(
                vtkDataSet *InputDataSet,
                vtkDoubleArray *DeformationValue,
                vtkDoubleArray *DeformationValueGradient);  // Output

        void ComputeDeformationValueGradient2(
                vtkDataSet *InputDataSet,
                vtkDoubleArray *DeformationValueGradient);  // Output

        void NormalizeVectorArray(vtkDoubleArray *VectorArray);

        inline void GetPointIndex(
                int *GridResolution,
                vtkIdType PointId,
                unsigned int *PointIndex);  // Output

        inline vtkIdType GetPointId(
                int *GridResolution,
                vtkIdType *PointIndex);

        inline void GetStencilIds(
                int *GridResolution,
                vtkIdType PointId,
                unsigned int TargetDirection,
                vtkIdType *StencilIds);  // Output

        void FindManifoldStructures(
                vtkDataSet *InputDataSet,
                vtkDoubleArray *DeformationValue,
                vtkDoubleArray *DeformationValueGradient,
                vtkDoubleArray *DeformationVector,
                vtkPolyData *OutputPolyData);   // Output

        inline vtkIdType FindMinimumOfIds(vtkIdList *CubeIds);

        void FindOrderedCubePointIds(
                int *GridResolution,
                vtkIdType AnchorVertexId,
                vtkIdList *OrderedCubeIds);   // Output

        void GetPointsCoordinatesOnCube(
                vtkDataSet *InputDataSet,
                vtkIdList *OrderedCubeIds,
                double PointsCoordinatesOnCube[CUBESIZE][DIMENSION]);   // Output
        
        void GetScalarDataOnCube(
                vtkDoubleArray *InputDoubleArray,
                vtkIdList *OrderedCubeIds,
                double ScalarDataOnCubei[CUBESIZE]);    // Output

        void GetVectorDataOnCube(
                vtkDoubleArray *InputDoubleArray,
                vtkIdList *OrderedCubeIds,
                double VectorDataOnCube[CUBESIZE][DIMENSION]);   // Output

        unsigned int GetGlobalEdgeId(
                int *GridResolution,
                vtkIdType AnchorVertexId,
                vtkIdType LocalEdgeId);   // Output

        void StoreIntersectionPoints(
                int *GridResolution,
                unsigned int AnchorVertexId,
                PointsListType IntersectionPointsListOnCube,
                PointsListType & IsoSurfacePointsList);   // Output

        void StoreIntersectionTriangles(
                int *GridResolution,
                unsigned int AnchorVertexId,
                TrianglesListType IntersectionTrianglesListOnCube,
                TrianglesListType & IsoSurfaceTrianglesList);   // Output

        void CreateOutputPolyData(
                PointsListType & IsoSurfacePointsList,
                TrianglesListType & IsoSurfaceTrianglesList,
                vtkPolyData *OutputPolyData);   // Output

        void OrganizeIsoSurfacePointIds(
                PointsListType & IsoSurfacePointsList);   // Input and Output

        void OrganizeIsoSurfaceTrianglePointIds(
                PointsListType & IsoSurfacePointsList,
                TrianglesListType & IsoSurfaceTrianglesList);    // Input and Output

        void AddNormalsArray(vtkPolyData *OutputPolyData);

        void AddDeformationVectorDeviation(vtkPolyData *OutputPolyData);

        // Member Data
        StructureType::StructureModeType StructureMode; 
        bool SmoothingStatus;
        unsigned int SmoothingKernelSize;
        bool DebugStatus;
        bool ProgressStatus;

    private:
        IdentifyStructures(const IdentifyStructures &);  // Not implemented.
        void operator=(const IdentifyStructures &);   // Not implemented.
};

#endif
