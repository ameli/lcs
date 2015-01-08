/*
 * =====================================================================================
 *
 *       Filename:  ImageDataDerivatives.cxx
 *
 *    Description:  Compute various derivatives of Image Data
 *
 *        Version:  1.0
 *        Created:  05/19/2014 03:18:33 PM
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

#include "ImageDataDerivatives.h"
#include <vtkDoubleArray.h>
#include <vtkStructuredData.h>
#include <cstdlib>   // for div

// ===========
// Definitions
// ===========

#ifndef DIMENSION
#define DIMENSION 3
#endif

// ===============
// Get Point Index
// ===============

// Decsription:
// Converts the list style ID if a point into coordinate-wise indices if the point
// in a structured grid. Grid Resolution is number of points in each direction.
// Note: PointId is 0-offset. So the first point is 0. Therefore, the PointIndex are
// also offset from 0. So the first point is (0,0,0) in 3D.

void ImageDataDerivatives::GetPointIndex(
        int *GridResolution,
        int PointId,
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

inline int ImageDataDerivatives::GetPointId(
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

void ImageDataDerivatives::GetStencilIds(
        int *GridResolution,
        int PointId,
        unsigned int TargetDirection,
        int *StencilIds)   // Output
{
    // Get point index
    unsigned int PointIndex[DIMENSION];
    ImageDataDerivatives::GetPointIndex(GridResolution,PointId,PointIndex);

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

// ========================================
// Compute Hessian Of Scalar Field At Point
// ========================================

// Description:
// This method computres Hessian matrix only for one given point. No smoothing filter 
// is used, whereas it's a simple finite differencing scheme is used.
// The output should be declared and initialized OUTSIDE, before calling this function.

void ImageDataDerivatives::ComputeHessianOfScalarFieldAtPoint(
        int *GridResolution,
        int PointId,
        vtkDoubleArray *InputScalarField,
        double **HessianAtPoint)   // Output
{
    // Iterate over Rows of Tensor
    for(unsigned int RowIterator = 0; RowIterator < DIMENSION; RowIterator++)
    {
        // Iterate over columns of Tensor
        for(unsigned int ColumnIterator = 0; ColumnIterator < DIMENSION; ColumnIterator++)
        {
            // Second order differentiation
            if(RowIterator == ColumnIterator)
            {
                // Get Stencil Ids
                int StencilIds[2];
                ImageDataDerivatives::GetStencilIds(GridResolution,PointId,RowIterator,StencilIds);

                #ifdef _OPENMP  // Thread safe region

                // Get Stencil point values
                double ForwardPointValue[1];
                double BackwardPointValue[1];
                double CentralPointValue[1];
                InputScalarField->GetTuple(StencilIds[1],ForwardPointValue);
                InputScalarField->GetTuple(StencilIds[0],BackwardPointValue);
                InputScalarField->GetTuple(PointId,CentralPointValue);

                // Finite Differencing
                HessianAtPoint[RowIterator][ColumnIterator] = 
                    BackwardPointValue[0] +        // Back   point of stencil
                    ForwardPointValue[0] +         // Front  point of stencil
                    (-2) * CentralPointValue[0];   // Center point of stencil

                #else  // Not thread safe region

                // Get Stencil point values
                double ForwardPointValue = InputScalarField->GetTuple1(StencilIds[1]);
                double BackwardPointValue = InputScalarField->GetTuple1(StencilIds[0]);
                double CentralPointValue = InputScalarField->GetTuple1(PointId);

                // Finite Differencing
                HessianAtPoint[RowIterator][ColumnIterator] = 
                    BackwardPointValue +        // Back   point of stencil
                    ForwardPointValue +         // Front  point of stencil
                    (-2) * CentralPointValue;   // Center point of stencil

                #endif

            }
            // First order differentiation in two directions
            else 
            {
                // Stencils in the middle row of the quad stencils
                int MiddleStencilIds[2];
                ImageDataDerivatives::GetStencilIds(GridResolution,PointId,RowIterator,MiddleStencilIds);

                // Four stencils of the right and left middle stencil in the transverse direction
                int RightStencilIds[2];
                int LeftStencilIds[2];
                ImageDataDerivatives::GetStencilIds(GridResolution,MiddleStencilIds[0],ColumnIterator,LeftStencilIds);
                ImageDataDerivatives::GetStencilIds(GridResolution,MiddleStencilIds[1],ColumnIterator,RightStencilIds);

                #ifdef _OPENMP // Thread Safe region

                // Get Stencil Point Values
                double UpperRightPointValue[1];
                double LowerRightPointValue[1];
                double UpperLeftPointValue[1];
                double LowerLeftPointValue[1];
                InputScalarField->GetTuple(RightStencilIds[1],UpperRightPointValue);  // Upper right point
                InputScalarField->GetTuple(RightStencilIds[0],LowerRightPointValue);  // Lower right point
                InputScalarField->GetTuple(LeftStencilIds[1],UpperLeftPointValue);    // Upper left  point
                InputScalarField->GetTuple(LeftStencilIds[0],LowerLeftPointValue);    // Lower left  point

                // Finite Differencing
                HessianAtPoint[RowIterator][ColumnIterator] =
                    (1/4)  * UpperRightPointValue[0] +   // Upper Right point
                    (-1/4) * LowerRightPointValue[0] +   // Lower Right point
                    (-1/4) * UpperLeftPointValue[0] +    // Upper Left  point
                    (1/4)  * LowerLeftPointValue[0];     // Lower Left  point

                #else  // Not thread safe region

                // Get Stencil Point Values
                double UpperRightPointValue = InputScalarField->GetTuple1(RightStencilIds[1]);  // Upper right point
                double LowerRightPointValue = InputScalarField->GetTuple1(RightStencilIds[0]);  // Lower right point
                double UpperLeftPointValue = InputScalarField->GetTuple1(LeftStencilIds[1]);    // Upper left point
                double LowerLeftPointValue = InputScalarField->GetTuple1(LeftStencilIds[0]);    // Lower left point

                // Finite Differencing
                HessianAtPoint[RowIterator][ColumnIterator] =
                    (1/4)  * UpperRightPointValue +   // Upper Right point
                    (-1/4) * LowerRightPointValue +   // Lower Right point
                    (-1/4) * UpperLeftPointValue +    // Upper Left  point
                    (1/4)  * LowerLeftPointValue;     // Lower Left  point

                #endif
            }
        }
    }
}
