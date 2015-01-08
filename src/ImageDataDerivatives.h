/*
 * =====================================================================================
 *
 *       Filename:  ImageDataDerivatives.h
 *
 *    Description:  Computes various derivatives of Image Data
 *
 *        Version:  1.0
 *        Created:  05/19/2014 03:18:58 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University Of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __ImageDataDerivatives_h
#define __ImageDataDerivatives_h

// ====================
// Forward Declarations
// ====================

class vtkDoubleArray;

// ======================
// Image Data Derivatives
// ======================

namespace ImageDataDerivatives
{
        inline void GetPointIndex(
                int *GridResolution,
                int PointId,
                unsigned int *PointIndex);  // Output

        inline int GetPointId(
                int *GridResolution,
                int *PointIndex);

        inline void GetStencilIds(
                int *GridResolution,
                int PointId,
                unsigned int TargetDirection,
                int *StencilIds);  // Output

    void ComputeHessianOfScalarFieldAtPoint(
            int *GridResolution,
            int PointId,
            vtkDoubleArray *InputScalarField,
            double **HesisanAtPoint);   // Output
}

#endif
