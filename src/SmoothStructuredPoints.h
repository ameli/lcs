/*
 * =====================================================================================
 *
 *       Filename:  SmoothStructuredPoints.h
 *
 *    Description:  Smooth Structured Points
 *
 *        Version:  1.0
 *        Created:  04/17/2014 10:32:49 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University Of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __SmoothStructuredPoints_h
#define __SmoothStructuredPoints_h

// ======
// Macros
// ======

// Debugging
#ifndef HERE
#define HERE std::cout << "\033[0;44m" << "DEBUG: " << __FILE__ << " at " << __LINE__ << "\033[0m" << std::endl;
#endif

// ====================
// Forward Declarations
// ====================

// Complete declarations
#include <vtkImageAlgorithm.h>

// Incomplete declarations
class vtkDoubleArray;

// ==============================
// Smooth Structured Points Class
// ==============================

class SmoothStructuredPoints : public vtkImageAlgorithm
{
    public:
        static SmoothStructuredPoints *New();
        vtkTypeRevisionMacro(SmoothStructuredPoints,vtkImageAlgorithm);
        virtual void PrintSelf(ostream &os,vtkIndent Indent);

        // Static Methods
        static void SmoothScalarField(
                unsigned int KernelSize,
                int *GridResolution,
                vtkDoubleArray *InputScalarField,
                vtkDoubleArray *OutputScalarField);   // Output

        static void SelfSmoothScalarField(
                unsigned int KernelSize,
                int *GridResolution,
                vtkDoubleArray *ScalarField);   // Input / Output

        static void SmoothScalarFieldUseVTK(
                unsigned int KernelSize,
                int *GridResolution,
                vtkDoubleArray *InputScalarField,
                vtkDoubleArray *OutputScalarField);   // Output

        static void SmoothUnitVectorField(
                unsigned int KernelSize,
                int *GridResolution,
                vtkDoubleArray *InputVectorField,
                vtkDoubleArray *OutputVectorField);   // Output

        static void SelfSmoothUnitVectorField(
                unsigned int KernelSize,
                int *GridResolution,
                vtkDoubleArray *VectorField);   // input / Output

        void SmoothVectorFieldUseVTK(
                vtkDataSet *InputDataSet,
                vtkDoubleArray *VectorField);  // Input and Output

        static void SmoothSymmetricTensorField(
                unsigned int KernelSize,
                int *GridResolution,
                vtkDoubleArray *InputTensorField,
                vtkDoubleArray *OutputTensorField);   // Output

        static void SelfSmoothSymmetricTensorField(
                unsigned int KernelSize,
                int *GridResolution,
                vtkDoubleArray *TensorField);   // Input / Output

        static void SmoothGradientOfScalarField(
                unsigned int KernelSize,
                int *GridResolution,
                vtkDoubleArray *InputScalarField,
                vtkDoubleArray *OutputGradientField);   // Output

        static void SmoothHessianAtPoint(
                unsigned int KernelSize,
                int *GridResolution,
                int PointId,
                vtkDoubleArray *InputScalarField,
                double **TensorAtPoint);   // Output

    protected:
        SmoothStructuredPoints();
        virtual ~SmoothStructuredPoints();

        // Pipeline Executives
        virtual int RequestData(
                vtkInformation *Request,
                vtkInformationVector **InputVector,
                vtkInformationVector *OutputVector);

    private:
        SmoothStructuredPoints(const SmoothStructuredPoints &);   // Not implemented
        void operator=(const SmoothStructuredPoints &);   // Not implemented
};

// ========================
// Structured Points Helper
// ========================

#ifndef __StructuredPointsHelper
#define __StructuredPointsHelper

namespace StructuredPointsHelper
{
    inline void GetPointIndex(
            int *GridResolution,
            int PointId,
            int *PointIndex);  // Output

    inline int GetPointId(
            int *GridResolution,
            int *PointIndex);

    inline void GetStencilIds(
            int *GridResolution,
            unsigned int PointId,
            unsigned int TargetDirection,
            unsigned int *StencilIds);  // Output

    inline void ClipIndexToBounds(
            int RightBound,
            int &Index);   // Input / Output

    inline void TileIndexToBounds(
            int RightBound,
            int &Index);   // Input / Output
}

#endif

// ================
// Smoothing Helper
// ================

namespace SmoothingHelper
{
    inline double GaussFunction(
            double Distance,
            double StandardDeviation,
            unsigned int Dimension);

    inline double GaussFunctionFirstDerivative(
            double Distance,
            double StandardDeviation,
            unsigned int Dimension);

    inline double GaussFunctionSecondDerivative(
            double Distance,
            double StandardDeviation,
            unsigned int Dimension);

    void InitializeSmoothingKernel(
            unsigned int KernelSize,
            double *SmoothingKernel);   // Output

    void InitializeSmoothingGradientKernel(
            unsigned int KernelSize,
            double *SmoothingGradientKernel);   // Output

    void InitializeSmoothingHessianKernel(
            unsigned int KernelSize,
            double *SmoothingHessianKernel);   // Output

    void ConvoluteKernelInOneDirection(
            double *Kernel,
            unsigned int KernelSize,
            unsigned int Direction,
            int *GridResolution, 
            vtkDoubleArray *OriginalScalarField,
            vtkDoubleArray *ConvolutedScalarField);   // Output

    void ConvertCartesianToSphericalCoordinates(
            double *CartesianVector,
            double SphericalVector[2]);   // Output

    void ConvertSphericalToCartesianCoordinates(
            double *SphericalVector,
            double CartesianVector[3]);   // Output

    double Convolute3DKernelAtPoint(
            unsigned int KernelSize,
            double **Kernels,
            int *GridResolution,
            int PointId,
            vtkDoubleArray *ScalarField);
}

#endif
