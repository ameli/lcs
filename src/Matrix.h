/*
 * =====================================================================================
 *
 *       Filename:  Matrix.h
 *
 *    Description:  Matrix operations
 *
 *        Version:  1.0
 *        Created:  02/27/2014 01:00:43 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli 
 *   Organization:  University Of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __Matrix_h
#define __Matrix_h

// =====================
// Foreward Declarations
// =====================

#include <iostream>

// ============
// Matrix Class
// ============

class Matrix
{
    public:
        Matrix();
        ~Matrix();

        // static methods
        static void InitializeMatrix(
                double **&InputMatrix,
                unsigned int NumRows,
                unsigned int NumColumns);

        static void CopyMatrix(
                double **OriginalMatrix,
                unsigned int NumRows,
                unsigned int NumColumns,
                double **CopiedMatrix);   // Output

        static void DeleteMatrix(
                double **InputMatrix,
                unsigned int NumRows,
                unsigned int NumColumns);

        static std::ostream & DisplayMatrix(
                double **InputMatrix,
                unsigned int NumRows,
                unsigned int NumColumns,
                std::ostream & os);   // Output

        static std::ostream & DisplayVector(
                double *InputVector,
                unsigned int InputVectorLength,
                std::ostream & os);  // Output

        static void ComputeCovarianceMatrix(
                double **InputMatrix,
                unsigned int InputNumRows,
                unsigned int InputNumColumns,
                double **&CovarianceMatrix);  // Output

        static void MatrixVectorMultiplication(
                double **InputMatrix,
                unsigned int NumRows,
                unsigned int NumColumns,
                double *InputVector,
                double *OutputVector);  // Output

        static void MatrixVectorMultiplication(
                double **InputMatrix,
                unsigned int InputMatrixSize,
                double *InputVector,
                double *OutputVector);   // Output

        static void VectorMatrixMultiplication(
                double *InputVector,
                double **InputMatrix,
                unsigned int NumRows,
                unsigned int NumColumns,
                double *OutputVector);   // Output

        static void VectorMatrixMultiplication(
                double *InputVector,
                double **InputMatrix,
                unsigned int InputMatrixSize,
                double *OutputVector);   // Output

        static double VectorDotProduct(
                double *Vector1,
                double *Vector2,
                unsigned int VectorLength);

        static double GetVectorNorm(
                double *InputVector,
                unsigned int InputVectorLength);

        static double NormalizeVector(
                double *InputVector,    // Input and output
                unsigned int InputVectorLength);

        static void SolveLinearSystem(
                double **InputMatrix,
                unsigned int InputMatrixSize,
                double *KnownVector,
                double *UnknownVector);   // Output

        static void FindPrincipalDirection(
                double **InputMatrix,
                unsigned int NumRows,
                unsigned int NumColumns,
                double *PrincipalVector);  // Output

        static bool PowerMethodIteration(
                double **InputMatrix,
                unsigned int InputMatrixSize,
                double *GuessEigenVector,
                double *EigenVector);   // Output

        static bool RayleighQuotientIteration(
                double **InputMatrix,
                unsigned int InputMatrixSize,
                double *GuessEigenVector,
                double *EigenVector);   // Output

    private:
        static double FindMaxRelativeError(
                double *OldVector,
                double *NewVector,
                unsigned int vectorLength);
};

#endif
