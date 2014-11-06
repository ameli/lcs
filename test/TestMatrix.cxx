/*
 * =====================================================================================
 *
 *       Filename:  TestMatrix.cxx
 *
 *    Description:  Test for Matrix.cxx
 *
 *        Version:  1.0
 *        Created:  02/27/2014 07:NumRows2:04 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University Of Caifornia, Berkeley
 *
 * =====================================================================================
 */

// =======
// Headers
// =======

#include "Matrix.h"
#include <iostream>
#include <math.h>

// ====
// Main
// ====

int main(int argc,char *argv[])
{
    double **Array;
    unsigned int NumRows = 6;
    unsigned int NumColumns = 4;
    Matrix::InitializeMatrix(Array,NumRows,NumColumns);

    // Hilbert Matrix
    for(unsigned int i=0;i<NumRows;i++)
    {
        for(unsigned int j=0;j<NumColumns;j++)
        {
            Array[i][j] = double(1.0/(2*i+j+1));
        }
    }

    // Display Matrix
    std::cout << "Array: " << std::endl;
    Matrix::DisplayMatrix(Array,NumRows,NumColumns,std::cout);

    // Covariance Matrix
    double **CovarianceArray;
    Matrix::ComputeCovarianceMatrix(Array,NumRows,NumColumns,CovarianceArray);

    std::cout << "\nCovariance Matrix: " << std::endl;
    Matrix::DisplayMatrix(CovarianceArray,NumColumns,NumColumns,std::cout);

    // Define a vector
    double *InputVector = new double[NumColumns];
    for(unsigned int i = 0; i < NumColumns; i++)
    {
        InputVector[i] = double(i);
    }
    std::cout << "\nVector1: " << std::endl;
    Matrix::DisplayVector(InputVector,NumColumns,std::cout);

    // Matrix Vector Multiplication
    double *OutputVector = new double[NumRows];
    Matrix::MatrixVectorMultiplication(Array,NumRows,NumColumns,InputVector,OutputVector);
    std::cout << "\nMultiplied Array * Vector1: " << std::endl;
    Matrix::DisplayVector(OutputVector,NumRows,std::cout);

    // Vector Matrix Multiplication
    double *InputVector2 = new double[NumRows];
    for(int i = 0; i<static_cast<int>(NumRows); i++)
    {
        InputVector2[i] = double(i);
    }
    double *OutputVector2 = new double[NumColumns];
    Matrix::VectorMatrixMultiplication(InputVector2,Array,NumRows,NumColumns,OutputVector2);
    std::cout << "\nVector2: " << std::endl;
    Matrix::DisplayVector(InputVector2,NumRows,std::cout);
    std::cout << "Multiplied Vector2 * Array: " << std::endl;
    Matrix::DisplayVector(OutputVector2,NumColumns,std::cout);

    // Vector Norm
    double Norm = Matrix::GetVectorNorm(InputVector,NumColumns);
    std::cout << "\nVector1 Norm: " << Norm << std::endl;

    // Normalize Vector
    Matrix::NormalizeVector(InputVector,NumColumns);
    std::cout << "\nNormalized Vector1: " << std::endl;
    Matrix::DisplayVector(InputVector,NumColumns,std::cout);

    // Principal Direction of matrix
    double *PrincipalDirection = new double[NumColumns];
    Matrix::FindPrincipalDirection(Array,NumRows,NumColumns,PrincipalDirection);
    std::cout << "\nPrincipal Direction of Array rows: " << std::endl;
    Matrix::DisplayVector(PrincipalDirection,NumColumns,std::cout);

    // Free memory
    Matrix::DeleteMatrix(Array,NumRows,NumColumns);
    delete [] InputVector;
    delete [] InputVector2;
    delete [] OutputVector;
    delete [] OutputVector2;
    delete [] PrincipalDirection;

    // Test Solve Linear System
    std::cout << "\nSolution of Ax = b: " << std::endl;

    int ArraySize = 5;
    double **Array2 = new double*[ArraySize];
    double *KnownVector = new double[ArraySize];
    double *UnknownVector = new double[ArraySize];

    // Initialize Array and Known Vector
    for(int i=0; i< ArraySize; i++)
    {
        Array2[i] = new double[ArraySize];
        for(int j = 0; j < ArraySize; j++)
        {
            // Create a well conditioned number array
            // Array2[i][j] = sin((i+1)*(j+1));

            // Create non-symetric ill-conditioned matrix
            Array2[i][j] = 1.0 / (2*i+j+1);
        }

        KnownVector[i] = double(i);
    }

    // Display Matrix and Known Vector
    std::cout << "\nArray2: " << std::endl;
    Matrix::DisplayMatrix(Array2,ArraySize,ArraySize,std::cout);

    std::cout << "\nKnown Vector: " << std::endl;
    Matrix::DisplayVector(KnownVector,ArraySize,std::cout);

    // Solve Linear System
    Matrix::SolveLinearSystem(Array2,ArraySize,KnownVector,UnknownVector);

    std::cout << "\nUnknown Vector: " << std::endl;
    Matrix::DisplayVector(UnknownVector,ArraySize,std::cout);

    // Multiply Solution for test
    double *Error = new double[ArraySize];
    Matrix::MatrixVectorMultiplication(Array2,ArraySize,ArraySize,UnknownVector,Error);

    // Subtract form actual Known Vector
    for(int i = 0; i<ArraySize; i++)
    {
        Error[i] -= KnownVector[i];
    }

    std::cout << "\nError of Solution: " << std::endl;
    Matrix::DisplayVector(Error,ArraySize,std::cout);

    // Free Memory
    Matrix::DeleteMatrix(Array2,ArraySize,ArraySize);
    delete [] KnownVector;
    delete [] UnknownVector;
    delete [] Error;

    // Test for finding PCA for matrix that has very close large eigenvalues
    std::cout << "\nTest array of close largest eigenvalues: " << std::endl;
    double **Array3;
    Matrix::InitializeMatrix(Array3,3,3);

    // The eigenvalues of Array3 has 3.00 and 3.05 with orthogonal eigenvectors
    Array3[0][0] =  1.30683; Array3[0][1] =  -0.61632; Array3[0][2] =  -0.39064;
    Array3[1][0] = -0.61632; Array3[1][1] =   2.83119; Array3[1][2] =  -0.13139;
    Array3[2][0] = -0.39064; Array3[2][1] =  -0.13139; Array3[2][2] =   2.91197;
    std::cout << "\nArray3: " << std::endl;
    Matrix::DisplayMatrix(Array3,3,3,std::cout);

    // Using Power method
    double GuessEigenVector31[3] = {1,0,0};
    double PrincipalVector31[3];
    Matrix::PowerMethodIteration(Array3,3,GuessEigenVector31,PrincipalVector31);
    std::cout << "\nPrincipal eigenvector with Power Method: " << std::endl;
    Matrix::DisplayVector(PrincipalVector31,3,std::cout);

    // Using Rayleigh Quotient Iteration
    double *GuessEigenVector32 = PrincipalVector31; // using previous result for guess
    double PrincipalVector32[3];
    Matrix::RayleighQuotientIteration(Array3,3,GuessEigenVector32,PrincipalVector32);
    std::cout << "\nPrincipal eigenvector with Rayleigh Quotient Iteration: " << std::endl;
    Matrix::DisplayVector(PrincipalVector32,3,std::cout);

    // The overall method
    std::cout << "\nPrincipal eigenvector with overall method: " << std::endl;
    double PrincipalEigenVector[3];
    Matrix::FindPrincipalDirection(Array3,3,3,PrincipalEigenVector);
    Matrix::DisplayVector(PrincipalEigenVector,3,std::cout);

    // Free memory
    Matrix::DeleteMatrix(Array3,3,3);

    // Test Another matrix
    std::cout << "\nTest another matrix: " << std::endl;
    double **Array4;
    Matrix::InitializeMatrix(Array4,8,3);
    double TempArray[8][3] = {
        {0.0176366, -0.990052, 0.13959},
        {0.534174, -0.828824, 0.16646},
        {0.584825, 0.316071, -0.747047},
        {0.59031, -0.802807, 0.0838722},
        {0.044219, 0.0216822, 0.998787},
        {0.0054937, 0.827568, -0.561338},
        {0.230108, 0.676782, -0.699297},
        {0.0985213, -0.852199, 0.513858}};
    for(int i =0;i<8;i++)
    {
        for(int j=0; j<3;j++)
        {
            Array4[i][j] = TempArray[i][j];
        }
    }
    Matrix::DisplayMatrix(Array4,8,3,std::cout);

    // Covariance Matrix
    double **CovarianceMatrix;
    Matrix::ComputeCovarianceMatrix(Array4,8,3,CovarianceMatrix);
    std::cout << "\nCovarianceMatrix: " << std::endl;
    Matrix::DisplayMatrix(CovarianceMatrix,3,3,std::cout);

    // Principal vector
    double PrincipalVector4[3];
    double GuessEigenVector[3] = {1,0,0};
    std::cout << "\nPrincipal direction using Power method: " << std::endl;
    Matrix::PowerMethodIteration(CovarianceMatrix,3,GuessEigenVector,PrincipalVector4);
    Matrix::DisplayVector(PrincipalVector4,3,std::cout);

    std::cout << "\nPrincipal direction using Rayleigh Iteration: " << std::endl;
    Matrix::RayleighQuotientIteration(CovarianceMatrix,3,PrincipalVector4,PrincipalVector4);
    Matrix::DisplayVector(PrincipalVector4,3,std::cout);

    Matrix::FindPrincipalDirection(Array4,8,3,PrincipalVector4);
    std::cout << "\nPrincipal Direction using overall method: " << std::endl;
    Matrix::DisplayVector(PrincipalVector4,3,std::cout);

    // Free Memory
    Matrix::DeleteMatrix(Array4,8,3);
    Matrix::DeleteMatrix(CovarianceMatrix,3,3);

    return 0;
}
