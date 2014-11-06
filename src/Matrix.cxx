/*
 * =====================================================================================
 *
 *       Filename:  Matrix.cxx
 *
 *    Description:  Matrix operations
 *
 *        Version:  1.0
 *        Created:  02/27/2014 01:00:22 AM
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

#include "Matrix.h"

// STL
#include <cstddef>  // for NULL
#include <cstdlib>  // for exit
#include <iomanip>  // for ostream format
#include <cmath>    // for sqrt, pow

// ======
// Macros
// ======

#define EPSILON 1e-5

// ===========
// Constructor
// ===========

Matrix::Matrix()
{
}

// ==========
// Destructor
// ==========

Matrix::~Matrix()
{
}

// =================
// Initialize Matrix
// =================

// Description:
// Allocate memory for 2D arrays.

void Matrix::InitializeMatrix(
        double **&InputMatrix,
        unsigned int NumRows,
        unsigned int NumColumns)
{
    InputMatrix = new double*[NumRows];

    for(unsigned int RowIterator = 0; RowIterator < NumRows; RowIterator++)
    {
        InputMatrix[RowIterator] = new double[NumColumns];
    }
}

// ===========
// Copy Matrix
// ===========

void Matrix::CopyMatrix(
        double **OriginalMatrix,
        unsigned int NumRows,
        unsigned int NumColumns,
        double **CopiedMatrix)
{
    for(unsigned int RowIterator = 0; RowIterator < NumRows; RowIterator++)
    {
        for(unsigned int ColumnIterator = 0; ColumnIterator < NumColumns; ColumnIterator++)
        {
            CopiedMatrix[RowIterator][ColumnIterator] = OriginalMatrix[RowIterator][ColumnIterator];
        }
    }
}

// =============
// Delete Matrix
// =============

void Matrix::DeleteMatrix(
        double **InputMatrix,
        unsigned int NumRows,
        unsigned int NumColumns)
{
    for(unsigned int RowIterator = 0; RowIterator < NumRows; RowIterator++)
    {
        delete [] InputMatrix[RowIterator];
    }

    delete [] InputMatrix;
    InputMatrix = NULL;
}

// ==============
// Display Matrix
// ==============

std::ostream & Matrix::DisplayMatrix(
        double **InputMatrix,
        unsigned int NumRows,
        unsigned int NumColumns,
        std::ostream & os)
{
    for(unsigned int RowIterator = 0; RowIterator < NumRows; RowIterator++)
    {
        for(unsigned int ColumnIterator = 0; ColumnIterator < NumColumns; ColumnIterator++)
        {
            os << std::setw(8) << std::right << InputMatrix[RowIterator][ColumnIterator];

            // Adding comma
            if(ColumnIterator < NumColumns-1)
            {
                os << std::setw(3) << std::left << ",";
            }
        }

        // New line
        os << std::endl;
    }

    return os;
}

// ==============
// Display Vector
// ==============

std::ostream & Matrix::DisplayVector(
        double *InputVector,
        unsigned int InputVectorLength,
        std::ostream & os)
{
    for(unsigned int Iterator = 0; Iterator < InputVectorLength; Iterator++)
    {
        os << std::setw(8) << std::right << InputVector[Iterator];

        // Add comma
        if(Iterator < InputVectorLength-1)
        {
            os << std::setw(3) << std::left << ",";
        }
    }

    // New line
    os << std::endl;

    return os;
}

// =========================
// Compute Covariance Matrix
// =========================

// Description:
// Covaiance Matrix has size of NumColumns * NumColumns;
// Output is allocated inside this function.

void Matrix::ComputeCovarianceMatrix(
        double **InputMatrix,
        unsigned int InputNumRows,
        unsigned int InputNumColumns,
        double **&CovarianceMatrix)
{
    // Covariane Matrix rows and columns
    unsigned int CovarianceNumRows = InputNumColumns;
    unsigned int CovarianceNumColumns = InputNumColumns;

    // Allocate Output
    Matrix::InitializeMatrix(CovarianceMatrix,CovarianceNumRows,CovarianceNumColumns);

    // Compute Each element
    for(unsigned int RowIterator = 0; RowIterator < CovarianceNumRows; RowIterator++)
    {
        for(unsigned int ColumnIterator = 0; ColumnIterator < CovarianceNumColumns; ColumnIterator++)
        {
            double Sum = 0;

            // Compute upper right half of symmetric matrix
            if(ColumnIterator >= RowIterator)
            {
                for(unsigned int DummyIterator = 0; DummyIterator < InputNumRows; DummyIterator++)
                {
                    Sum +=  InputMatrix[DummyIterator][RowIterator] * InputMatrix[DummyIterator][ColumnIterator];
                }
                CovarianceMatrix[RowIterator][ColumnIterator] = Sum;
            }

            // Use upper right half of matrix for lower left half
            else
            {
                CovarianceMatrix[RowIterator][ColumnIterator] = CovarianceMatrix[ColumnIterator][RowIterator];
            }
        }
    }
}

// ===========================
// Matrix Vector Multipication
// ===========================

// Rectanglular matrices
void Matrix::MatrixVectorMultiplication(
        double **InputMatrix,
        unsigned int NumRows,
        unsigned int NumColumns,
        double *InputVector,
        double *OutputVector)
{
    // Multiplication
    for(unsigned int RowIterator = 0; RowIterator < NumRows; RowIterator++)
    {
        OutputVector[RowIterator] = 0;
        for(unsigned int ColumnIterator = 0; ColumnIterator < NumColumns; ColumnIterator++)
        {
            OutputVector[RowIterator] += InputMatrix[RowIterator][ColumnIterator] * InputVector[ColumnIterator];
        }
    }
}

// Square matrices
void Matrix::MatrixVectorMultiplication(
        double **InputMatrix,
        unsigned int InputMatrixSize,
        double *InputVector,
        double *OutputVector)   // Output
{
    Matrix::MatrixVectorMultiplication(
            InputMatrix,
            InputMatrixSize,
            InputMatrixSize,
            InputVector,
            OutputVector);
}

// ============================
// Vector Matrix Multiplication
// ============================

// Rectangular Matrices
void Matrix::VectorMatrixMultiplication(
        double *InputVector,
        double **InputMatrix,
        unsigned int NumRows,
        unsigned int NumColumns,
        double *OutputVector)
{
    for(unsigned int ColumnIterator = 0; ColumnIterator < NumColumns; ColumnIterator++)
    {
        OutputVector[ColumnIterator] = 0;
        for(unsigned int RowIterator = 0; RowIterator < NumRows; RowIterator++)
        {
            OutputVector[ColumnIterator] += InputVector[RowIterator] * InputMatrix[RowIterator][ColumnIterator];
        }
    }
}

// Square Matrices
void Matrix::VectorMatrixMultiplication(
        double *InputVector,
        double **InputMatrix,
        unsigned int InputMatrixSize,
        double *OutputVector)   // Output
{
    Matrix::VectorMatrixMultiplication(
            InputVector,
            InputMatrix,
            InputMatrixSize,
            InputMatrixSize,
            OutputVector);
}

// ==================
// Vector Dot Product
// ==================

double Matrix::VectorDotProduct(
        double *Vector1,
        double *Vector2,
        unsigned int VectorLength)
{
    double Sum = 0;
    for(unsigned int Iterator = 0; Iterator < VectorLength; Iterator++)
    {
        Sum += Vector1[Iterator] * Vector2[Iterator];
    }

    return Sum;
}

// ===============
// Get Vector Norm
// ===============

double Matrix::GetVectorNorm(
        double *InputVector,
        unsigned int InputVectorLength)
{
    // Using dot product
    double SquareNorm = Matrix::VectorDotProduct(InputVector,InputVector,InputVectorLength);

    return sqrt(SquareNorm);
}

// ================
// Normalize Vector
// ================

// Description:
// This function over-writes the input vector.

double Matrix::NormalizeVector(
        double *InputVector,
        unsigned int InputVectorLength)
{
    // Norm
    double Norm = Matrix::GetVectorNorm(InputVector,InputVectorLength);

    // Normalize
    for(unsigned int Iterator = 0; Iterator < InputVectorLength; Iterator++)
    {
        InputVector[Iterator] /= Norm;
    }

    return Norm;
}

// ===================
// Solve Linear System
// ===================

// Description:
// Solves the linear system of equations using Gaussian Elimination with complete Pivoting.
// This method is only efficient for small matrices. The computational cost is O(n^3).
// This method changes InputMatrix and KnownVector.
// Use this method for small matrics, e.g. 3x3.

void Matrix::SolveLinearSystem(
        double **InputMatrix,
        unsigned int InputMatrixSize,
        double *KnownVector,
        double *UnknownVector)   // Output
{
    // Ordering of indecides of output vector due to re-arrangement of pivoting
    unsigned int *IndexOrder = new unsigned int[InputMatrixSize];
    for(unsigned int Iterator = 0; Iterator < InputMatrixSize; Iterator++)
    {
        IndexOrder[Iterator] = Iterator;
    }

    // Computing upper triangular matrix (PLU decomosition)
    // loop over Rows
    for(unsigned int RowIterator = 0; RowIterator < InputMatrixSize-1; RowIterator++)
    {
        // Find max of first column for pivoting
        unsigned int PivotRow = RowIterator;
        double MaxPivot = InputMatrix[PivotRow][RowIterator];

        // Search over first column of each iteration
        for(unsigned int RowIterator2 = RowIterator+1; RowIterator2 < InputMatrixSize; RowIterator2++)
        {
            if(fabs(InputMatrix[RowIterator2][RowIterator]) > fabs(MaxPivot))
            {
                // Update max value and index
                PivotRow = RowIterator2;
                MaxPivot = InputMatrix[PivotRow][RowIterator];
            }
        }

        // Switch Rows for pivoting
        if(PivotRow != RowIterator)
        {
            // Switch pointers of InputMatrix Rows
            double *TempMatrixRow = InputMatrix[RowIterator];
            InputMatrix[RowIterator] = InputMatrix[PivotRow];
            InputMatrix[PivotRow] = TempMatrixRow;

            // Switch values Known Vector entry
            double TempVectorEntry = KnownVector[RowIterator];
            KnownVector[RowIterator] = KnownVector[PivotRow];
            KnownVector[PivotRow] = TempVectorEntry;

            // Keep record of swapping rows
            IndexOrder[RowIterator] = PivotRow;
            IndexOrder[PivotRow] = RowIterator;
        }

        // Iterator over rest of rows for Gaussian elimination
        for(unsigned int RowIterator3 = RowIterator+1; RowIterator3 < InputMatrixSize; RowIterator3++)
        {
            double RowRatio = double(InputMatrix[RowIterator3][RowIterator]) / InputMatrix[RowIterator][RowIterator];

            // Iterate over columns
            for(unsigned int ColumnIterator = RowIterator; ColumnIterator < InputMatrixSize; ColumnIterator++)
            {
                InputMatrix[RowIterator3][ColumnIterator] -= RowRatio * InputMatrix[RowIterator][ColumnIterator];
            }

            // Update Known vector
            KnownVector[RowIterator3] -= RowRatio * KnownVector[RowIterator];
        }
    }

    // Backward substitution for obtaing Unknown vector (solution)
    for(int RowIterator = InputMatrixSize-1; RowIterator >= 0 ; RowIterator--)
    {
        // Sum over all previous UnknownVectors
        double Sum = 0;
        if(RowIterator < static_cast<int>(InputMatrixSize - 1))
        {
            for(int ColumnIterator = RowIterator+1;
                ColumnIterator < static_cast<int>(InputMatrixSize);
                ColumnIterator++)
            {
                Sum += InputMatrix[RowIterator][ColumnIterator] * UnknownVector[ColumnIterator];
            }
        }

        // Solution
        UnknownVector[RowIterator] = 
            double((KnownVector[RowIterator] - Sum)) / InputMatrix[RowIterator][RowIterator];
    }

    // Free memory
    delete [] IndexOrder;
}

// ========================
// Find Principal Direction
// ========================

// Description:
// Principal Component Analysis to find the "first" principal direction for
// some vectors. The vectors are stored as ROWS of InputMatrix.
// InputNumRows: number of vectors.
// InputNumColumn: dimension of each vector

void Matrix::FindPrincipalDirection(
        double **InputMatrix,
        unsigned int InputNumRows,
        unsigned int InputNumColumns,
        double *PrincipalVector)   // Output
{
    // Covariance Matrix
    double **CovarianceMatrix;
    Matrix::ComputeCovarianceMatrix(InputMatrix,InputNumRows,InputNumColumns,CovarianceMatrix);
    unsigned int CovarianceMatrixSize = InputNumColumns;

    // Arbitrary Vector
    double *GuessEigenVector = new double[CovarianceMatrixSize];
    GuessEigenVector[0] = 1.0;
    for(unsigned int Iterator = 1; Iterator < CovarianceMatrixSize; Iterator++)
    {
        GuessEigenVector[Iterator] = 0.0;
    }

    // First try: Power  Method Iteration
    bool ConvergenceStatus = Matrix::PowerMethodIteration(
            CovarianceMatrix,
            CovarianceMatrixSize,
            GuessEigenVector,
            PrincipalVector);

    // Use Output of Power method for Guess EigenVector
    MatrixVectorMultiplication(CovarianceMatrix,CovarianceMatrixSize,PrincipalVector,GuessEigenVector);

    // Rayleight Quotient method (Cubic vonvergence)
    ConvergenceStatus = Matrix::RayleighQuotientIteration(
            CovarianceMatrix,
            CovarianceMatrixSize,
            GuessEigenVector,
            PrincipalVector);

    if(ConvergenceStatus == false)
    {
        std::cerr << "Warning: Eigenvalue computation can not converge." << std::endl;
    }

    // Free Memory
    Matrix::DeleteMatrix(CovarianceMatrix,CovarianceMatrixSize,CovarianceMatrixSize);
    delete [] GuessEigenVector;
}

// ======================
// Power Method Iteration
// ======================

bool Matrix::PowerMethodIteration(
        double **InputMatrix,
        unsigned int InputMatrixSize,
        double *GuessEigenVector,
        double *EigenVector)
{
    // Old EigenVector
    double *OldEigenVector = new double[InputMatrixSize];
    for(unsigned int Iterator = 0; Iterator < InputMatrixSize; Iterator++)
    {
        OldEigenVector[Iterator] = GuessEigenVector[Iterator];
    }
    Matrix::NormalizeVector(OldEigenVector,InputMatrixSize);

    // Variables for iteration
    double *NewEigenVector = new double[InputMatrixSize];
    double MaxRelativeError = 1;
    unsigned int TerminationCounter = 0;
    unsigned int TerminationLimit = 10;
    bool ConvergenceStatus = true;

    // Iteration loop
    while(MaxRelativeError > EPSILON)
    {
        // Terminate if not converge after 50 iterations
        TerminationCounter++;
        if(TerminationCounter > TerminationLimit)
        {
            ConvergenceStatus = false;
            break;
        }

        // Multiply Matrix to old vector
        Matrix::MatrixVectorMultiplication(
                InputMatrix,
                InputMatrixSize,
                OldEigenVector,
                NewEigenVector);

        // Normalize new vector
        Matrix::NormalizeVector(NewEigenVector,InputMatrixSize);
        
        // Find Max Relative Error
        MaxRelativeError = Matrix::FindMaxRelativeError(OldEigenVector,NewEigenVector,InputMatrixSize);

        // Update EigenVector
        double *TempEigenVector = OldEigenVector;
        OldEigenVector = NewEigenVector;
        NewEigenVector = TempEigenVector;
    }

    // Output
    for(unsigned int Iterator = 0; Iterator < InputMatrixSize; Iterator++)
    {
        // Use old Eigenvector for output, becasue New and old were swapped.
        EigenVector[Iterator] = OldEigenVector[Iterator];
    }

    // Free Memory
    delete [] OldEigenVector;
    delete [] NewEigenVector;

    return ConvergenceStatus;

}

// =======================
// Find Max Relative Error
// =======================

// Description:
// This is not actually the max relative error. The Max absolute error is computed
// and then its relative error is returned.

double Matrix::FindMaxRelativeError(
        double *OldVector,
        double *NewVector,
        unsigned int VectorLength)
{
    double MaxAbsoluteError = 0;
    unsigned int MaxAbsoluteErrorIndex = 0;

    for(unsigned int Iterator = 0; Iterator < VectorLength; Iterator++)
    {
        double TempAbsoluteError = fabs(NewVector[Iterator] - OldVector[Iterator]);

        // Compare error with previous component error
        if(TempAbsoluteError > MaxAbsoluteError)
        {
            MaxAbsoluteError = TempAbsoluteError;
            MaxAbsoluteErrorIndex = Iterator;
        }
    }

    // Relative Error
    // Note: Relative error is w.r.t Vector"Old", not VectorNew.
    double MaxRelativeError = MaxAbsoluteError / (OldVector[MaxAbsoluteErrorIndex] + EPSILON);

    return MaxRelativeError;
}


// ===========================
// Rayleigh Quotient Iteration
// ===========================

// Description:
// This method is cubically fast convergent for symmetric matrices.
// Since It uses Gaussian elimination for solving Linear system, this is 
// not recommended for large matrices.
// Note: Tghis method follows the GuessEigenVector, i.e. it converges to the
// Eigenvector that is closer to the guess. Use the output of Power Iterartion
// as a good guess for Rayleight Quotient Method.

bool Matrix::RayleighQuotientIteration(
        double **InputMatrix,
        unsigned int InputMatrixSize,
        double *GuessEigenVector,
        double *EigenVector)
{
    // Old Rayleigh Quotient (EigenValue)
    Matrix::NormalizeVector(GuessEigenVector,InputMatrixSize);
    double *TempVector = new double[InputMatrixSize];
    Matrix::MatrixVectorMultiplication(InputMatrix,InputMatrixSize,GuessEigenVector,TempVector);
    double OldRayleighQuotient = Matrix::VectorDotProduct(GuessEigenVector,TempVector,InputMatrixSize);
    delete [] TempVector;

    // Old EigenVector
    double *OldEigenVector = new double[InputMatrixSize];
    for(unsigned int Iterator = 0; Iterator < InputMatrixSize; Iterator++)
    {
        OldEigenVector[Iterator] = GuessEigenVector[Iterator];
    }
    Matrix::NormalizeVector(OldEigenVector,InputMatrixSize);

    // Variables for iteration
    double NewRayleighQuotient;
    double *NewEigenVector = new double[InputMatrixSize];
    double RelativeError = 1;
    unsigned int TerminationCounter = 0;
    unsigned int TerminationLimit = 10;
    bool ConvergenceStatus = true;
    double **CopiedMatrix;
    Matrix::InitializeMatrix(CopiedMatrix,InputMatrixSize,InputMatrixSize);

    // Iteration loop
    while(RelativeError > EPSILON)
    {
        // Terminate long iterations
        TerminationCounter++;
        if(TerminationCounter > TerminationLimit)
        {
            ConvergenceStatus = false;
            break;
        } 

        // Create Matrix for Linear system
        Matrix::CopyMatrix(InputMatrix,InputMatrixSize,InputMatrixSize,CopiedMatrix);

        // Subtract diagonals from Rayleigh quotient
        for(unsigned int Iterator = 0; Iterator < InputMatrixSize; Iterator++)
        {
            CopiedMatrix[Iterator][Iterator] -= OldRayleighQuotient;
        }

        // Solve linear system for new eigenvector
        Matrix::SolveLinearSystem(CopiedMatrix,InputMatrixSize,OldEigenVector,NewEigenVector);
        Matrix::NormalizeVector(NewEigenVector,InputMatrixSize);

        // New Rayleigh quotient
        double *TempVector = new double[InputMatrixSize];
        Matrix::MatrixVectorMultiplication(InputMatrix,InputMatrixSize,NewEigenVector,TempVector);
        NewRayleighQuotient = Matrix::VectorDotProduct(NewEigenVector,TempVector,InputMatrixSize);
        delete [] TempVector;

        // Relative Error
        double AbsoluteError = fabs(NewRayleighQuotient - OldRayleighQuotient);
        RelativeError = AbsoluteError / fabs(OldRayleighQuotient + EPSILON);

        // Update Rayleigh Quotient
        OldRayleighQuotient = NewRayleighQuotient;

        // Update EigenVector (Swap pointer of old and new)
        double *TempEigenVector = OldEigenVector;
        OldEigenVector = NewEigenVector;
        NewEigenVector = TempEigenVector;
    }

    // Output
    for(unsigned int Iterator = 0; Iterator < InputMatrixSize; Iterator++)
    {
        // Use old Eigenvector for output, becasue New and old were swapped.
        EigenVector[Iterator] = OldEigenVector[Iterator];
    }

    // Delete Memory
    Matrix::DeleteMatrix(CopiedMatrix,InputMatrixSize,InputMatrixSize);
    delete [] OldEigenVector;
    delete [] NewEigenVector;

    return ConvergenceStatus;
}
