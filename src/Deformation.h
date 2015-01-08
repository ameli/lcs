/*
 * =====================================================================================
 *
 *       Filename:  Deformation.h
 *
 *    Description:  Computes Deformation of flowmap
 *
 *        Version:  1.0
 *        Created:  02/12/2014 05:20:13 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University Of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __Deformation_h
#define __Deformation_h

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

// =====================
// Foreward Declarations
// =====================

// Complete declarations
#include <vtkDataSetAlgorithm.h>

// Incomplete declarations
class vtkDataArrayCollection;
class vtkDoubleArray;

// =================
// Deformation Class
// =================

class Deformation : public vtkDataSetAlgorithm
{
    public:
        static Deformation *New();
        vtkTypeRevisionMacro(Deformation,vtkDataSetAlgorithm);
        virtual void PrintSelf(ostream &os,vtkIndent Indent);

        // Accessors / Mutators
        vtkGetMacro(SmoothingStatus,bool);
        vtkSetMacro(SmoothingStatus,bool);
        void SmoothingOn() { this->SmoothingStatus = true; };
        void SmoothingOff() { this->SmoothingStatus = false; };

        vtkGetMacro(SmoothingKernelSize,unsigned int);
        vtkSetMacro(SmoothingKernelSize,unsigned int);

        vtkGetMacro(ProgressStatus,bool);
        vtkSetMacro(ProgressStatus,bool);
        void ProgressOn() { this->ProgressStatus = true; };
        void ProgressOff() { this->ProgressStatus = false; };

    protected:
        Deformation();
        virtual ~Deformation();

        // Pipeline Executives
        virtual int FillOutputPortInformation(
                int port,
                vtkInformation *info);

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
        void ComputeEigenValuesAndVectors(
                vtkDataSet *InputDataSet,
                vtkDoubleArray *EigenValues,                      // Output
                vtkDataArrayCollection *EigenVectorsCollection);  // Output

        void SortEigenValuesAndVectors(
                double *TensorEigenValues,    // Input and Output
                double **TensorEigenVectors,  // Input and Output
                unsigned int TensorSize);

        void ComputeDeformationValuesAndVectors(
                vtkDoubleArray *EigenValues,
                vtkDataArrayCollection *EigenVectorsCollection,
                vtkDoubleArray *MaxStrainValues,     // Output
                vtkDoubleArray *MaxStrainVectors,    // Output
                vtkDoubleArray *MinStrainValues,     // Output
                vtkDoubleArray *MinStrainVectors,    // Output
                vtkDoubleArray *ShearValues,         // Output
                vtkDoubleArray *ShearVectors1,       // Output
                vtkDoubleArray *ShearVectors2);      // Output

        void GenerateOutputDataSet(
                vtkDataSet *InputDataSet,
                vtkDataSet *OutputDataSet,           // Output
                vtkDoubleArray *MaxStrainValues,
                vtkDoubleArray *MaxStrainVectors,
                vtkDoubleArray *MinStrainValues,
                vtkDoubleArray *MinStrainVectors,
                vtkDoubleArray *ShearValues,
                vtkDoubleArray *ShearVectors1,
                vtkDoubleArray *ShearVectors2);

        // Member Data
        bool SmoothingStatus;
        unsigned int SmoothingKernelSize;
        bool ProgressStatus;
        
    private:
        Deformation(const Deformation &);   // Not implemented.
        void operator=(const Deformation &);   // Not implemented.
};

// ==================
// Deformation Helper
// ==================

namespace DeformationHelper
{
    double ** CreateMatrix(unsigned int NumberOfRows,unsigned int NumberOfColumns);
    void FreeMatrix(double **Matrix,unsigned int NumberOfColumns);
    void TransposeSquareMatrix(double **SquareMatrix,unsigned int MatrixSize);
}

// =========
// Templates
// =========

// Description:
// Defining an object that pairs each eigenvalue and it's correspondent
// eigenvector into a single object, including their original index
// in array. This object will be used for sorting eigenvalues.

template<class Type>
class EigenPair
{
    public:
        Type EigenValue;
        Type *EigenVector;
        unsigned int EigenIndexMap;
};

// Description:
// Overwrite the comparison operator based on accesnding eigenvalues.
// This is used for sorting eigenvalues in an accending manner.

template<class Type >
bool operator<(
        const EigenPair<Type> & EigenPairObject1,
        const EigenPair<Type> & EigenPairObject2)
{
    return EigenPairObject1.EigenValue < EigenPairObject2.EigenValue;
}

// Description:
// Given an eigenvalue and an array of eigenvector for that eigenvalue,
// this function creates the EigenPair object. This is used for sorting
// eigenvalues.

template<class Type>
EigenPair<Type> CreateEigenPairObject(
        const Type & InputEigenValue,
        Type *InputEigenVector,
        unsigned int EigenIndexMap)
{
    // Declare a local object
    EigenPair<Type> EigenPairObject;

    // Assign member data
    EigenPairObject.EigenValue = InputEigenValue;
    EigenPairObject.EigenVector = InputEigenVector;
    EigenPairObject.EigenIndexMap = EigenIndexMap;

    return EigenPairObject;
}

#endif
