/*
 * =====================================================================================
 *
 *       Filename:  Vector.cxx
 *
 *    Description:  Vector class for general purpose vectors
 *
 *        Version:  1.0
 *        Created:  02/21/2014 11:36:18 AM
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

#include "Vector.h"

// STL
#include <cstddef>  // for NULL
#include <cstdlib>  // for exit
#include <cmath>    // for sqrt, pow

// ======
// Macros
// ======

#ifndef DIMENSION
#define DIMENSION 3
#endif

// ===========
// Constructor
// ===========

// Constructor
Vector::Vector()
{
    // Member Data
    this->Length = DIMENSION;
    this->Components = NULL;

    // Initialize Componments
    this->InitializeVector();
}

// Constructor with one argument
Vector::Vector(unsigned int ArrayLength)
{
    // Length
    this->Length = ArrayLength;

    // Components
    this->Components = NULL;

    // Initialize Components
    this->InitializeVector();
}

// Constructor with arguments
Vector::Vector(double *InputArray,unsigned int ArrayLength)
{
    // Length
    this->Length = ArrayLength;

    // Components
    this->Components = new double[this->Length];
    for(unsigned int Iterator = 0; Iterator < this->Length; Iterator++)
    {
        this->Components[Iterator] = InputArray[Iterator];
    }
}

// Copy Constructor
Vector::Vector(const Vector &rhsVector)
{
    this->Components = NULL;
    this->Clone(rhsVector);
}

// Clone
void Vector::Clone(const Vector & rhsVector)
{
    // Clear current vector
    this->DeleteVector();

    // Copy new vector
    this->Length = rhsVector.GetLength();
    this->InitializeVector();

    // Deep copy values
    for(unsigned int Index = 0; Index < this->Length; Index++)
    {
        this->Components[Index] = rhsVector.GetComponent(Index);
    }
}

// ==========
// Destructor
// ==========

Vector::~Vector()
{
    // Components
    this->DeleteVector();

    // Length
    this->Length = 0;
}

// =================
// Initialize Vector
// =================

void Vector::InitializeVector()
{
    // Length
    if(this->Length == 0)
    {
        this->Length = DIMENSION;
    }

    // Assign new memory
    if(this->Components != NULL)
    {
        this->DeleteVector();
    }

    this->Components = new double[this->Length];

    // Set default values to zero
    for(unsigned int Index = 0; Index < this->Length; Index++)
    {
        this->Components[Index] = 0;
    }
}

// =============
// Delete Vector
// =============

void Vector::DeleteVector()
{
    // Components
    if(this->Components != NULL)
    {
        delete [] this->Components;
        this->Components = NULL;
    }
}

// ====================
// Accessors / Mutators
// ====================

// Get Components
double * Vector::GetComponents() const
{
    return this->Components;
}

// Get Component
double Vector::GetComponent(unsigned int Index) const
{
    return this->Components[Index];
}

// Set Components (one argument)
void Vector::SetComponents(double *InputComponents)
{
    this->SetComponents(InputComponents,this->Length);
}

// Set Components (two arguments)
void Vector::SetComponents(double *InputComponents,unsigned int InputLength)
{
    // Delete previous vector and re-initialize
    this->DeleteVector();
    this->InitializeVector();

    // Deep copy from input
    for(unsigned int Index = 0; Index < InputLength; Index++)
    {
        this->Components[Index] = InputComponents[Index];
    }
}

// Set Component
void Vector::SetComponent(unsigned int Index,double ComponentValue)
{
    // Check is vector was initialized
    if(this->Components == NULL)
    {
        this->InitializeVector();
    }

    // Set value
    this->Components[Index] = ComponentValue;
}

// Get Length
unsigned int Vector::GetLength() const
{
    return this->Length;
}

// Set Length
void Vector::SetLength(unsigned int InputLength)
{
    this->Length = InputLength;
}

// =============
// Cross Product
// =============

Vector Vector::CrossProduct(const Vector & Vector1, const Vector & Vector2)
{
    // Check Length
    Vector::CheckLength(Vector1,Vector2);

    // Length
    unsigned int Length = Vector1.GetLength();

    // This method works for vector length 3
    if(Length != 3)
    {
        std::cerr << "Cross product workds for 3D vectors." << std::endl;
    }

    // Output
    Vector OutputVector;
    OutputVector.SetLength(Length);

    OutputVector.SetComponent(0,Vector1[1]*Vector2[2] - Vector1[2]*Vector2[1]);
    OutputVector.SetComponent(1,Vector1[2]*Vector2[0] - Vector1[0]*Vector2[2]);
    OutputVector.SetComponent(2,Vector1[0]*Vector2[1] - Vector1[1]*Vector2[0]);

    return OutputVector;
}

// ===========
// Dot Product
// ===========

double Vector::DotProduct(const Vector & Vector1,const Vector & Vector2)
{
    // Check Lengths
    Vector::CheckLength(Vector1,Vector2);

    // Sum over component products
    double Sum = 0;
    for(unsigned int Index = 0; Index < Vector1.GetLength(); Index++)
    {
        Sum += Vector1.GetComponent(Index) * Vector2.GetComponent(Index);
    }

    return Sum;
}

// ========
// Get Norm
// ========

double Vector::GetNorm() const
{
    double Sum = 0;
    for(unsigned int Index = 0; Index < this->Length; Index++)
    {
        Sum += pow(this->Components[Index],2);
    }

    return sqrt(Sum);
}

// ========
// Normaliz
// ========

// Note: this method over-writes the "this" object.
double Vector::Normalize()
{
    double VectorNorm = this->GetNorm();
    Vector TempVector = (*this)/VectorNorm;
    this->Clone(TempVector);

    return VectorNorm;
}

// ========
// Opposite
// ========

// Note: this method over-writes the "this" object.
Vector & Vector::Opposite()
{
    for(unsigned int Index = 0; Index < this->Length; Index++)
    {
        this->Components[Index] = -this->Components[Index];
    }

    return *this;
}

// ============
// Check Length
// ============

void Vector::CheckLength(const Vector & Vector1,const Vector & Vector2)
{
    if(Vector1.GetLength() != Vector2.GetLength())
    {
        std::cerr << "Vectors lengths do not match." << std::endl;
        exit(0);
    }
}

// ===============================
// Convert Vectors Array To Matrix
// ===============================

// Description:
// Matrix allocation is done inside function.

void Vector::ConvertVectorsArrayToMatrix(
        Vector *VectorsArray,
        unsigned int NumberOfVectors,
        double **&Matrix)
{
    // Get Number of Components for each vector
    unsigned int NumberOfComponents = VectorsArray[0].GetLength();

    //  Allocate Matrix
    Matrix = new double*[NumberOfVectors];
    for(unsigned int RowIterator = 0; RowIterator < NumberOfVectors; RowIterator++)
    {
        Matrix[RowIterator] = new double[NumberOfComponents];
    }

    // Deep copy
    for(unsigned int VectorIterator = 0; VectorIterator < NumberOfVectors; VectorIterator++)
    {
        for(unsigned int ComponentIterator = 0; ComponentIterator < NumberOfComponents; ComponentIterator++)
        {
            Matrix[VectorIterator][ComponentIterator] = VectorsArray[VectorIterator].GetComponent(ComponentIterator);
        }
    }
}

// =========
// Operators
// =========

// +
Vector operator+(const Vector & Vector1,const Vector & Vector2)
{
    // Check Lengths
    Vector::CheckLength(Vector1,Vector2);

    // Output Vector
    Vector OutputVector;

    // Update output
    for(unsigned int Index = 0; Index < Vector1.GetLength(); Index++)
    {
        OutputVector.SetComponent(Index,Vector1[Index] + Vector2[Index]);
    }

    return OutputVector;
}

// -
Vector operator-(const Vector & Vector1,const Vector & Vector2)
{
    // Check Lengths
    Vector::CheckLength(Vector1,Vector2);

    // Output Vector
    Vector OutputVector;

    // Update output
    for(unsigned int Index = 0; Index < Vector1.GetLength(); Index++)
    {
        OutputVector.SetComponent(Index,Vector1[Index] - Vector2[Index]);
    }

    return OutputVector;
}

// * (Dot Product)
double operator*(const Vector & Vector1,const Vector & Vector2)
{
    return Vector::DotProduct(Vector1,Vector2);
}

// ==
bool operator==(const Vector & Vector1,const Vector & Vector2)
{
    bool Status = true;

    // Check lengths
    if(Vector1.GetLength() != Vector2.GetLength())
    {
        Status = false;
        return Status;
    }

    // Check components
    for(unsigned int Index = 0; Index < Vector1.GetLength(); Index++)
    {

        if(abs(Vector1[Index] - Vector2[Index]) > 1e-15)
        {
            Status = false;
            return Status;
        }
    }

    return Status;
}

// !=
bool operator!=(const Vector & Vector1,const Vector & Vector2)
{
    return !(Vector1 == Vector2);
}

// []
double Vector::operator[](unsigned int Index) const
{
    return this->GetComponent(Index);
}

// =
Vector & Vector::operator=(const Vector & rhsVector)
{
    if(this == &rhsVector)
    {
        return *this;
    }

    // Update this
    this->Clone(rhsVector);

    return *this;
}

// +=
Vector & Vector::operator+=(const Vector & rhsVector)
{
    // Operation on a temporary object
    Vector TempVector = *this + rhsVector;

    // Update "this" object
    this->Clone(TempVector);

    // Output
    return *this;
}

// -=
Vector & Vector::operator-=(const Vector & rhsVector)
{
    // Operation on a temporary object
    Vector TempVector = *this - rhsVector;

    // Update "this" object
    this->Clone(TempVector);

    // Output
    return *this;
}

// <<
std::ostream & operator<<(std::ostream & os,const Vector & rhsVector)
{
    unsigned int rhsVectorLength = rhsVector.GetLength();

    for(unsigned int Index = 0; Index < rhsVectorLength; Index++)
    {
        // Print component value
        os << rhsVector[Index];

        // Adding comma
        if(Index < rhsVectorLength - 1)
        {
            os << ", ";
        }
    }

    return os;
}
