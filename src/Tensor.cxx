/*
 * =====================================================================================
 *
 *       Filename:  Tensor.cxx
 *
 *    Description:  Container class for Tensors
 *
 *        Version:  1.0
 *        Created:  03/19/2014 11:37:58 AM
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

#include "Tensor.h"

// STL
#include <cstddef>   // for NULL
#include <cstdlib>   // for exit

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
Tensor::Tensor()
{
    // Member data
    this->Size = DIMENSION;
    this->Components = NULL;

    // Initialize Components
    this->InitializeTensor();
}

// Constructor with one argument
Tensor::Tensor(unsigned int InputSize)
{
    // Size
    this->Size = InputSize;

    // Components
    this->Components = NULL;

    // Initialize Components
    this->InitializeTensor();
}

// Constructor with Arguments
Tensor::Tensor(double **InputTensor,unsigned int InputSize)
{
    // Size
    this->Size = InputSize;

    // Components
    this->Components = new double*[this->Size];
    for(unsigned int RowIterator = 0; RowIterator < this->Size; RowIterator++)
    {
        this->Components[RowIterator] = new double[this->Size];

        for(unsigned int ColumnIterator = 0; ColumnIterator < this->Size; ColumnIterator++)
        {
            this->Components[RowIterator][ColumnIterator] = InputTensor[RowIterator][ColumnIterator];
        }
    }
}

// Copy Constructor
Tensor::Tensor(const Tensor &rhsTensor)
{
    this->Components = NULL;
    this->Clone(rhsTensor);
}

// Clone
void Tensor::Clone(const Tensor & rhsTensor)
{
    // Clear current Tensor
    this->DeleteTensor();

    // Copy new Tensor
    this->Size = rhsTensor.GetSize();
    this->InitializeTensor();

    // Deep copy values
    for(unsigned int RowIterator = 0; RowIterator < this->Size; RowIterator++)
    {
        for(unsigned int ColumnIterator = 0; ColumnIterator < this->Size; ColumnIterator++)
        {
            this->Components[RowIterator][ColumnIterator] = rhsTensor[RowIterator][ColumnIterator];
        }
    }
}

// ==========
// Destructor
// ==========

Tensor::~Tensor()
{
    // Components
    this->DeleteTensor();

    // Size
    this->Size = 0;
}

// =================
// Initialize Tensor
// =================

void Tensor::InitializeTensor()
{
    // Size
    if(this->Size == 0)
    {
        this->Size = DIMENSION;
    }

    // Assign new memory
    if(this->Components != NULL)
    {
        this->DeleteTensor();
    }

    this->Components = new double*[this->Size];

    // Assign columns of each row
    for(unsigned int RowIterator = 0; RowIterator < this->Size; RowIterator++)
    {
        this->Components[RowIterator] = new double[this->Size];

        // Set zero as default for all elements
        for(unsigned int ColumnIterator = 0; ColumnIterator < this->Size; ColumnIterator++)
        {
            this->Components[RowIterator][ColumnIterator] = 0;
        }
    }
}

// =============
// Delete Tensor
// =============

void Tensor::DeleteTensor()
{
    // Components
    if(this->Components != NULL)
    {
        // Delete Rows of each row
        for(unsigned int RowIterator = 0; RowIterator < this->Size; RowIterator++)
        {
            delete [] this->Components[RowIterator];
            this->Components[RowIterator] = NULL;
        }

        // Delete first column
        delete [] this->Components;
        this->Components = NULL;
    }
}

// ====================
// Accessors / Mutators
// ====================

// Get Components
double ** Tensor::GetComponents() const
{
    return this->Components;
}

// Get Component
double Tensor::GetComponent(unsigned int RowIndex,unsigned int ColumnIndex)
{
    return this->Components[RowIndex][ColumnIndex];
}

// Set Components
void Tensor::SetComponents(double **InputComponents)
{
    this->SetComponents(InputComponents,this->Size);
}

// Set Components (two arguments)
void Tensor::SetComponents(double **InputComponents,unsigned int Size)
{
    // Delete previous tensor and re-initialize
    this->DeleteTensor();
    this->InitializeTensor();

    // Deep copy from input
    for(unsigned int RowIterator = 0; RowIterator < this->Size; RowIterator++)
    {
        for(unsigned int ColumnIterator = 0; ColumnIterator < this->Size; ColumnIterator++)
        {
            this->Components[RowIterator][ColumnIterator] = InputComponents[RowIterator][ColumnIterator];
        }
    }
}

// Set Component
void Tensor::SetComponent(unsigned int RowIndex,unsigned int ColumnIndex,double ComponentValue)
{
    // Check if Tensor was initialized
    if(this->Components == NULL)
    {
        this->InitializeTensor();
    }

    // Set value
    this->Components[RowIndex][ColumnIndex] = ComponentValue;
}

// Get Size
unsigned int Tensor::GetSize() const
{
    return this->Size;
}

// Set Size
void Tensor::SetSize(unsigned int InputSize)
{
    this->Size = InputSize;
}


// ==========
// Check Size
// ==========

void Tensor::CheckSize(const Tensor & Tensor1,const Tensor & Tensor2)
{
    if(Tensor1.GetSize() != Tensor2.GetSize())
    {
        std::cerr << "Tensors sizes do not match." << std::endl;
        exit(0);
    }
}

// =========
// Operators
// =========

// +
Tensor operator+(const Tensor & Tensor1,const Tensor & Tensor2)
{
    // Check sizes
    Tensor::CheckSize(Tensor1,Tensor2);

    // Output Tensor
    Tensor OutputTensor;

    // Update output
    for(unsigned int RowIterator = 0; RowIterator < Tensor1.GetSize(); RowIterator++)
    {
        for(unsigned int ColumnIterator = 0; ColumnIterator < Tensor1.GetSize(); ColumnIterator++)
        {
            double OutputTensorComponent =
                Tensor1[RowIterator][ColumnIterator] + 
                Tensor2[RowIterator][ColumnIterator];
            OutputTensor.SetComponent(RowIterator,ColumnIterator,OutputTensorComponent);
        }
    }

    return OutputTensor;
}

// -
Tensor operator-(const Tensor & Tensor1,const Tensor & Tensor2)
{
    // Check sizes
    Tensor::CheckSize(Tensor1,Tensor2);

    // Output Tensor
    Tensor OutputTensor;

    // Update output
    for(unsigned int RowIterator = 0; RowIterator < Tensor1.GetSize(); RowIterator++)
    {
        for(unsigned int ColumnIterator = 0; ColumnIterator < Tensor1.GetSize(); ColumnIterator++)
        {
            double OutputTensorComponent =
                Tensor1[RowIterator][ColumnIterator] -
                Tensor2[RowIterator][ColumnIterator];
            OutputTensor.SetComponent(RowIterator,ColumnIterator,OutputTensorComponent);
        }
    }

    return OutputTensor;
}

// Tensor-Vector *
Vector Tensor::operator*(const Vector & rhsVector)
{
    Vector OutputVector;
    OutputVector.SetLength(this->Size);

    // Matrix-Vector product
    for(unsigned int RowIterator = 0; RowIterator < this->Size; RowIterator++)
    {
        double Sum = 0;
        for(unsigned int DummyIterator = 0; DummyIterator < this->Size; DummyIterator++)
        {
            Sum += this->Components[RowIterator][DummyIterator] * rhsVector[DummyIterator];
        }

        OutputVector.SetComponent(RowIterator,Sum);
    }
    
    return OutputVector;
}

// ==
bool operator==(const Tensor & Tensor1,const Tensor & Tensor2)
{
    bool Status = true;

    // Check size
    if(Tensor1.GetSize()!= Tensor2.GetSize())
    {
        Status = false;
        return Status;
    }

    // Check components
    for(unsigned int RowIterator = 0; RowIterator < Tensor1.GetSize(); RowIterator++)
    {
        for(unsigned int ColumnIterator = 0; ColumnIterator < Tensor1.GetSize(); ColumnIterator++)
        {
            double Component1 = Tensor1[RowIterator][ColumnIterator];
            double Component2 = Tensor2[RowIterator][ColumnIterator];

            if(abs(Component1 - Component2) > 1e-15)
            {
                Status = false;
                return Status;
            }
        }
    }

    return Status;
}

// !=
bool operator!=(const Tensor & Tensor1,const Tensor & Tensor2)
{
    return !(Tensor1 == Tensor2);
}

// []
double* Tensor::operator[](unsigned int RowIndex) const
{
    return this->Components[RowIndex];
}

// =
Tensor & Tensor::operator=(const Tensor & rhsTensor)
{
    if(this == &rhsTensor)
    {
        return *this;
    }

    // Update this
    this->Clone(rhsTensor);

    return *this;
}

// +=
Tensor & Tensor::operator+=(const Tensor & rhsTensor)
{
    // Operation on a temporary object
    Tensor TempTensor = *this + rhsTensor;

    // Update "this" object
    this->Clone(TempTensor);

    // Output
    return *this;
}

// -=
Tensor & Tensor::operator-=(const Tensor & rhsTensor)
{
    // Operation on a temporary object
    Tensor TempTensor = *this - rhsTensor;

    // Update "this" object
    this->Clone(TempTensor);

    // Output
    return *this;
}

// <<
std::ostream & operator<<(std::ostream & os,const Tensor & rhsTensor)
{
    for(unsigned int RowIterator = 0; RowIterator < rhsTensor.GetSize(); RowIterator++)
    {
        for(unsigned int ColumnIterator = 0; ColumnIterator < rhsTensor.GetSize(); ColumnIterator++)
        {
            // Print component
            os << rhsTensor[RowIterator][ColumnIterator];

            // Add comma
            if(ColumnIterator < rhsTensor.GetSize()-1)
            {
                os << ", ";
            }
        }

        // New line
        os << std::endl;
    }

    return os;
}

// lhs Vector *
Vector operator*(const Vector & lhsVector,const Tensor & rhsTensor)
{
    // Check size
    if(lhsVector.GetLength() != rhsTensor.GetSize())
    {
        std::cerr << "Size dimension does not match." << std::endl;
        exit(0);
    }

    Vector OutputVector;
    OutputVector.SetLength(rhsTensor.GetSize());

    // Vector-Matrix product
    for(unsigned int ColumnIterator = 0; ColumnIterator < rhsTensor.GetSize(); ColumnIterator++)
    {
        double Sum = 0;
        for(unsigned int DummyIterator = 0; DummyIterator < rhsTensor.GetSize(); DummyIterator++)
        {
            Sum += lhsVector[DummyIterator] * rhsTensor[DummyIterator][ColumnIterator];
        }

        OutputVector.SetComponent(ColumnIterator,Sum);
    }
    
    return OutputVector;
}
