/*
 * =====================================================================================
 *
 *       Filename:  Tensor.h
 *
 *    Description:  Container class for Tensors
 *
 *        Version:  1.0
 *        Created:  03/19/2014 11:38:20 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University Of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __Tensor_h
#define __Tensor_h

// =====================
// Foreward Declarations
// =====================

#include "Vector.h"

// ============
// Tensor Class
// ============

class Tensor
{
    public:
        // Constructor . Destructor
        Tensor();
        Tensor(unsigned int InputSize);
        Tensor(double **InputTensor,unsigned int InputSize);
        Tensor(const Tensor & rhsTensor);
        void Clone(const Tensor & rhsTensor);
        ~Tensor();

        // Accessors / Mutators
        double **GetComponents() const;
        double GetComponent(unsigned int RowIndex,unsigned int ColumnIndex);
        void SetComponents(double **InputComponents);
        void SetComponents(double **InputComponents,unsigned int Size);
        void SetComponent(unsigned int RowIndex,unsigned int CoumnIndex,double ComponentValue);
       
        unsigned int GetSize() const;
        void SetSize(unsigned int InputSize);

        // Public methods
        static void CheckSize(const Tensor & Tensor1,const Tensor & Tensor2);
       
        // RHS Operators
        friend Tensor operator+(const Tensor & Tensor1,const Tensor & Tensor2);
        friend Tensor operator-(const Tensor & Tensor1,const Tensor & Tensor2);
        template<class Type> Tensor operator*(const Type rhsScalar);   // Tensor-Scalar product
        Vector operator*(const Vector &rhsVector);                     // Tesnor-Vector product
        friend bool operator==(const Tensor & Tensor1,const Tensor & Tensor2);
        friend bool operator!=(const Tensor & Tensor1,const Tensor & Tensor2);
        double * operator[](unsigned int RowIndex) const;
        Tensor & operator=(const Tensor & rhsTensor);
        Tensor & operator+=(const Tensor & rhsTensor);
        Tensor & operator-=(const Tensor & rhsTensor);

        // LHS Operators
        friend std::ostream & operator<<(std::ostream & os,const Tensor & rhsTensor);
        friend Vector operator*(const Vector & lhsVector,const Tensor & rhsTensor);
        template<class Type> friend Tensor operator*(const Type lhsScalar,const Tensor & rhsTensor);
        
    protected:
        void InitializeTensor();
        void DeleteTensor();

    private:
        // member data
        double **Components;
        unsigned int Size;
};

// ==================
// Template Operators
// ==================

// rhs scalar *
template<class Type> Tensor Tensor::operator*(const Type rhsScalar)
{
    Tensor OutputTensor(*this);

    for(unsigned int RowIterator = 0; RowIterator < this->Size; RowIterator++)
    {
        for(unsigned int ColumnIterator = 0; ColumnIterator < this->Size; ColumnIterator++)
        {
            double Value = this->Components[RowIterator][ColumnIterator] * rhsScalar;
            OutputTensor.SetComponent(RowIterator,ColumnIterator,Value);
        }
    }

    return OutputTensor;
}

// lhs scalar *
template<class Type> Tensor operator*(const Type lhsScalar,const Tensor & rhsTensor)
{
    return Tensor(rhsTensor) * lhsScalar;
}

#endif
