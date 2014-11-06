/*
 * =====================================================================================
 *
 *       Filename:  Vector.h
 *
 *    Description:  Vector class for general purpose vectors
 *
 *        Version:  1.0
 *        Created:  02/21/2014 11:40:11 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University Of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __Vector_h
#define __Vector_h

// ======
// Macros
// ======

// Debugging
#ifndef HERE
#define HERE std::cout << "\033[0;44m" << "DEBUG: " << __FILE__ << " at " << __LINE__ << "\033[0m" << std::endl;
#endif

// ======================
// Foreward Declaratrions
// ======================

#include <iostream>

// ============
// Vector Class
// ============

class Vector
{
    public:
        // Constructor / Destructor
        Vector();
        Vector(unsigned int ArrayLength);
        Vector(double *InputArray,unsigned int ArrayLength);
        Vector(const Vector & rhsVector);
        void Clone(const Vector & rhsVector);
        ~Vector();

        // Accessors / Mutators
        double * GetComponents() const;
        double GetComponent(unsigned int Index) const;
        void SetComponents(double *InputComponents);
        void SetComponents(double *InputComponents, unsigned int InputLength);
        void SetComponent(unsigned int Index,double ComponentValue);

        unsigned int GetLength() const;
        void SetLength(unsigned int InputLength);

        // Methods
        static double DotProduct(const Vector &Vector1,const Vector &Vector2);
        static Vector CrossProduct(const Vector &Vector1,const Vector &Vector2);
        double GetNorm() const;
        double Normalize();
        Vector & Opposite();
        static void CheckLength(const Vector & Vector1,const Vector & Vector2);
        static void ConvertVectorsArrayToMatrix(
                Vector *VectorsArray,unsigned int NumberOfVectors,double **&Matrix);

        // RHS Operators
        friend Vector operator+(const Vector & Vector1,const Vector & Vector2);
        friend Vector operator-(const Vector & Vector1,const Vector & Vector2);
        friend double operator*(const Vector & Vector1,const Vector & Vector2);   // Dot Product
        friend bool operator==(const Vector & Vector1,const Vector & Vector2);
        friend bool operator!=(const Vector & Vector1,const Vector & Vector2);
        double operator[](unsigned int Index) const;
        Vector & operator=(const Vector & rhsVector);
        Vector & operator+=(const Vector & rhsVector);
        Vector & operator-=(const Vector & rhsVector);

        // RHS Template Operators
        template <class Type> Vector operator*(const Type rhsScalar) const;
        template <class Type>  Vector operator/(const Type rhsScalar) const;

        // LHS Template Operators
        friend std::ostream & operator<<(std::ostream & os,const Vector & rhsVector);
        template<class Type> friend Vector operator*(const Type & lhsScalar,const Vector & rhsVector);

    protected:
        void InitializeVector();
        void DeleteVector();

    private:
        double *Components;
        unsigned int Length;
};

// ==================
// Template Operators
// ==================

// rhs *
template <class Type> Vector Vector::operator*(const Type rhsScalar) const
{
    Vector OutputVector(*this);

    for(unsigned int Index = 0; Index < this->Length; Index++)
    {
        OutputVector.SetComponent(Index,this->GetComponent(Index)*rhsScalar);
    }

    return OutputVector;
}

// rhs /
template <class Type> Vector Vector::operator/(const Type rhsScalar) const
{
    return this->operator*(1/rhsScalar);
}

// lhs *
template <class Type> Vector operator*(const Type & lhsScalar,const Vector & rhsVector)
{
    return Vector(rhsVector) *  lhsScalar;
}

#endif
