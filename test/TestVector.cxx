/*
 * =====================================================================================
 *
 *       Filename:  TestVector.cxx
 *
 *    Description:  Test for Vector.cxx
 *
 *        Version:  1.0
 *        Created:  02/26/2014 03:13:36 PM
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

#include "Vector.h"
#include <iostream>

// ====
// Main
// ====

int main(int argc,char *argv[])
{
    double Array1[3] = {2.3,3.4,4.5};
    double Array2[3] = {1,2,4};

    // Define a V1
    Vector V1;
    V1.SetComponents(Array1,3);
    std::cout << "V1: " << V1 << std::endl;
    V1.SetComponent(1,0.2);
    std::cout << "V1: " << V1 << std::endl;

    // Define V2
    Vector V2;
    V2.SetComponents(Array2,3);
    std::cout << "V2: " << V2 << std::endl; 

    // Arithmetic operations
    Vector V3 = V1 + V2;
    std::cout << "V3: = V1 + V2: " << V3 << std::endl;
    std::cout << "V1 - V2: " << V1 - V2 << std::endl;
    std::cout << "V1 == V2: " << static_cast<int>(V1 == V2) << std::endl;
    std::cout << "V1 != V2: " << static_cast<int>(V1 != V2) << std::endl;
    std::cout << "Dot Product: V1.V2: " << V1 * V2 << std::endl;
    std::cout << "Cross Product: V1^V2: " << Vector::CrossProduct(V1,V2) << std::endl;
    V3 = V1;
    std::cout << "V3 = V1: " << V3 << std::endl;
    V3 += V2;
    std::cout << "V3 += V2: " << V3 << std::endl;
    V3 -= V2;
    std::cout << "V3 -= V2: " << V3 << std::endl;
    std::cout << "V3 * 2.5: " << V3 * 2.5 << std::endl;
    std::cout << "2.5 * V3: " << 2.5 * V3<< std::endl;
    std::cout << "V3 / 2.5: " << V3 / 2.5 << std::endl;
    std::cout << "V3.Opposite(): " << V3.Opposite() << std::endl;
    std::cout << "V3 Norm: " << V3.GetNorm() << std::endl;
    std::cout << "Normalize V3 norm: " << V3.Normalize() << ". V3: " << V3 << std::endl;

    // Convert Array of vectors to Matrix
    Vector *VectorsArray = new Vector[2];
    VectorsArray[0] = V1;
    VectorsArray[1] = V2;
    double **Matrix;
    // double **Matrix = new double*[2];
    // for(int i = 0; i < 2; i++)
    // {
    //     Matrix[i] = new double[3];
    // }
    Vector::ConvertVectorsArrayToMatrix(VectorsArray,2,Matrix);

    for(int i = 0; i < 2; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            std::cout << "Matrix[" << i <<"][" << j << "] = " << Matrix[i][j] << ", ";
        }
        std::cout << std::endl;
    }

    // Delete Matrix
    for(int i=0;i<2;i++)
    {
        delete [] Matrix[i];
    }
    delete [] Matrix;
    Matrix = NULL;

    return 0;
}
