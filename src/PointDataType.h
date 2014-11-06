/*
 * =====================================================================================
 *
 *       Filename:  PointDataType.h
 *
 *    Description:  A container class for Iso-Surface Points
 *
 *        Version:  1.0
 *        Created:  04/02/2014 03:23:57 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University Of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __PointDataType_h
#define __PointDataType_h

// ====================
// Forward Declarations
// ====================

#include "Vector.h"

// =====================
// Point Data Type Class
// =====================


// Container class for points of Iso-surface manifold
class PointDataType
{
    public:
        // Constructor
        PointDataType();
        PointDataType(const PointDataType & rhs);
        ~PointDataType();

        // Accessors
        Vector GetCoordinates() const;
        void SetCoordinates(Vector InputCoordinates);

        Vector GetDeformationVector() const;
        void SetDeformationVector(Vector InputDeformationVector);

        double GetDeformationValue() const;
        void SetDeformationValue(double InputDeformationValue);

        double GetDeformationValueCurvature() const;
        void SetDeformationValueCurvature(double InputDeformationValueCurvature);

        unsigned int GetOrganizedPointId() const;
        void SetOrganizedPointId(unsigned int InpiutOrganizedPointId);

    private:
        // Member Data
        Vector Coordinates;
        Vector DeformationVector;
        double DeformationValue;
        double DeformationValueCurvature;
        unsigned int OrganizedPointId;
};

#endif
