/*
 * =====================================================================================
 *
 *       Filename:  PointDataType.cxx
 *
 *    Description:  A container class for Iso-Surface Points
 *
 *        Version:  1.0
 *        Created:  04/02/2014 03:23:17 PM
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

#include "PointDataType.h"
#include "Vector.h"

// ===========
// Constructor
// ===========

PointDataType::PointDataType():
    OrganizedPointId(0)
{
}

// ================
// Copy Constructor
// ================

PointDataType::PointDataType(const PointDataType & rhs)
{
    this->Coordinates = rhs.GetCoordinates();
    this->DeformationVector = rhs.GetDeformationVector();
    this->DeformationValue = rhs.GetDeformationValue();
    this->DeformationValueCurvature = rhs.GetDeformationValueCurvature();
    this->OrganizedPointId = rhs.GetOrganizedPointId();
}

// ==========
// Destructor
// ==========

PointDataType::~PointDataType()
{
}

// ====================
// Accessors / Mutators
// ====================

// Get Coordinates
Vector PointDataType::GetCoordinates() const
{
    return this->Coordinates;
}

// Set Coordinates
void PointDataType::SetCoordinates(Vector InputCoordinates)
{
    this->Coordinates = InputCoordinates;
}

// Get Deformation Vector
Vector PointDataType::GetDeformationVector() const
{
    return this->DeformationVector;
}

// Set Deformation Vector
void PointDataType::SetDeformationVector(Vector InputDeformationVector)
{
    this->DeformationVector = InputDeformationVector;
}

// Get Deformation Value
double PointDataType::GetDeformationValue() const
{
    return this->DeformationValue;
}

// Set Deformation Value
void PointDataType::SetDeformationValue(double InputDeformationValue)
{
    this->DeformationValue = InputDeformationValue;
}

// Get Deformation Value Curvature
double PointDataType::GetDeformationValueCurvature() const
{
    return this->DeformationValueCurvature;
}

// Set Deformation Value Curvature
void PointDataType::SetDeformationValueCurvature(double InputDeformationValueCurvature)
{
    this->DeformationValueCurvature = InputDeformationValueCurvature;
}

// Get Organized Point Id
unsigned int PointDataType::GetOrganizedPointId() const
{
    return this->OrganizedPointId;
}

// Set Organized Point Id
void PointDataType::SetOrganizedPointId(unsigned int InputOrganizedPointId)
{
    this->OrganizedPointId = InputOrganizedPointId;
}
