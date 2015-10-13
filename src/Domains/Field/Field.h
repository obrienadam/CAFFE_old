/**
 * @file    Field.h
 * @author  Adam O'Brien <obrienadam89@gmail.com>
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * https://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 *
 * This file contains the interface for class Field, which is used for
 * storing various types of 3-Dimensional fields and their boundaries.
 */

#ifndef FIELD_H
#define FIELD_H

#include <string>
#include <iostream>
#include <functional>

#include "Array3D.h"
#include "HexaFvmMesh.h"

template <class T>
class Field
{
public:

    enum FieldType{CONSERVED, AUXILLARY};
    enum BoundaryPatch{FIXED, ZERO_GRADIENT, EMPTY};

    Field(const HexaFvmMesh &mesh, FieldType type, std::string name);
    Field(const Field& other);

    Field<T>& operator=(const Field<T>& other);

    int sizeI() const { return fieldData_.sizeI(); }
    int sizeJ() const { return fieldData_.sizeJ(); }
    int sizeK() const { return fieldData_.sizeK(); }
    int size() const { return fieldData_.size(); }

    //- Access
    T& operator()(int i, int j, int k);
    const T& operator()(int i, int j, int k) const;

    T& operator()(int i, int j, int k, int faceNo);
    const T& operator()(int i, int j, int k, int faceNo) const;

    T& operator()(int k);
    const T& operator()(int k) const;

    const HexaFvmMesh& getMesh() const { return mesh_; }

    void getSubfield(int iLower, int iUpper,
                     int jLower, int jUpper,
                     int kLower, int kUpper,
                     Array3D<T> &subField) const;

    void setValues(int iLower, int iUpper,
                   int jLower, int jUpper,
                   int kLower, int kUpper,
                   T value);

    void setValue(T value);
    void setInitialCondition(const std::function<T (Point3D)> &icFunction);

    void setFixedBoundaryPatches(const T *refValues);
    void setFixedBoundaryPatches(const T &refValue);
    void setAll(T value);

    void setEastFacesFromPatch();
    void setWestFacesFromPatch();
    void setNorthFacesFromPatch();
    void setSouthFacesFromPatch();
    void setTopFacesFromPatch();
    void setBottomFacesFromPatch();

    void setZeroGradientBoundaryEast();
    void setZeroGradientBoundaryWest();
    void setZeroGradientBoundaryNorth();
    void setZeroGradientBoundarySouth();
    void setZeroGradientBoundaryTop();
    void setZeroGradientBoundaryBottom();

    T& faceE(int i, int j, int k){ return facesI_(i + 1, j, k); }
    T& faceW(int i, int j, int k){ return facesI_(i, j, k); }
    T& faceN(int i, int j, int k){ return facesJ_(i, j + 1, k); }
    T& faceS(int i, int j, int k){ return facesJ_(i, j, k); }
    T& faceT(int i, int j, int k){ return facesK_(i, j, k + 1); }
    T& faceB(int i, int j, int k){ return facesK_(i, j, k); }
    T& face(int i, int j, int k, int faceNo);

    T& faceI(int i, int j, int k){ return facesI_(i, j, k); }
    T& faceJ(int i, int j, int k){ return facesJ_(i, j, k); }
    T& faceK(int i, int j, int k){ return facesK_(i, j, k); }

    const T& faceE(int i, int j, int k) const { return facesI_(i + 1, j, k); }
    const T& faceW(int i, int j, int k) const { return facesI_(i, j, k); }
    const T& faceN(int i, int j, int k) const { return facesJ_(i, j + 1, k); }
    const T& faceS(int i, int j, int k) const { return facesJ_(i, j, k); }
    const T& faceT(int i, int j, int k) const { return facesK_(i, j, k + 1); }
    const T& faceB(int i, int j, int k) const { return facesK_(i, j, k); }
    const T& face(int i, int j, int k, int faceNo) const;

    const T& faceI(int i, int j, int k) const { return facesI_(i, j, k); }
    const T& faceJ(int i, int j, int k) const { return facesJ_(i, j, k); }
    const T& faceK(int i, int j, int k) const { return facesK_(i, j, k); }

    T maxNeighbour(int i, int j, int k);
    T minNeighbour(int i, int j, int k);

    const Array3D<T>* cellData() const { return &fieldData_; }
    void print();

    Array3D<T> eastBoundaryPatch;
    Array3D<T> westBoundaryPatch;
    Array3D<T> northBoundaryPatch;
    Array3D<T> southBoundaryPatch;
    Array3D<T> topBoundaryPatch;
    Array3D<T> bottomBoundaryPatch;

    std::string name;

protected:

    void allocate(int nCellsI, int nCellsJ, int nCellsK);

    //- Main data
    Array3D<T> fieldData_;

    //- Reference to the mesh, required for constructing a field
    const HexaFvmMesh &mesh_;
    FieldType type_;

    //- Face nodes for primitive and conserved fields
    Array3D<T> facesI_;
    Array3D<T> facesJ_;
    Array3D<T> facesK_;
};

#include "Field.tpp"

#endif
