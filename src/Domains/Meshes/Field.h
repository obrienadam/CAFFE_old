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

#include "Array3D.h"
#include "HexaFvmMesh.h"

template <class T>
class Field : public Array3D<T>
{
public:

    enum FieldType{CONSERVED, AUXILLARY};
    enum BoundaryPatch{FIXED, ZERO_GRADIENT, EMPTY};

    Field(const HexaFvmMesh &mesh, FieldType type, std::string name);
    Field(const Field& other);

    Field<T>& operator=(const Field<T>& other);

    //- Access
    T& operator()(int i, int j, int k);
    const T& operator()(int i, int j, int k) const;

    T& operator()(int i, int j, int k, int faceNo);
    const T& operator()(int i, int j, int k, int faceNo) const;

    T& operator()(int k);
    const T& operator()(int k) const;

    const HexaFvmMesh& getMesh() const { return mesh_; }

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

    void setEastBoundary(const std::string &boundaryType, T boundaryValue);
    void setWestBoundary(const std::string &boundaryType, T boundaryValue);
    void setNorthBoundary(const std::string &BoundaryType, T boundaryValue);
    void setSouthBoundary(const std::string &boundaryType, T boundaryValue);
    void setTopBoundary(const std::string &boundaryType, T boundaryValue);
    void setBottomBoundary(const std::string &boundaryType, T boundaryValue);

    void setBoundaryFields();

    void setEastBoundaryField();
    void setWestBoundaryField();
    void setNorthBoundaryField();
    void setSouthBoundaryField();
    void setTopBoundaryField();
    void setBottomBoundaryField();

    BoundaryPatch getEastBoundaryPatch(){ return eastBoundaryPatchId_; }
    BoundaryPatch getWestBoundaryPatch(){ return westBoundaryPatchId_; }
    BoundaryPatch getNorthBoundaryPatch(){ return northBoundaryPatchId_; }
    BoundaryPatch getSouthBoundaryPatch(){ return southBoundaryPatchId_; }
    BoundaryPatch getTopBoundaryPatch(){ return topBoundaryPatchId_; }
    BoundaryPatch getBottomBoundaryPatch(){ return bottomBoundaryPatchId_; }

    Array3D<T> eastBoundaryPatch_;
    Array3D<T> westBoundaryPatch_;
    Array3D<T> northBoundaryPatch_;
    Array3D<T> southBoundaryPatch_;
    Array3D<T> topBoundaryPatch_;
    Array3D<T> bottomBoundaryPatch_;

    void setImplicitBoundaryCoeffs(int i, int j, int k, double *a, T &b);

    const Array3D<T>* cellData() const { return this; }

    /**
     * @brief Print the field to the console, for debugging purposes.
     */
    void print();

    std::string name;

protected:

    void allocate(int nCellsI, int nCellsJ, int nCellsK);

    //- Reference to the mesh, required for constructing a field
    const HexaFvmMesh &mesh_;
    FieldType type_;

    //- Face nodes for primitive and conserved fields
    Array3D<T> facesI_;
    Array3D<T> facesJ_;
    Array3D<T> facesK_;

    //- Boundary patches
    BoundaryPatch eastBoundaryPatchId_;
    BoundaryPatch westBoundaryPatchId_;
    BoundaryPatch northBoundaryPatchId_;
    BoundaryPatch southBoundaryPatchId_;
    BoundaryPatch topBoundaryPatchId_;
    BoundaryPatch bottomBoundaryPatchId_;
};

#include "FieldI.h"

#endif
