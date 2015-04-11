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

#include "Array2D.h"
#include "Array3D.h"

enum BoundaryPatch{FIXED, ZERO_GRADIENT};
enum{CONSERVED, AUXILLARY, PRIMITIVE};

template <class T>
class Field : public Array3D<T>
{
private:

    //- Face nodes for primitive and conserved fields

    Array3D<T> facesI_;
    Array3D<T> facesJ_;
    Array3D<T> facesK_;

    //- Face fluxes conserved fields (units/m^2)

    Array3D<T> faceFluxesI_;
    Array3D<T> faceFluxesJ_;
    Array3D<T> faceFluxesK_;

    //- Boundary patches

    BoundaryPatch eastBoundaryPatch_;
    BoundaryPatch westBoundaryPatch_;
    BoundaryPatch northBoundaryPatch_;
    BoundaryPatch southBoundaryPatch_;
    BoundaryPatch topBoundaryPatch_;
    BoundaryPatch bottomBoundaryPatch_;

public:

    Field(std::string name = "UnnamedField", int type = AUXILLARY);
    Field(int nI, int nJ, int nK, std::string name = "UnnamedField", int type = AUXILLARY);
    Field(const Field& other);

    /**
     * @brief A CONSERVED type allocates memory for both face values and face flux values. A PRIMITIVE type allocates
     * memory only for face values. An AUXILLARY type only allocates the centered values.
     */
    int type;
    std::string name;

    void allocate(int nI, int nJ, int nK);

    //- The "type" determines whether or not tranport equations need to be solved for this field

    //- Access

    T& operator()(int i, int j, int k);
    T& operator()(int k);

    T& faceE(int i, int j, int k){ return facesI_(i + 1, j, k); }
    T& faceW(int i, int j, int k){ return facesI_(i, j, k); }
    T& faceN(int i, int j, int k){ return facesJ_(i, j + 1, k); }
    T& faceS(int i, int j, int k){ return facesJ_(i, j, k); }
    T& faceT(int i, int j, int k){ return facesK_(i, j, k + 1); }
    T& faceB(int i, int j, int k){ return facesK_(i, j, k); }

    T& faceI(int i, int j, int k){ return facesI_(i, j, k); }
    T& faceJ(int i, int j, int k){ return facesJ_(i, j, k); }
    T& faceK(int i, int j, int k){ return facesK_(i, j, k); }

    T& fluxI(int i, int j, int k){ return faceFluxesI_(i, j, k); }
    T& fluxJ(int i, int j, int k){ return faceFluxesJ_(i, j, k); }
    T& fluxK(int i, int j, int k){ return faceFluxesK_(i, j, k); }

    T fluxE(int i, int j, int k){ return faceFluxesI_(i + 1, j, k); }
    T fluxW(int i, int j, int k){ return -faceFluxesI_(i, j, k); }
    T fluxN(int i, int j, int k){ return faceFluxesJ_(i, j + 1, k); }
    T fluxS(int i, int j, int k){ return -faceFluxesJ_(i, j, k); }
    T fluxT(int i, int j, int k){ return faceFluxesK_(i, j, k + 1); }
    T fluxB(int i, int j, int k){ return -faceFluxesK_(i, j, k); }
    T sumFluxes(int i, int j, int k);

    /**
     * @brief Get a stencil from the specified cell
     * @return A 3D array containing the stencil.
     */
    Array3D<T> getStencil(int i, int j, int k);

    /**
     * @brief Place a stencil in the array provided from the specified cell
     * @param stencil A 3D array to contain the stencil.
     */
    void getStencil(int i, int j, int k, Array3D<T>& stencil);

    /**
     * @brief Assign boundaries a patch that identifies the type, as well as a reference value.
     */
    void setAllBoundaries(BoundaryPatch eastBoundaryType, T eastBoundaryValue,
                          BoundaryPatch westBoundaryType, T westBoundaryValue,
                          BoundaryPatch northBoundaryType, T northBoundaryValue,
                          BoundaryPatch southBoundaryType, T southBoundaryValue,
                          BoundaryPatch topBoundaryType, T topBoundaryValue,
                          BoundaryPatch bottomBoundaryType, T bottomBoundaryValue);

    void setEastBoundary(BoundaryPatch boundaryType, T boundaryValue);
    void setWestBoundary(BoundaryPatch boundaryType, T boundaryValue);
    void setNorthBoundary(BoundaryPatch BoundaryType, T boundaryValue);
    void setSouthBoundary(BoundaryPatch BoundaryType, T boundaryValue);
    void setTopBoundary(BoundaryPatch BoundaryType, T boundaryValue);
    void setBottomBoundary(BoundaryPatch BoundaryType, T boundaryValue);

    void setBoundaryFields();

    void setEastBoundaryField();
    void setWestBoundaryField();
    void setNorthBoundaryField();
    void setSouthBoundaryField();
    void setTopBoundaryField();
    void setBottomBoundaryField();

    /**
     * @brief Print the field to the console, for debugging purposes.
     */
    void print();
};

#include "FieldI.h"

#endif
