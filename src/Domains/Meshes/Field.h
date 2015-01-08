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

    //- Face fluxes conserved fields (units/m^2)

    Array3D<T> faceFluxesI_;
    Array3D<T> faceFluxesJ_;
    Array3D<T> faceFluxesK_;

    //- Boundary patches

    Array2D<BoundaryPatch> eastBoundaryPatch_;
    Array2D<BoundaryPatch> westBoundaryPatch_;
    Array2D<BoundaryPatch> northBoundaryPatch_;
    Array2D<BoundaryPatch> southBoundaryPatch_;
    Array2D<BoundaryPatch> topBoundaryPatch_;
    Array2D<BoundaryPatch> bottomBoundaryPatch_;

    //- Boundary fields

    Array2D<T> eastBoundaryField_;
    Array2D<T> westBoundaryField_;
    Array2D<T> northBoundaryField_;
    Array2D<T> southBoundaryField_;
    Array2D<T> topBoundaryField_;
    Array2D<T> bottomBoundaryField_;

public:

    Field(std::string name = "UnnamedField", int type = AUXILLARY);
    Field(int nI, int nJ, int nK, std::string name = "UnnamedField", int type = AUXILLARY);
    Field(const Field& other);

    int type;
    std::string name;

    void resize(int nI, int nJ, int nK);

    //- The "type" determines whether or not tranport equations need to be solved for this field

    //- Access

    T& fluxE(int i, int j, int k){ return faceFluxesI_(i + 1, j, k); }
    T& fluxW(int i, int j, int k){ return faceFluxesI_(i, j, k); }
    T& fluxN(int i, int j, int k){ return faceFluxesJ_(i, j + 1, k); }
    T& fluxS(int i, int j, int k){ return faceFluxesJ_(i, j, k); }
    T& fluxT(int i, int j, int k){ return faceFluxesK_(i, j, k + 1); }
    T& fluxB(int i, int j, int k){ return faceFluxesK_(i, j, k); }

    //- Boundary related methods

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

    //- Debug

    void print();

};

#include "FieldI.h"

#endif
