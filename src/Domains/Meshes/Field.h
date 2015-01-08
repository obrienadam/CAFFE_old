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

    //- Input

    void setAllBoundaries(BoundaryPatch eastBoundaryType,
                          BoundaryPatch westBoundaryType,
                          BoundaryPatch northBoundaryType,
                          BoundaryPatch southBoundaryType,
                          BoundaryPatch topBoundaryType,
                          BoundaryPatch bottomBoundaryType);

    void setEastBoundary(BoundaryPatch boundaryType);
    void setWestBoundary(BoundaryPatch boundaryType);
    void setNorthBoundary(BoundaryPatch BoundaryType);
    void setSouthBoundary(BoundaryPatch BoundaryType);
    void setTopBoundary(BoundaryPatch BoundaryType);
    void setBottomBoundary(BoundaryPatch BoundaryType);

    //- Debug

    void print();

};

#include "FieldI.h"

#endif
