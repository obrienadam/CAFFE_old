#ifndef FIELD_H
#define FIELD_H

#include <string>

#include "Array3D.h"

enum{CONSERVED, AUXILLARY};

template <class T>
class Field : public Array3D<T>
{

private:

    //- Face fluxes conserved fields (units/m^2)

    Array3D<T> faceFluxesI_;
    Array3D<T> faceFluxesJ_;
    Array3D<T> faceFluxesK_;

public:

    Field(std::string name = "UnnamedField", int type = AUXILLARY)
        :
          name(name),
          type(type)
    {

    }

    Field(int nI, int nJ, int nK, std::string name = "UnnamedField", int type = AUXILLARY)
        :
          Array3D<T>(nI, nJ, nK),
          name(name),
          type(type)
    {

        if (type == CONSERVED)
        {

            faceFluxesI_.allocate(nI + 1, nJ, nK);
            faceFluxesJ_.allocate(nI, nJ + 1, nK);
            faceFluxesK_.allocate(nI, nJ, nK + 1);

        }

    }

    void resize(int nI, int nJ, int nK)
    {

        Array3D<T>::allocate(nI, nJ, nK);

        if (type == CONSERVED)
        {

            faceFluxesI_.allocate(nI + 1, nJ, nK);
            faceFluxesJ_.allocate(nI, nJ + 1, nK);
            faceFluxesK_.allocate(nI, nJ, nK + 1);

        }

    }

    //- The "type" determines whether or not tranport equations need to be solved for this field

    int type;
    std::string name;

    //- Access

    T& fluxE(int i, int j, int k){ return faceFluxesI_(i + 1, j, k); }
    T& fluxW(int i, int j, int k){ return faceFluxesI_(i, j, k); }
    T& fluxN(int i, int j, int k){ return faceFluxesJ_(i, j + 1, k); }
    T& fluxS(int i, int j, int k){ return faceFluxesJ_(i, j, k); }
    T& fluxT(int i, int j, int k){ return faceFluxesK_(i, j, k + 1); }
    T& fluxB(int i, int j, int k){ return faceFluxesK_(i, j, k); }

};

#endif
