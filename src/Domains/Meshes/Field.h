#ifndef FIELD_H
#define FIELD_H

#include <string>
#include <iostream>

#include "Array3D.h"

enum{CONSERVED, AUXILLARY, PRIMITIVE};

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

            faceFluxesI_.allocate(Array3D<T>::nI_ + 1, Array3D<T>::nJ_, Array3D<T>::nK_);
            faceFluxesJ_.allocate(Array3D<T>::nI_, Array3D<T>::nJ_ + 1, Array3D<T>::nK_);
            faceFluxesK_.allocate(Array3D<T>::nI_, Array3D<T>::nJ_, Array3D<T>::nK_ + 1);

        }

    }

    Field(const Field& other)
        :
          Array3D<T>(other.nI_, other.nJ_, other.nK_),
          name(other.name),
          type(other.type)
    {

        int i, j, k;

        if (type == CONSERVED)
        {

            faceFluxesI_.allocate(Array3D<T>::nI_ + 1, Array3D<T>::nJ_, Array3D<T>::nK_);
            faceFluxesJ_.allocate(Array3D<T>::nI_, Array3D<T>::nJ_ + 1, Array3D<T>::nK_);
            faceFluxesK_.allocate(Array3D<T>::nI_, Array3D<T>::nJ_, Array3D<T>::nK_ + 1);

        }

        for(k = 0; k < Array3D<T>::nK_; ++k)
        {

            for(j = 0; j < Array3D<T>::nJ_; ++j)

            {

                for(i = 0; i < Array3D<T>::nI_; ++i)

                {

                    Array3D<T>::data_[i][j][k] = other.Array3D<T>::data_[i][j][k];

                } // end for i
            } // end for j
        } // end for k

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

    //- Debug

    void print()
    {
        using namespace std;

        int i, j, k;

        cout << "Field name = " << "\"" << name << "\"\n";

        for(k = 0; k < Array3D<T>::nK_; ++k)
        {

            for(j = 0; j < Array3D<T>::nJ_; ++j)
            {

                for(i = 0; i < Array3D<T>::nI_; ++i)
                {

                    cout << Array3D<T>::data_[i][j][k] << " ";

                } // end for i

                cout << endl;

            } // end for j

            cout << endl;

        } // end for k

    }

};

#endif
