/**
 * @file    FieldI.h
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
 * This file contains the templated method implementations for class
 * Field.
 */

#include "Field.h"
#include "Output.h"

// ************* Constructors and Destructors *************

template <class T>
Field<T>::Field(std::string name, int type)
    :
      name(name),
      type(type)
{

}

template<class T>
Field<T>::Field(int nI, int nJ, int nK, std::string name, int type)
    :
      name(name),
      type(type)
{
    allocate(nI, nJ, nK);
}

template<class T>
Field<T>::Field(const Field &other)
    :
      Field(other.nI_, other.nJ_, other.nK_, other.name, other.type)
{
    int i, j, k;

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

// ************* Public Methods *************

template<class T>
void Field<T>::allocate(int nI, int nJ, int nK)
{
    Array3D<T>::allocate(nI, nJ, nK);

    if (type == CONSERVED || type == PRIMITIVE)
    {
        facesI_.allocate(nI + 1, nJ, nK);
        facesJ_.allocate(nI, nJ + 1, nK);
        facesK_.allocate(nI, nJ, nK + 1);
    }

    if (type == CONSERVED)
    {
        faceFluxesI_.allocate(nI + 1, nJ, nK);
        faceFluxesJ_.allocate(nI, nJ + 1, nK);
        faceFluxesK_.allocate(nI, nJ, nK + 1);
    }

    eastBoundaryField_.allocate(Array3D<T>::nJ_, Array3D<T>::nK_);
    westBoundaryField_.allocate(Array3D<T>::nJ_, Array3D<T>::nK_);
    northBoundaryField_.allocate(Array3D<T>::nI_, Array3D<T>::nK_);
    southBoundaryField_.allocate(Array3D<T>::nI_, Array3D<T>::nK_);
    topBoundaryField_.allocate(Array3D<T>::nI_, Array3D<T>::nJ_);
    bottomBoundaryField_.allocate(Array3D<T>::nI_, Array3D<T>::nJ_);
}

template<class T>
Array3D<T> Field<T>::getStencil(int i, int j, int k)
{
    Array3D<T> stencil(3, 3, 3);
    int l, m, n;

    for(n = -1; n < 2; ++n)
    {
        for(m = -1; m < 2; ++m)
        {
            for(l = -1; l < 2; ++l)
            {
                stencil(l, m, n) = operator ()(i + l, j + m, k + n);
            }
        }
    }
}

template<class T>
void Field<T>::getStencil(int i, int j, int k, Array3D<T> &stencil)
{
    int l, m, n;

    for(n = -1; n < 2; ++n)
    {
        for(m = -1; m < 2; ++m)
        {
            for(l = -1; l < 2; ++l)
            {
                stencil(l, m, n) = operator ()(i + l, j + m, k + n);
            }
        }
    }
}

template<class T>
T& Field<T>::operator()(int i, int j, int k)
{
    if(i >= 0 && j >= 0 && k >= 0 &&
            i < Array3D<T>::nI_ && j < Array3D<T>::nJ_ && k < Array3D<T>::nK_)
    {
        return Array3D<T>::data_[i][j][k];
    }

    // Access to the boundary fields

    if(i < 0)
    {
        return westBoundaryField_(j, k);
    }
    else if (i >= Array3D<T>::nI_)
    {
        return eastBoundaryField_(j, k);
    }

    if(j < 0)
    {
        return southBoundaryField_(i, k);
    }
    else if (j >= Array3D<T>::nJ_)
    {
        return northBoundaryField_(i, k);
    }

    if(k < 0)
    {
        return bottomBoundaryField_(i, j);
    }
    else if (k >= Array3D<T>::nK_)
    {
        return topBoundaryField_(i, j);
    }

    // Just to get rid of the compiler warning

    return Array3D<T>::data_[i][j][k];
}

template<class T>
T Field<T>::sumFluxes(int i, int j, int k)
{
    return faceFluxesI_(i + 1, j, k) - faceFluxesI_(i, j, k) + faceFluxesJ_(i, j + 1, k) - faceFluxesJ_(i, j, k) + faceFluxesK_(i, j, k + 1) - faceFluxesK_(i, j, k);
}

template<class T>
void Field<T>::print()
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

template<class T>
void Field<T>::setAllBoundaries(BoundaryPatch eastBoundaryType, T eastBoundaryValue,
                                BoundaryPatch westBoundaryType, T westBoundaryValue,
                                BoundaryPatch northBoundaryType, T northBoundaryValue,
                                BoundaryPatch southBoundaryType, T southBoundaryValue,
                                BoundaryPatch topBoundaryType, T topBoundaryValue,
                                BoundaryPatch bottomBoundaryType, T bottomBoundaryValue)
{
    setEastBoundary(eastBoundaryType, eastBoundaryValue);
    setWestBoundary(westBoundaryType, westBoundaryValue);
    setNorthBoundary(northBoundaryType, northBoundaryValue);
    setSouthBoundary(southBoundaryType, southBoundaryValue);
    setTopBoundary(topBoundaryType, topBoundaryValue);
    setBottomBoundary(bottomBoundaryType, bottomBoundaryValue);

    Output::print("Field", "boundaries for field \"" + name + "\" have been set.");
}

template<class T>
void Field<T>::setEastBoundary(BoundaryPatch boundaryType, T boundaryValue)
{
    int j, k;

    eastBoundaryPatch_ = boundaryType;

    for(k = 0; k < Array3D<T>::nK_; ++k)
    {
        for(j = 0; j < Array3D<T>::nJ_; ++j)
        {
            eastBoundaryField_(j, k) = boundaryValue;
        } // end for j
    } //  end for k
}

template<class T>
void Field<T>::setWestBoundary(BoundaryPatch boundaryType, T boundaryValue)
{
    int j, k;

    westBoundaryPatch_ = boundaryType;

    for(k = 0; k < Array3D<T>::nK_; ++k)
    {
        for(j = 0; j < Array3D<T>::nJ_; ++j)
        {
            westBoundaryField_(j, k) = boundaryValue;
        } // end for j
    } //  end for k
}

template<class T>
void Field<T>::setNorthBoundary(BoundaryPatch boundaryType, T boundaryValue)
{
    int i, k;

    northBoundaryPatch_ = boundaryType;

    for(k = 0; k < Array3D<T>::nK_; ++k)
    {
        for(i = 0; i < Array3D<T>::nI_; ++i)
        {
            northBoundaryField_(i, k) = boundaryValue;
        } // end for i
    } //  end for k
}

template<class T>
void Field<T>::setSouthBoundary(BoundaryPatch boundaryType, T boundaryValue)
{
    int i, k;

    southBoundaryPatch_ = boundaryType;

    for(k = 0; k < Array3D<T>::nK_; ++k)
    {
        for(i = 0; i < Array3D<T>::nI_; ++i)
        {
            southBoundaryField_(i, k) = boundaryValue;
        } // end for i
    } //  end for k
}

template<class T>
void Field<T>::setTopBoundary(BoundaryPatch boundaryType, T boundaryValue)
{
    int i, j;

    topBoundaryPatch_ = boundaryType;

    for(j = 0; j < Array3D<T>::nJ_; ++j)
    {
        for(i = 0; i < Array3D<T>::nI_; ++i)
        {
            topBoundaryField_(i, j) = boundaryValue;
        } // end for i
    } //  end for j
}

template<class T>
void Field<T>::setBottomBoundary(BoundaryPatch boundaryType, T boundaryValue)
{
    int i, j;

    bottomBoundaryPatch_ = boundaryType;

    for(j = 0; j < Array3D<T>::nJ_; ++j)
    {
        for(i = 0; i < Array3D<T>::nI_; ++i)
        {
            bottomBoundaryField_(i, j) = boundaryValue;
        } // end for i
    } //  end for j
}

template<class T>
void Field<T>::setBoundaryFields()
{
    setEastBoundaryField();
    setWestBoundaryField();
    setNorthBoundaryField();
    setSouthBoundaryField();
    setTopBoundaryField();
    setBottomBoundaryField();
}

template<class T>
void Field<T>::setEastBoundaryField()
{
    int j, k, nJ(eastBoundaryField_.sizeI()), nK(eastBoundaryField_.sizeJ()), maxI(Array3D<T>::nI_ - 1);

    for(k = 0; k < nK; ++k)
    {
        for(j = 0; j < nJ; ++j)
        {
            switch(eastBoundaryPatch_)
            {
            case FIXED:

                break;

            case ZERO_GRADIENT:

                eastBoundaryField_(j, k) = Array3D<T>::data_[maxI][j][k];

                break;
            };
        } // end for j
    } // end for k
}

template<class T>
void Field<T>::setWestBoundaryField()
{
    int j, k, nJ(westBoundaryField_.sizeI()), nK(westBoundaryField_.sizeJ());

    for(k = 0; k < nK; ++k)
    {
        for(j = 0; j < nJ; ++j)
        {
            switch(westBoundaryPatch_)
            {

            case FIXED:

                break;

            case ZERO_GRADIENT:

                westBoundaryField_(j, k) = Array3D<T>::data_[0][j][k];

                break;
            };
        } // end for j
    } // end for k
}

template<class T>
void Field<T>::setNorthBoundaryField()
{
    int i, k, nI(northBoundaryField_.sizeI()), nK(northBoundaryField_.sizeJ()), maxJ(Array3D<T>::nJ_ - 1);

    for(k = 0; k < nK; ++k)
    {
        for(i = 0; i < nI; ++i)
        {
            switch(northBoundaryPatch_)
            {

            case FIXED:

                break;

            case ZERO_GRADIENT:

                northBoundaryField_(i, k) = Array3D<T>::data_[i][maxJ][k];

                break;
            };
        } // end for i
    } // end for k

}

template<class T>
void Field<T>::setSouthBoundaryField()
{
    int i, k, nI(southBoundaryField_.sizeI()), nK(southBoundaryField_.sizeJ());

    for(k = 0; k < nK; ++k)
    {

        for(i = 0; i < nI; ++i)
        {

            switch(southBoundaryPatch_)
            {
            case FIXED:

                break;

            case ZERO_GRADIENT:

                southBoundaryField_(i, k) = Array3D<T>::data_[i][0][k];

                break;
            };
        } // end for i
    } // end for k
}

template<class T>
void Field<T>::setTopBoundaryField()
{
    int i, j, nI(topBoundaryField_.sizeI()), nJ(topBoundaryField_.sizeJ()), maxK(Array3D<T>::nK_ - 1);

    for(j = 0; j < nJ; ++j)
    {
        for(i = 0; i < nI; ++i)
        {
            switch(topBoundaryPatch_)
            {
            case FIXED:

                break;

            case ZERO_GRADIENT:

                topBoundaryField_(i, j) = Array3D<T>::data_[i][j][maxK];

                break;
            };
        } // end for i
    } // end for j

}

template<class T>
void Field<T>::setBottomBoundaryField()
{
    int i, j, nI(bottomBoundaryField_.sizeI()), nJ(bottomBoundaryField_.sizeJ());

    for(j = 0; j < nJ; ++j)
    {
        for(i = 0; i < nI; ++i)
        {
            switch(bottomBoundaryPatch_)
            {
            case FIXED:

                break;

            case ZERO_GRADIENT:

                bottomBoundaryField_(i, j) = Array3D<T>::data_[i][j][0];

                break;
            };
        } // end for i
    } // end for j
}
