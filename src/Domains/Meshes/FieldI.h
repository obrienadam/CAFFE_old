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
    int k;

    for(k = 0; k < Array3D<T>::n_; ++k)
    {
        Array3D<T>::data_[k] = other.Array3D<T>::data_[k];
    }
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
        return Array3D<T>::operator ()(i, j, k);
    }

    // Access to the boundaries

    if(i < 0)
    {
        return facesI_(0, j, k);
    }
    else if (i >= Array3D<T>::nI_)
    {
        return facesI_(Array3D<T>::nI_, j, k);
    }

    if(j < 0)
    {
        return facesJ_(i, 0, k);
    }
    else if (j >= Array3D<T>::nJ_)
    {
        return facesJ_(i, Array3D<T>::nJ_, k);
    }

    if(k < 0)
    {
        return facesK_(i, j, 0);
    }
    else if (k >= Array3D<T>::nK_)
    {
        return facesK_(i, j, Array3D<T>::nK_);
    }

    // Just to get rid of the compiler warning

    return Array3D<T>::data_[0];
}

template<class T>
T& Field<T>::operator ()(int i, int j, int k, int faceNo)
{
    switch(faceNo)
    {
    case 0:
        return operator ()(i + 1, j, k);
    case 1:
        return operator ()(i - 1, j, k);
    case 2:
        return operator ()(i, j + 1, k);
    case 3:
        return operator ()(i, j - 1, k);
    case 4:
        return operator ()(i, j, k + 1);
    case 5:
        return operator ()(i, j, k - 1);
    default:
        Output::raiseException("Field", "operator()", "invalid face number specified.");
    };

    return operator ()(i, j, k);
}

template<class T>
T& Field<T>::operator ()(int k)
{
    if(k < 0 || k >= Array3D<T>::n_)
        Output::raiseException("Field", "operator()", "Attempted to access an element outside the bounds of the field.");

    return Array3D<T>::data_[k];
}

template<class T>
T Field<T>::sumFluxes(int i, int j, int k)
{
    return faceFluxesI_(i + 1, j, k) - faceFluxesI_(i, j, k)
            + faceFluxesJ_(i, j + 1, k) - faceFluxesJ_(i, j, k)
            + faceFluxesK_(i, j, k + 1) - faceFluxesK_(i, j, k);
}

template<class T>
T Field<T>::maxNeighbour(int i, int j, int k)
{
    using namespace std;

    T maxNb = operator ()(i, j, k);

    maxNb = max(maxNb, operator()(i + 1, j, k));
    maxNb = max(maxNb, operator()(i - 1, j, k));
    maxNb = max(maxNb, operator()(i, j + 1, k));
    maxNb = max(maxNb, operator()(i, j - 1, k));
    maxNb = max(maxNb, operator()(i, j, k + 1));
    maxNb = max(maxNb, operator()(i, j, k - 1));

    return maxNb;
}

template<class T>
T Field<T>::minNeighbour(int i, int j, int k)
{
    using namespace std;

    T minNb = operator ()(i, j, k);

    minNb = min(minNb, operator()(i + 1, j, k));
    minNb = min(minNb, operator()(i - 1, j, k));
    minNb = min(minNb, operator()(i, j + 1, k));
    minNb = min(minNb, operator()(i, j - 1, k));
    minNb = min(minNb, operator()(i, j, k + 1));
    minNb = min(minNb, operator()(i, j, k - 1));

    return minNb;
}

template<class T>
void Field<T>::print()
{
    using namespace std;

    int i, j, k, l;

    cout << "Field name = " << "\"" << name << "\"\n";

    for(k = 0, l = 0; k < Array3D<T>::nK_; ++k)
    {
        for(j = 0; j < Array3D<T>::nJ_; ++j)
        {
            for(i = 0; i < Array3D<T>::nI_; ++i, ++l)
            {
                cout << Array3D<T>::data_[l] << " ";
            } // end for i

            cout << endl;
        } // end for j

        cout << endl;
    } // end for k
}

template<class T>
T& Field<T>::face(int i, int j, int k, int faceNo)
{
    switch(faceNo)
    {
    case 0:
        return facesI_(i + 1, j, k);
    case 1:
        return facesI_(i, j, k);
    case 2:
        return facesJ_(i, j + 1, k);
    case 3:
        return facesJ_(i, j, k);
    case 4:
        return facesK_(i, j, k + 1);
    case 5:
        return facesK_(i, j, k);
    default:
        Output::raiseException("Field", "face", "attempted to access a face that does not exist.");
    };
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
            facesI_(Array3D<T>::nI_, j, k) = boundaryValue;
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
            facesI_(0, j, k) = boundaryValue;
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
            facesJ_(i, Array3D<T>::nJ_, k) = boundaryValue;
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
            facesJ_(i, 0, k) = boundaryValue;
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
            facesK_(i, j, Array3D<T>::nK_) = boundaryValue;
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
            facesK_(i, j, 0) = boundaryValue;
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
    int j, k, nJ, nK, maxI;

    switch(eastBoundaryPatch_)
    {

    case FIXED: case EXTRAPOLATE:

        break;

    case ZERO_GRADIENT:

        nJ = Array3D<T>::nJ_;
        nK = Array3D<T>::nK_;
        maxI = Array3D<T>::nI_ - 1;

        for(k = 0; k < nK; ++k)
        {
            for(j = 0; j < nJ; ++j)
            {
                facesI_(Array3D<T>::nI_, j, k) = Array3D<T>::operator ()(maxI, j, k);
            }
        }

        break;
    };
}

template<class T>
void Field<T>::setWestBoundaryField()
{
    int j, k, nJ, nK;

    switch(westBoundaryPatch_)
    {

    case FIXED: case EXTRAPOLATE:

        break;

    case ZERO_GRADIENT:

        nJ = Array3D<T>::nJ_;
        nK = Array3D<T>::nK_;

        for(k = 0; k < nK; ++k)
        {
            for(j = 0; j < nJ; ++j)
            {
                facesI_(0, j, k) = Array3D<T>::operator ()(0, j, k);
            }
        }

        break;
    };
}

template<class T>
void Field<T>::setNorthBoundaryField()
{
    int i, k, nI, nK, maxJ;

    switch(northBoundaryPatch_)
    {

    case FIXED: case EXTRAPOLATE:

        break;

    case ZERO_GRADIENT:

        nI = Array3D<T>::nI_;
        nK = Array3D<T>::nK_;
        maxJ = Array3D<T>::nJ_ - 1;

        for(k = 0; k < nK; ++k)
        {
            for(i = 0; i < nI; ++i)
            {
                facesJ_(i, Array3D<T>::nJ_, k) = Array3D<T>::operator ()(i, maxJ, k);
            }
        }

        break;
    };
}

template<class T>
void Field<T>::setSouthBoundaryField()
{
    int i, k, nI, nK;

    switch(southBoundaryPatch_)
    {
    case FIXED: case EXTRAPOLATE:

        break;

    case ZERO_GRADIENT:

        nI = Array3D<T>::nI_;
        nK = Array3D<T>::nK_;

        for(k = 0; k < nK; ++k)
        {
            for(i = 0; i < nI; ++i)
            {
                facesJ_(i, 0, k) = Array3D<T>::operator ()(i, 0, k);
            }
        }

        break;
    };
}

template<class T>
void Field<T>::setTopBoundaryField()
{
    int i, j, nI, nJ, maxK;

    switch(topBoundaryPatch_)
    {
    case FIXED: case EXTRAPOLATE:

        break;

    case ZERO_GRADIENT:

        nI = Array3D<T>::nI_;
        nJ = Array3D<T>::nJ_;
        maxK = Array3D<T>::nK_ - 1;

        for(j = 0; j < nJ; ++j)
        {
            for(i = 0; i < nI; ++i)
            {
                facesK_(i, j, Array3D<T>::nK_) = Array3D<T>::operator ()(i, j, maxK);
            }
        }

        break;
    };
}

template<class T>
void Field<T>::setBottomBoundaryField()
{
    int i, j, nI, nJ;

    switch(bottomBoundaryPatch_)
    {
    case FIXED: case EXTRAPOLATE:

        break;

    case ZERO_GRADIENT:

        nI = Array3D<T>::nI_;
        nJ = Array3D<T>::nJ_;

        for(j = 0; j < nJ; ++j)
        {
            for(i = 0; i < nI; ++i)
            {
                facesK_(i, j, 0) = Array3D<T>::operator ()(i, j, 0);
            }
        }

        break;
    };
}
