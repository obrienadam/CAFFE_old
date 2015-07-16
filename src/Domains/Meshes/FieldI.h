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

template<class T>
Field<T>::Field(const HexaFvmMesh &mesh, FieldType type, std::string name)
    :
      mesh_(mesh),
      type_(type),
      name(name)
{
    allocate(mesh_.nCellsI(), mesh_.nCellsJ(), mesh_.nCellsK());
}

template<class T>
Field<T>::Field(const Field &other)
    :
      mesh_(other.mesh_),
      type_(other.type_),
      name(other.name)
{
    allocate(mesh_.nCellsI(), mesh_.nCellsJ(), mesh_.nCellsK());
}

// ************* Public Methods *************

template<class T>
Field<T>& Field<T>::operator =(const Field<T>& other)
{
    int i;

    //- Check to make sure this assignment is legal. Fields cannont be equated to other fields that reference a different mesh
    if(&mesh_ != &other.mesh_)
        Output::raiseException("Field", "operator=", "cannot assign a field to another field if the mesh references are different.");

    for(i = 0; i < Array3D<T>::n_; ++i)
        Array3D<T>::data_[i] = other.data_[i];

    if(type_ == CONSERVED && other.type_ == CONSERVED)
    {
        for(i = 0; i < facesI_.size(); ++i)
            facesI_(i) = other.facesI_(i);

        for(i = 0; i < facesJ_.size(); ++i)
            facesJ_(i) = other.facesJ_(i);

        for(i = 0; i < facesK_.size(); ++i)
            facesK_(i) = other.facesK_(i);
    }

    return *this;
}

template<class T>
void Field<T>::allocate(int nCellsI, int nCellsJ, int nCellsK)
{
    Array3D<T>::allocate(nCellsI, nCellsJ, nCellsK);

    if (type_ == CONSERVED)
    {
        facesI_.allocate(nCellsI + 1, nCellsJ, nCellsK);
        facesJ_.allocate(nCellsI, nCellsJ + 1, nCellsK);
        facesK_.allocate(nCellsI, nCellsJ, nCellsK + 1);
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
const T& Field<T>::operator()(int i, int j, int k) const
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
const T& Field<T>::operator ()(int i, int j, int k, int faceNo) const
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
const T& Field<T>::operator ()(int k) const
{
    if(k < 0 || k >= Array3D<T>::n_)
        Output::raiseException("Field", "operator()", "Attempted to access an element outside the bounds of the field.");

    return Array3D<T>::data_[k];
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
void Field<T>::setEastBoundary(const std::string &boundaryType, T boundaryValue)
{
    int j, k;

    if(boundaryType == "fixed")
        eastBoundaryPatch_ = FIXED;
    else if(boundaryType == "zeroGradient")
        eastBoundaryPatch_ = ZERO_GRADIENT;
    else if(boundaryType == "empty")
        eastBoundaryPatch_ = EMPTY;
    else
        Output::raiseException("Field", "setEastBoundary", "unrecognized boundary condition type \"" + boundaryType + "\".");

    for(k = 0; k < Array3D<T>::nK_; ++k)
    {
        for(j = 0; j < Array3D<T>::nJ_; ++j)
        {
            facesI_(Array3D<T>::nI_, j, k) = boundaryValue;
        } // end for j
    } //  end for k
}

template<class T>
void Field<T>::setWestBoundary(const std::string &boundaryType, T boundaryValue)
{
    int j, k;

    if(boundaryType == "fixed")
        westBoundaryPatch_ = FIXED;
    else if(boundaryType == "zeroGradient")
        westBoundaryPatch_ = ZERO_GRADIENT;
    else if(boundaryType == "empty")
        westBoundaryPatch_ = EMPTY;
    else
        Output::raiseException("Field", "setWestBoundary", "unrecognized boundary condition type \"" + boundaryType + "\".");

    for(k = 0; k < Array3D<T>::nK_; ++k)
    {
        for(j = 0; j < Array3D<T>::nJ_; ++j)
        {
            facesI_(0, j, k) = boundaryValue;
        } // end for j
    } //  end for k
}

template<class T>
void Field<T>::setNorthBoundary(const std::string &boundaryType, T boundaryValue)
{
    int i, k;

    if(boundaryType == "fixed")
        northBoundaryPatch_ = FIXED;
    else if(boundaryType == "zeroGradient")
        northBoundaryPatch_ = ZERO_GRADIENT;
    else if(boundaryType == "empty")
        northBoundaryPatch_ = EMPTY;
    else
        Output::raiseException("Field", "setNorthBoundary", "unrecognized boundary condition type \"" + boundaryType + "\".");

    for(k = 0; k < Array3D<T>::nK_; ++k)
    {
        for(i = 0; i < Array3D<T>::nI_; ++i)
        {
            facesJ_(i, Array3D<T>::nJ_, k) = boundaryValue;
        } // end for i
    } //  end for k
}

template<class T>
void Field<T>::setSouthBoundary(const std::string &boundaryType, T boundaryValue)
{
    int i, k;

    if(boundaryType == "fixed")
        southBoundaryPatch_ = FIXED;
    else if(boundaryType == "zeroGradient")
        southBoundaryPatch_ = ZERO_GRADIENT;
    else if(boundaryType == "empty")
        southBoundaryPatch_ = EMPTY;
    else
        Output::raiseException("Field", "setSouthBoundary", "unrecognized boundary condition type \"" + boundaryType + "\".");

    for(k = 0; k < Array3D<T>::nK_; ++k)
    {
        for(i = 0; i < Array3D<T>::nI_; ++i)
        {
            facesJ_(i, 0, k) = boundaryValue;
        } // end for i
    } //  end for k
}

template<class T>
void Field<T>::setTopBoundary(const std::string &boundaryType, T boundaryValue)
{
    int i, j;

    if(boundaryType == "fixed")
        topBoundaryPatch_ = FIXED;
    else if(boundaryType == "zeroGradient")
        topBoundaryPatch_ = ZERO_GRADIENT;
    else if(boundaryType == "empty")
        topBoundaryPatch_ = EMPTY;
    else
        Output::raiseException("Field", "setTopBoundary", "unrecognized boundary condition type \"" + boundaryType + "\".");

    for(j = 0; j < Array3D<T>::nJ_; ++j)
    {
        for(i = 0; i < Array3D<T>::nI_; ++i)
        {
            facesK_(i, j, Array3D<T>::nK_) = boundaryValue;
        } // end for i
    } //  end for j
}

template<class T>
void Field<T>::setBottomBoundary(const std::string &boundaryType, T boundaryValue)
{
    int i, j;

    if(boundaryType == "fixed")
        bottomBoundaryPatch_ = FIXED;
    else if(boundaryType == "zeroGradient")
        bottomBoundaryPatch_ = ZERO_GRADIENT;
    else if(boundaryType == "empty")
        bottomBoundaryPatch_ = EMPTY;
    else
        Output::raiseException("Field", "setBottomBoundary", "unrecognized boundary condition type \"" + boundaryType + "\".");

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

    case FIXED:

        break;

    case ZERO_GRADIENT: case EMPTY:

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

    case FIXED:

        break;

    case ZERO_GRADIENT: case EMPTY:

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

    case FIXED:

        break;

    case ZERO_GRADIENT: case EMPTY:

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
    case FIXED:

        break;

    case ZERO_GRADIENT: case EMPTY:

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
    case FIXED:

        break;

    case ZERO_GRADIENT: case EMPTY:

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
    case FIXED:

        break;

    case ZERO_GRADIENT: case EMPTY:

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

template<class T>
void Field<T>::setImplicitBoundaryCoeffs(int i, int j, int k, double *a, T &b)
{
    if(i == mesh_.uCellI())
    {
        switch(eastBoundaryPatch_)
        {
        case FIXED:
            b += -a[1]*operator ()(i + 1, j, k);
            break;
        case ZERO_GRADIENT:
            a[0] += a[1];
            break;
        case EMPTY:
            a[1] = 0.;
            break;
        }
    }

    if(i == 0)
    {
        switch(westBoundaryPatch_)
        {
        case FIXED:
            b += -a[2]*operator ()(i - 1, j, k);
            break;
        case ZERO_GRADIENT:
            a[0] += a[2];
            break;
        case EMPTY:
            a[2] = 0.;
            break;
        }
    }

    if(j == mesh_.uCellJ())
    {
        switch(northBoundaryPatch_)
        {
        case FIXED:
            b += -a[3]*operator ()(i, j + 1, k);
            break;
        case ZERO_GRADIENT:
            a[0] += a[3];
            break;
        case EMPTY:
            a[3] = 0.;
            break;
        }
    }

    if(j == 0)
    {
        switch(southBoundaryPatch_)
        {
        case FIXED:
            b += -a[4]*operator ()(i, j - 1, k);
            break;
        case ZERO_GRADIENT:
            a[0] += a[4];
            break;
        case EMPTY:
            a[4] = 0.;
            break;
        }
    }

    if(k == mesh_.uCellK())
    {
        switch(topBoundaryPatch_)
        {
        case FIXED:
            b += -a[5]*operator ()(i, j, k + 1);
            break;
        case ZERO_GRADIENT:
            a[0] += a[5];
            break;
        case EMPTY:
            a[5] = 0.;
            break;
        }
    }

    if(k == 0)
    {
        switch(bottomBoundaryPatch_)
        {
        case FIXED:
            b += -a[6]*operator ()(i, j, k - 1);
            break;
        case ZERO_GRADIENT:
            a[0] += a[6];
            break;
        case EMPTY:
            a[6] = 0.;
            break;
        }
    }
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
