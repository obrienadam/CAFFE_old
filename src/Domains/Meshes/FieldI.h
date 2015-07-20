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
      Field<T>::Field(other.mesh_, other.type_, other.name + "_copy")
{
    operator =(other);
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

    eastBoundaryPatch_.allocate(1, mesh_.nCellsJ(), mesh_.nCellsK());
    westBoundaryPatch_.allocate(1, mesh_.nCellsJ(), mesh_.nCellsK());
    northBoundaryPatch_.allocate(mesh_.nCellsI(), 1, mesh_.nCellsK());
    southBoundaryPatch_.allocate(mesh_.nCellsI(), 1, mesh_.nCellsK());
    topBoundaryPatch_.allocate(mesh_.nCellsI(), mesh_.nCellsJ(), 1);
    bottomBoundaryPatch_.allocate(mesh_.nCellsI(), mesh_.nCellsJ(), 1);
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
        return westBoundaryPatch_(abs(i) - 1, j, k);
    }
    else if (i >= Array3D<T>::nI_)
    {
        return eastBoundaryPatch_(i - Array3D<T>::nI_, j, k);
    }

    if(j < 0)
    {
        return southBoundaryPatch_(i, abs(j) - 1, k);
    }
    else if (j >= Array3D<T>::nJ_)
    {
        return northBoundaryPatch_(i, j - Array3D<T>::nJ_, k);
    }

    if(k < 0)
    {
        return bottomBoundaryPatch_(i, j, abs(k) - 1);
    }
    else if (k >= Array3D<T>::nK_)
    {
        return topBoundaryPatch_(i, j, k - Array3D<T>::nK_);
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
        return westBoundaryPatch_(abs(i) - 1, j, k);
    }
    else if (i >= Array3D<T>::nI_)
    {
        return eastBoundaryPatch_(i - Array3D<T>::nI_, j, k);
    }

    if(j < 0)
    {
        return southBoundaryPatch_(i, abs(j) - 1, k);
    }
    else if (j >= Array3D<T>::nJ_)
    {
        return northBoundaryPatch_(i, j - Array3D<T>::nJ_, k);
    }

    if(k < 0)
    {
        return bottomBoundaryPatch_(i, j, abs(k) - 1);
    }
    else if (k >= Array3D<T>::nK_)
    {
        return topBoundaryPatch_(i, j, k - Array3D<T>::nK_);
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
        eastBoundaryPatchId_ = FIXED;
    else if(boundaryType == "zeroGradient")
        eastBoundaryPatchId_ = ZERO_GRADIENT;
    else if(boundaryType == "empty")
        eastBoundaryPatchId_ = EMPTY;
    else
        Output::raiseException("Field", "setEastBoundary", "unrecognized boundary condition type \"" + boundaryType + "\".");

    for(k = 0; k < eastBoundaryPatch_.sizeK(); ++k)
    {
        for(j = 0; j < eastBoundaryPatch_.sizeJ(); ++j)
        {
            eastBoundaryPatch_(0, j, k) = boundaryValue;
            facesI_(mesh_.nCellsI(), j, k) = eastBoundaryPatch_(0, j, k);
        } // end for j
    } //  end for k
}

template<class T>
void Field<T>::setWestBoundary(const std::string &boundaryType, T boundaryValue)
{
    int j, k;

    if(boundaryType == "fixed")
        westBoundaryPatchId_ = FIXED;
    else if(boundaryType == "zeroGradient")
        westBoundaryPatchId_ = ZERO_GRADIENT;
    else if(boundaryType == "empty")
        westBoundaryPatchId_ = EMPTY;
    else
        Output::raiseException("Field", "setWestBoundary", "unrecognized boundary condition type \"" + boundaryType + "\".");

    for(k = 0; k < westBoundaryPatch_.sizeK(); ++k)
    {
        for(j = 0; j < westBoundaryPatch_.sizeJ(); ++j)
        {
            westBoundaryPatch_(0, j, k) = boundaryValue;
            facesI_(0, j, k) = westBoundaryPatch_(0, j, k);
        } // end for j
    } //  end for k
}

template<class T>
void Field<T>::setNorthBoundary(const std::string &boundaryType, T boundaryValue)
{
    int i, k;

    if(boundaryType == "fixed")
        northBoundaryPatchId_ = FIXED;
    else if(boundaryType == "zeroGradient")
        northBoundaryPatchId_ = ZERO_GRADIENT;
    else if(boundaryType == "empty")
        northBoundaryPatchId_ = EMPTY;
    else
        Output::raiseException("Field", "setNorthBoundary", "unrecognized boundary condition type \"" + boundaryType + "\".");

    for(k = 0; k < northBoundaryPatch_.sizeK(); ++k)
    {
        for(i = 0; i < northBoundaryPatch_.sizeI(); ++i)
        {
            northBoundaryPatch_(i, 0, k) = boundaryValue;
            facesJ_(i, mesh_.nCellsJ(), k) = northBoundaryPatch_(i, 0, k);
        } // end for i
    } //  end for k
}

template<class T>
void Field<T>::setSouthBoundary(const std::string &boundaryType, T boundaryValue)
{
    int i, k;

    if(boundaryType == "fixed")
        southBoundaryPatchId_ = FIXED;
    else if(boundaryType == "zeroGradient")
        southBoundaryPatchId_ = ZERO_GRADIENT;
    else if(boundaryType == "empty")
        southBoundaryPatchId_ = EMPTY;
    else
        Output::raiseException("Field", "setSouthBoundary", "unrecognized boundary condition type \"" + boundaryType + "\".");

    for(k = 0; k < southBoundaryPatch_.sizeK(); ++k)
    {
        for(i = 0; i < southBoundaryPatch_.sizeI(); ++i)
        {
            southBoundaryPatch_(i, 0, k) = boundaryValue;
            facesJ_(i, 0, k) = southBoundaryPatch_(i, 0, k);
        } // end for i
    } //  end for k
}

template<class T>
void Field<T>::setTopBoundary(const std::string &boundaryType, T boundaryValue)
{
    int i, j;

    if(boundaryType == "fixed")
        topBoundaryPatchId_ = FIXED;
    else if(boundaryType == "zeroGradient")
        topBoundaryPatchId_ = ZERO_GRADIENT;
    else if(boundaryType == "empty")
        topBoundaryPatchId_ = EMPTY;
    else
        Output::raiseException("Field", "setTopBoundary", "unrecognized boundary condition type \"" + boundaryType + "\".");

    for(j = 0; j < topBoundaryPatch_.sizeJ(); ++j)
    {
        for(i = 0; i < topBoundaryPatch_.sizeI(); ++i)
        {
            topBoundaryPatch_(i, j, 0) = boundaryValue;
            facesK_(i, j, mesh_.nCellsK()) = topBoundaryPatch_(i, j, 0);
        } // end for i
    } //  end for j
}

template<class T>
void Field<T>::setBottomBoundary(const std::string &boundaryType, T boundaryValue)
{
    int i, j;

    if(boundaryType == "fixed")
        bottomBoundaryPatchId_ = FIXED;
    else if(boundaryType == "zeroGradient")
        bottomBoundaryPatchId_ = ZERO_GRADIENT;
    else if(boundaryType == "empty")
        bottomBoundaryPatchId_ = EMPTY;
    else
        Output::raiseException("Field", "setBottomBoundary", "unrecognized boundary condition type \"" + boundaryType + "\".");

    for(j = 0; j < bottomBoundaryPatch_.sizeJ(); ++j)
    {
        for(i = 0; i < bottomBoundaryPatch_.sizeI(); ++i)
        {
            bottomBoundaryPatch_(i, j, 0) = boundaryValue;
            facesK_(i, j, 0) = bottomBoundaryPatch_(i, j, 0);
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
    int maxI, j, k;

    switch(eastBoundaryPatchId_)
    {

    case FIXED:

        break;

    case ZERO_GRADIENT: case EMPTY:

        maxI = Array3D<T>::nI_ - 1;

        for(k = 0; k < eastBoundaryPatch_.sizeK(); ++k)
        {
            for(j = 0; j < westBoundaryPatch_.sizeJ(); ++j)
            {
                eastBoundaryPatch_(0, j, k) = Array3D<T>::operator ()(maxI, j, k);
                facesI_(mesh_.nCellsI(), j, k) = eastBoundaryPatch_(0, j, k);
            }
        }

        break;
    };
}

template<class T>
void Field<T>::setWestBoundaryField()
{
    int j, k;

    switch(westBoundaryPatchId_)
    {

    case FIXED:

        break;

    case ZERO_GRADIENT: case EMPTY:

        for(k = 0; k < westBoundaryPatch_.sizeK(); ++k)
        {
            for(j = 0; j < westBoundaryPatch_.sizeJ(); ++j)
            {
                westBoundaryPatch_(0, j, k) = Array3D<T>::operator ()(0, j, k);
                facesI_(0, j, k) = westBoundaryPatch_(0, j, k);
            }
        }

        break;
    };
}

template<class T>
void Field<T>::setNorthBoundaryField()
{
    int i, maxJ, k;

    switch(northBoundaryPatchId_)
    {

    case FIXED:

        break;

    case ZERO_GRADIENT: case EMPTY:

        maxJ = Array3D<T>::nJ_ - 1;

        for(k = 0; k < northBoundaryPatch_.sizeK(); ++k)
        {
            for(i = 0; i < northBoundaryPatch_.sizeI(); ++i)
            {
                northBoundaryPatch_(i, 0, k) = Array3D<T>::operator ()(i, maxJ, k);
                facesJ_(i, mesh_.nCellsJ(), k) = northBoundaryPatch_(i, 0, k);
            }
        }

        break;
    };
}

template<class T>
void Field<T>::setSouthBoundaryField()
{
    int i, k;

    switch(southBoundaryPatchId_)
    {
    case FIXED:

        break;

    case ZERO_GRADIENT: case EMPTY:

        for(k = 0; k < southBoundaryPatch_.sizeK(); ++k)
        {
            for(i = 0; i < southBoundaryPatch_.sizeI(); ++i)
            {
                southBoundaryPatch_(i, 0, k) = Array3D<T>::operator ()(i, 0, k);
                facesJ_(i, 0, k) = southBoundaryPatch_(i, 0, k);
            }
        }

        break;
    };
}

template<class T>
void Field<T>::setTopBoundaryField()
{
    int i, j, maxK;

    switch(topBoundaryPatchId_)
    {
    case FIXED:

        break;

    case ZERO_GRADIENT: case EMPTY:

        maxK = Array3D<T>::nK_ - 1;

        for(j = 0; j < topBoundaryPatch_.sizeJ(); ++j)
        {
            for(i = 0; i < topBoundaryPatch_.sizeI(); ++i)
            {
                topBoundaryPatch_(i, j, 0) = Array3D<T>::operator ()(i, j, maxK);
                facesK_(i, j, mesh_.nCellsK()) = topBoundaryPatch_(i, j, 0);
            }
        }

        break;
    };
}

template<class T>
void Field<T>::setBottomBoundaryField()
{
    int i, j;

    switch(bottomBoundaryPatchId_)
    {
    case FIXED:

        break;

    case ZERO_GRADIENT: case EMPTY:

        for(j = 0; j < bottomBoundaryPatch_.sizeJ(); ++j)
        {
            for(i = 0; i < bottomBoundaryPatch_.sizeI(); ++i)
            {
                bottomBoundaryPatch_(i, j, 0) = Array3D<T>::operator ()(i, j, 0);
                facesK_(i, j, 0) = bottomBoundaryPatch_(i, j, 0);
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
        switch(eastBoundaryPatchId_)
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
        switch(westBoundaryPatchId_)
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
        switch(northBoundaryPatchId_)
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
        switch(southBoundaryPatchId_)
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
        switch(topBoundaryPatchId_)
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
        switch(bottomBoundaryPatchId_)
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
