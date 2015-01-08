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
void Field<T>::resize(int nI, int nJ, int nK)
{

    Array3D<T>::allocate(nI, nJ, nK);

    if (type == CONSERVED)
    {

        faceFluxesI_.allocate(nI + 1, nJ, nK);
        faceFluxesJ_.allocate(nI, nJ + 1, nK);
        faceFluxesK_.allocate(nI, nJ, nK + 1);

    }

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
void Field<T>::setAllBoundaries(BoundaryPatch eastBoundaryType,
                                BoundaryPatch westBoundaryType,
                                BoundaryPatch northBoundaryType,
                                BoundaryPatch southBoundaryType,
                                BoundaryPatch topBoundaryType,
                                BoundaryPatch bottomBoundaryType)
{

    setEastBoundary(eastBoundaryType);
    setWestBoundary(westBoundaryType);
    setNorthBoundary(northBoundaryType);
    setSouthBoundary(southBoundaryType);
    setTopBoundary(topBoundaryType);
    setBottomBoundary(bottomBoundaryType);

    Output::printToScreen("Field", "boundaries for field \"" + name + "\" have been set.");

}

template<class T>
void Field<T>::setEastBoundary(BoundaryPatch boundaryType)
{

    int j, k;

    eastBoundaryPatch_.allocate(Array3D<T>::nJ_, Array3D<T>::nK_);

    for(k = 0; k < Array3D<T>::nK_; ++k)
    {

        for(j = 0; j < Array3D<T>::nJ_; ++j)
        {

            eastBoundaryPatch_(j, k) = boundaryType;

        } // end for j
    } //  end for k

}

template<class T>
void Field<T>::setWestBoundary(BoundaryPatch boundaryType)
{

    int j, k;

    westBoundaryPatch_.allocate(Array3D<T>::nJ_, Array3D<T>::nK_);

    for(k = 0; k < Array3D<T>::nK_; ++k)
    {

        for(j = 0; j < Array3D<T>::nJ_; ++j)
        {

            westBoundaryPatch_(j, k) = boundaryType;

        } // end for j
    } //  end for k

}

template<class T>
void Field<T>::setNorthBoundary(BoundaryPatch boundaryType)
{

    int i, k;

    northBoundaryPatch_.allocate(Array3D<T>::nI_, Array3D<T>::nK_);

    for(k = 0; k < Array3D<T>::nK_; ++k)
    {

        for(i = 0; i < Array3D<T>::nI_; ++i)
        {

            northBoundaryPatch_(i, k) = boundaryType;

        } // end for i
    } //  end for k

}

template<class T>
void Field<T>::setSouthBoundary(BoundaryPatch boundaryType)
{

    int i, k;

    southBoundaryPatch_.allocate(Array3D<T>::nI_, Array3D<T>::nK_);

    for(k = 0; k < Array3D<T>::nK_; ++k)
    {

        for(i = 0; i < Array3D<T>::nI_; ++i)
        {

            southBoundaryPatch_(i, k) = boundaryType;

        } // end for i
    } //  end for k

}

template<class T>
void Field<T>::setTopBoundary(BoundaryPatch boundaryType)
{

    int i, j;

    topBoundaryPatch_.allocate(Array3D<T>::nI_, Array3D<T>::nJ_);

    for(j = 0; j < Array3D<T>::nJ_; ++j)
    {

        for(i = 0; i < Array3D<T>::nI_; ++i)
        {

            topBoundaryPatch_(i, j) = boundaryType;

        } // end for i
    } //  end for j

}

template<class T>
void Field<T>::setBottomBoundary(BoundaryPatch boundaryType)
{

    int i, j;

    bottomBoundaryPatch_.allocate(Array3D<T>::nI_, Array3D<T>::nJ_);

    for(j = 0; j < Array3D<T>::nJ_; ++j)
    {

        for(i = 0; i < Array3D<T>::nI_; ++i)
        {

            bottomBoundaryPatch_(i, j) = boundaryType;

        } // end for i
    } //  end for j

}
