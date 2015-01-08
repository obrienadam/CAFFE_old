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

    Output::printToScreen("Field", "boundaries for field \"" + name + "\" have been set.");

}

template<class T>
void Field<T>::setEastBoundary(BoundaryPatch boundaryType, T boundaryValue)
{

    int j, k;

    eastBoundaryPatch_.allocate(Array3D<T>::nJ_, Array3D<T>::nK_);
    eastBoundaryField_.allocate(Array3D<T>::nJ_, Array3D<T>::nK_);

    for(k = 0; k < Array3D<T>::nK_; ++k)
    {

        for(j = 0; j < Array3D<T>::nJ_; ++j)
        {

            eastBoundaryPatch_(j, k) = boundaryType;
            eastBoundaryField_(j, k) = boundaryValue;

        } // end for j
    } //  end for k

}

template<class T>
void Field<T>::setWestBoundary(BoundaryPatch boundaryType, T boundaryValue)
{

    int j, k;

    westBoundaryPatch_.allocate(Array3D<T>::nJ_, Array3D<T>::nK_);
    westBoundaryField_.allocate(Array3D<T>::nJ_, Array3D<T>::nK_);

    for(k = 0; k < Array3D<T>::nK_; ++k)
    {

        for(j = 0; j < Array3D<T>::nJ_; ++j)
        {

            westBoundaryPatch_(j, k) = boundaryType;
            westBoundaryField_(j, k) = boundaryValue;

        } // end for j
    } //  end for k

}

template<class T>
void Field<T>::setNorthBoundary(BoundaryPatch boundaryType, T boundaryValue)
{

    int i, k;

    northBoundaryPatch_.allocate(Array3D<T>::nI_, Array3D<T>::nK_);
    northBoundaryField_.allocate(Array3D<T>::nI_, Array3D<T>::nK_);

    for(k = 0; k < Array3D<T>::nK_; ++k)
    {

        for(i = 0; i < Array3D<T>::nI_; ++i)
        {

            northBoundaryPatch_(i, k) = boundaryType;
            northBoundaryField_(i, k) = boundaryValue;

        } // end for i
    } //  end for k

}

template<class T>
void Field<T>::setSouthBoundary(BoundaryPatch boundaryType, T boundaryValue)
{

    int i, k;

    southBoundaryPatch_.allocate(Array3D<T>::nI_, Array3D<T>::nK_);
    southBoundaryField_.allocate(Array3D<T>::nI_, Array3D<T>::nK_);

    for(k = 0; k < Array3D<T>::nK_; ++k)
    {

        for(i = 0; i < Array3D<T>::nI_; ++i)
        {

            southBoundaryPatch_(i, k) = boundaryType;
            southBoundaryField_(i, k) = boundaryValue;

        } // end for i
    } //  end for k

}

template<class T>
void Field<T>::setTopBoundary(BoundaryPatch boundaryType, T boundaryValue)
{

    int i, j;

    topBoundaryPatch_.allocate(Array3D<T>::nI_, Array3D<T>::nJ_);
    topBoundaryField_.allocate(Array3D<T>::nI_, Array3D<T>::nJ_);

    for(j = 0; j < Array3D<T>::nJ_; ++j)
    {

        for(i = 0; i < Array3D<T>::nI_; ++i)
        {

            topBoundaryPatch_(i, j) = boundaryType;
            topBoundaryField_(i, j) = boundaryValue;

        } // end for i
    } //  end for j

}

template<class T>
void Field<T>::setBottomBoundary(BoundaryPatch boundaryType, T boundaryValue)
{

    int i, j;

    bottomBoundaryPatch_.allocate(Array3D<T>::nI_, Array3D<T>::nJ_);
    bottomBoundaryField_.allocate(Array3D<T>::nI_, Array3D<T>::nJ_);

    for(j = 0; j < Array3D<T>::nJ_; ++j)
    {

        for(i = 0; i < Array3D<T>::nI_; ++i)
        {

            bottomBoundaryPatch_(i, j) = boundaryType;
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

            switch(eastBoundaryPatch_(j, k))
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

            switch(westBoundaryPatch_(j, k))
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

            switch(northBoundaryPatch_(i, k))
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

            switch(southBoundaryPatch_(i, k))
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

            switch(topBoundaryPatch_(i, j))
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

            switch(bottomBoundaryPatch_(i, j))
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
