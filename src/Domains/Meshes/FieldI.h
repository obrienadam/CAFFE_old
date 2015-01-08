#include "Field.h"

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
