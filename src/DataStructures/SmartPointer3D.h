#ifndef SMART_POINTER_3D_H
#define SMART_POINTER_3D_H

#include <iterator>

template <class T>
class SmartPointer3D
{

private:

    int nI_, nJ_, nK_, n_;
    T*** data_;

public:

    //- Constructors and destructors

    SmartPointer3D();
    SmartPointer3D(int nI, int nJ, int nK);
    ~SmartPointer3D();

    //- Memory management

    void allocate(int nI, int nJ, int nK);
    void deallocate();

    //- Return the container sizes

    int sizeI(){return nI_;}
    int sizeJ(){return nJ_;}
    int sizeK(){return nK_;}
    int size(){return n_;}

    inline T& operator()(int i, int j, int k);

    //- Iterators

    class iterator
    {

    private:

        int i_, j_, k_;
        T* dataPtr_;
        SmartPointer3D<T>* objectPtr_;

    public:

        typedef std::bidirectional_iterator_tag iterator_category;

        iterator();
        iterator(T* dataPtr,
                 SmartPointer3D<T>* objectPtr,
                 int i,
                 int j,
                 int k);

        iterator& operator++();

        T& operator*();

        bool operator!=(const iterator& rhs);

    };

    iterator begin();
    iterator end();

};

#include "SmartPointer3DI.h"

#endif
