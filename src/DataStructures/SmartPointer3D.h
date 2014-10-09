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
        T* ptr_;

    public:

        typedef std::bidirectional_iterator_tag iterator_category;

        iterator();
        iterator(T* ptr);

        iterator operator++();
        T& operator*();
        T* operator->();
        bool operator==(const iterator& rhs);
        bool operator!=(const iterator& rhs);

    };

    iterator begin();
    iterator end();

};

#include "SmartPointer3DI.h"

#endif
