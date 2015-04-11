#ifndef ARRAY_2D_H
#define ARRAY_2D_H

#include <iterator>

template <class T>
class Array2D
{

protected:

    int nI_, nJ_, n_;
    T* data_;

public:

    //- Constructors and destructors

    Array2D();
    Array2D(int nI, int nJ);
    ~Array2D();

    //- Memory management

    void allocate(int nI, int nJ);
    void deallocate();

    //- Return the container sizes

    int sizeI(){return nI_;}
    int sizeJ(){return nJ_;}
    int size(){return n_;}

    inline T& operator()(int i, int j);

    //- Iterators

    class iterator
    {

    private:

        int i_, j_;
        T* dataPtr_;
        Array2D<T>* objectPtr_;

    public:

        typedef std::bidirectional_iterator_tag iterator_category;

        iterator();
        iterator(T* dataPtr,
                 Array2D<T>* objectPtr,
                 int i,
                 int j);

        iterator& operator++();

        T& operator*();

        bool operator!=(const iterator& rhs);

    };

    iterator begin();
    iterator end();
};

#include "Array2DI.h"

#endif
