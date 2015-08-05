#ifndef ARRAY_3D_H
#define ARRAY_3D_H

#include <iterator>

template <class T>
class Array3D
{

protected:

    int nI_, nJ_, nK_, n_;
    int nInJ_, nInK_, nJnK_;
    T* data_;

public:

    //- Constructors and destructors

    Array3D();
    Array3D(int nI, int nJ, int nK);
    ~Array3D();

    //- Memory management

    virtual void allocate(int nI, int nJ, int nK);
    void deallocate();

    //- Return the container sizes

    int sizeI() const {return nI_;}
    int sizeJ() const {return nJ_;}
    int sizeK() const {return nK_;}
    int size() const {return n_;}
    int elementNo(int i, int j, int k) const { return k*nInJ_ + j*nI_ + i; }

    virtual inline T& operator()(int i, int j, int k);
    virtual inline const T& operator()(int i, int j, int k) const;
   	inline T& operator()(int k);
    inline const T& operator()(int k) const;
    T* data(){ return data_; }

    virtual void setAll(const T &value);

    //- Iterators

    class iterator
    {

    private:

        int k_;
        T* dataPtr_;
        Array3D<T>* objectPtr_;

    public:

        typedef std::bidirectional_iterator_tag iterator_category;

        iterator();
        iterator(T* dataPtr,
                 Array3D<T>* objectPtr,
                 int k);

        iterator& operator++();

        T& operator*();

        bool operator!=(const iterator& rhs);
    };

    iterator begin();
    iterator end();

};

#include "Array3DI.h"

#endif
