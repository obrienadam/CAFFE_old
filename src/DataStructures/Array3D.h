#ifndef ARRAY_3D_H
#define ARRAY_3D_H

#include <vector>

template <class T>
class Array3D : public std::vector<T>
{
public:

    Array3D(int sizeI = 0, int sizeJ = 0, int sizeK = 0);

    void resize(int sizeI, int sizeJ, int sizeK);
    void clear();
    void assign(const T& val);
    void add(const T& val);

    int sizeI() const { return sizeI_; }
    int sizeJ() const { return sizeJ_; }
    int sizeK() const { return sizeK_; }
    int upperI() const { return sizeI_ - 1; }
    int upperJ() const { return sizeJ_ - 1; }
    int upperK() const { return sizeK_ - 1; }

    T& operator()(int i, int j, int k);
    const T& operator ()(int i, int j, int k) const;

protected:

    void setSizes(int sizeI, int sizeJ, int sizeK);

    int sizeI_, sizeJ_, sizeK_, sizeIJ_;
};

#include "Array3D.tpp"

#endif
