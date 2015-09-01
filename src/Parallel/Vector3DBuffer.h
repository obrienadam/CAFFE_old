#ifndef VECTOR_3D_BUFFER_H
#define VECTOR_3D_BUFFER_H

#include <vector>

#include "Vector3D.h"

class Vector3DBuffer
{
public:

    Vector3DBuffer(int nBuffers = 6, int bufferSize = 1000);

    const double* initSendBuffer(const std::vector<Vector3D> &vecs);
    double* initRecvBuffer(std::vector<Vector3D> &vecs);
    double* initCollectiveBuffer(std::vector<Vector3D> &vecs);

    //- This method is typically used upon terminating non-blocking communications
    void freeBuffers();

    //- These methods are typically used for blocking communications
    void freeLastSendBuffer();
    void freeLastRecvBuffer();

private:

    std::vector< std::vector<double> > sendBuffers_, recvBuffers_;
    std::vector< std::vector<double> >::iterator sendBufferIter_, recvBufferIter_;

    std::vector< std::vector<Vector3D>* > vecPtrs_;
    std::vector< std::vector<Vector3D>* >::iterator vecPtrIter_;

    const int MINIMUM_BUFFER_SIZE_;

};

#endif
