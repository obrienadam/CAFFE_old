#ifndef PARALLEL_BUFFERS_H
#define PARALLEL_BUFFERS_H

#include <vector>
#include "Vector3D.h"

class ParallelBuffers
{
public:

    ParallelBuffers();

    double* initRecvBuffer(std::vector<double> &doubles);
    double* initRecvBuffer(std::vector<Vector3D> &vecs);
    double* initCollectiveBuffer(bool copyToRecvBuff, std::vector<Vector3D> &vecs);

    const double* initSendBuffer(const std::vector<Vector3D> &vecs);

    //- These should only be called if all outstanding MPI requests are resolved
    void freeAllBuffers();
    void freeAllDoubleBuffers();
    void freeAllVector3DBuffers();

    //- These are primarily used for blocking comms, collective comms etc
    void freeLastVector3DSendBuffer();
    void freeLastVector3DRecvBuffer();

private:

    static std::vector<double>& getNextAvailableBuffer(std::vector< std::vector<double> > &buffer,
                                                       std::vector< std::vector<double> >::iterator &bufferIter);

    static void copyVec3DVecToBuff(const std::vector<Vector3D> &vecs, double *buff);
    static void copyBuffToVec3DVec(const double *buff, std::vector<Vector3D> &vecs);

    std::vector< std::vector<double> > doubleRecvBuffs_;
    std::vector< std::vector<double> > vec3DSendBuffs_, vec3DRecvBuffs_;

    std::vector< std::vector<double> >::iterator nextDoubleRecvBuff_;
    std::vector< std::vector<double> >::iterator nextVec3DSendBuff_, nextVec3DRecvBuff_;

    std::vector< std::vector<double>* > recvDoublePtrs_;
    std::vector< std::vector<Vector3D>* > recvVec3DPtrs_;
};

#endif
