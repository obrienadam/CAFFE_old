#include "ParallelBuffers.h"

ParallelBuffers::ParallelBuffers()
{
    using namespace std;

    doubleRecvBuffs_.resize(1, vector<double>(100));
    vec3DRecvBuffs_.resize(1, vector<double>(300));
    vec3DSendBuffs_.resize(1, vector<double>(300));

    nextDoubleRecvBuff_ = doubleRecvBuffs_.begin();
    nextVec3DRecvBuff_ = vec3DRecvBuffs_.begin();
    nextVec3DSendBuff_ = vec3DSendBuffs_.begin();
}

double* ParallelBuffers::initRecvBuffer(std::vector<double> &doubles)
{
    using namespace std;

    vector<double> &recvBuffer = getNextAvailableBuffer(doubleRecvBuffs_, nextDoubleRecvBuff_);

    if(recvBuffer.size() < doubles.size())
        recvBuffer.resize(doubles.size());

    recvDoublePtrs_.push_back(&doubles);

    return recvBuffer.data();
}

double* ParallelBuffers::initRecvBuffer(std::vector<Vector3D> &vecs)
{
    using namespace std;

    vector<double> &recvBuffer = getNextAvailableBuffer(vec3DRecvBuffs_, nextVec3DRecvBuff_);

    if(recvBuffer.size() < 3*vecs.size())
        recvBuffer.resize(3*vecs.size());

    recvVec3DPtrs_.push_back(&vecs);

    return recvBuffer.data();
}

double* ParallelBuffers::initCollectiveBuffer(bool copyToRecvBuff, std::vector<Vector3D> &vecs)
{
    using namespace std;

    double* recvBuff = initRecvBuffer(vecs);

    if(copyToRecvBuff)
        copyVec3DVecToBuff(vecs, recvBuff);

    return recvBuff;
}

const double* ParallelBuffers::initSendBuffer(const std::vector<Vector3D> &vecs)
{
    using namespace std;

    vector<double> &sendBuff = getNextAvailableBuffer(vec3DSendBuffs_, nextVec3DSendBuff_);

    if(sendBuff.size() < 3*vecs.size())
        sendBuff.resize(3*vecs.size());

    copyVec3DVecToBuff(vecs, sendBuff.data());

    return sendBuff.data();
}

void ParallelBuffers::freeAllBuffers()
{
    freeAllDoubleBuffers();
    freeAllVector3DBuffers();
}

void ParallelBuffers::freeAllDoubleBuffers()
{
    using namespace std;

    if(recvDoublePtrs_.size() == 0)
        return;

    auto it = doubleRecvBuffs_.begin();
    for(int doubleNo = 0, end = recvDoublePtrs_.size(); doubleNo < end; ++doubleNo, ++it)
    {
        const vector<double> &recvBuff = *it;
        vector<double> &doubleVec = *recvDoublePtrs_[doubleNo];

        for(int i = 0, end = doubleVec.size(); i < end; ++i)
        {
            doubleVec[i] = recvBuff[i];
        }
    }

    recvDoublePtrs_.clear();
    nextDoubleRecvBuff_ = doubleRecvBuffs_.begin();
}

void ParallelBuffers::freeAllVector3DBuffers()
{
    using namespace std;

    auto it = vec3DRecvBuffs_.begin();
    for(int vec3DNo = 0, end = recvVec3DPtrs_.size(); vec3DNo < end; ++vec3DNo, ++it)
    {
        auto &recvBuff = *it;
        vector<Vector3D> &vec3DVec = *(recvVec3DPtrs_[vec3DNo]);

        copyBuffToVec3DVec(recvBuff.data(), vec3DVec);
    }

    recvVec3DPtrs_.clear();
    nextVec3DRecvBuff_ = vec3DRecvBuffs_.begin();
    nextVec3DSendBuff_ = vec3DSendBuffs_.begin();
}

void ParallelBuffers::freeLastVector3DSendBuffer()
{
    using namespace std;

    if(nextVec3DSendBuff_ == vec3DSendBuffs_.begin())
        return;

    --nextVec3DSendBuff_;
}

void ParallelBuffers::freeLastVector3DRecvBuffer()
{
    using namespace std;

    if(recvVec3DPtrs_.empty())
        return;

    const vector<double> &recvBuff = *(nextVec3DRecvBuff_ - 1);
    vector<Vector3D> &vec3DVec = *recvVec3DPtrs_.back();

    copyBuffToVec3DVec(recvBuff.data(), vec3DVec);

    --nextVec3DRecvBuff_;
    recvVec3DPtrs_.pop_back();
}

//*********************** Private helper methods **************************************

std::vector<double>& ParallelBuffers::getNextAvailableBuffer(std::vector<std::vector<double> > &buffer, std::vector<std::vector<double> >::iterator &bufferIter)
{
    using namespace std;

    ++bufferIter;
    if(bufferIter == buffer.end())
    {
        buffer.push_back(vector<double>(100));
        bufferIter = buffer.end() - 1;
    }

    return *(bufferIter - 1);
}

void ParallelBuffers::copyVec3DVecToBuff(const std::vector<Vector3D> &vecs, double *buff)
{
    for(int i = 0, j = 0, end = vecs.size(); i < end; ++i, j += 3)
    {
        buff[j] = vecs[i].x;
        buff[j + 1] = vecs[i].y;
        buff[j + 2] = vecs[i].z;
    }
}

void ParallelBuffers::copyBuffToVec3DVec(const double *buff, std::vector<Vector3D> &vecs)
{
    for(int i = 0, j = 0, end = vecs.size(); i < end; ++i, j += 3)
    {
        vecs[i].x = buff[j];
        vecs[i].y = buff[j + 1];
        vecs[i].z = buff[j + 2];
    }
}
