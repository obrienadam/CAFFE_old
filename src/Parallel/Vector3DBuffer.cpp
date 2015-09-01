#include "Vector3DBuffer.h"

Vector3DBuffer::Vector3DBuffer(int nBuffers, int bufferSize)
    :
      MINIMUM_BUFFER_SIZE_(bufferSize),
      sendBuffers_(nBuffers, std::vector<double>(MINIMUM_BUFFER_SIZE_)),
      recvBuffers_(nBuffers, std::vector<double>(MINIMUM_BUFFER_SIZE_)),
      sendBufferIter_(sendBuffers_.begin()),
      recvBufferIter_(recvBuffers_.begin()),
      vecPtrs_(nBuffers, nullptr),
      vecPtrIter_(vecPtrs_.begin())
{

}

const double *Vector3DBuffer::initSendBuffer(const std::vector<Vector3D> &vecs)
{
    if(sendBufferIter_ == sendBuffers_.end())
    {
        sendBuffers_.push_back(std::vector<double>(MINIMUM_BUFFER_SIZE_));
        sendBufferIter_ = std::prev(sendBuffers_.end(), -1);
    }

    std::vector<double> &sendBuffer = *sendBufferIter_;

    if(sendBuffer.size() < 3*vecs.size())
        sendBuffer.resize(3*vecs.size());

    for(int i = 0, j = 0, end = vecs.size(); i < end; ++i, j += 3)
    {
        sendBuffer[j] = vecs[i].x;
        sendBuffer[j + 1] = vecs[i].y;
        sendBuffer[j + 2] = vecs[i].z;
    }

    return sendBufferIter_++->data();
}

double *Vector3DBuffer::initRecvBuffer(std::vector<Vector3D> &vecs)
{
    if(recvBufferIter_ == recvBuffers_.end())
    {
        recvBuffers_.push_back(std::vector<double>(MINIMUM_BUFFER_SIZE_));
        recvBufferIter_ = std::prev(recvBuffers_.end(), -1);
        vecPtrs_.push_back(nullptr);
        vecPtrIter_ = std::prev(vecPtrs_.end(), -1);
    }

    if(recvBufferIter_->size() < 3*vecs.size())
        recvBufferIter_->resize(3*vecs.size());

    *vecPtrIter_ = &vecs;

    ++vecPtrIter_;

    return recvBufferIter_++->data();
}

double* Vector3DBuffer::initCollectiveBuffer(std::vector<Vector3D> &vecs)
{
    if(recvBufferIter_ == recvBuffers_.end())
    {
        recvBuffers_.push_back(std::vector<double>(MINIMUM_BUFFER_SIZE_));
        recvBufferIter_ = std::prev(recvBuffers_.end(), -1);
        vecPtrs_.push_back(nullptr);
        vecPtrIter_ = std::prev(vecPtrs_.end(), -1);
    }

    if(recvBufferIter_->size() < 3*vecs.size())
        recvBufferIter_->resize(3*vecs.size());

    *vecPtrIter_ = &vecs;
    std::vector<double> &buff = *recvBufferIter_;

    for(int i = 0, j = 0, end = vecs.size(); i < end; ++i, j += 3)
    {
        buff[j] = vecs[i].x;
        buff[j + 1] = vecs[i].y;
        buff[j + 2] = vecs[i].z;
    }

    ++vecPtrIter_;

    return recvBufferIter_++->data();
}

void Vector3DBuffer::freeBuffers()
{
    auto recvIt = recvBuffers_.begin();
    auto vecPtrIt = vecPtrs_.begin();

    //- Unload all recv buffers into the original vectors
    for(; recvIt != recvBufferIter_; ++recvIt, ++vecPtrIt)
    {
        std::vector<Vector3D> &vec = **vecPtrIt;
        std::vector<double> &buff = *recvIt;

        for(int i = 0, j = 0, end = vec.size(); i < end; ++i, j += 3)
        {
            vec[i].x = buff[j];
            vec[i].y = buff[j + 1];
            vec[i].z = buff[j + 2];
        }
    }

    //- Reset all iterators, essentially marking the buffers as available
    vecPtrs_.assign(vecPtrs_.size(), nullptr);
    sendBufferIter_ = sendBuffers_.begin();
    recvBufferIter_ = recvBuffers_.begin();
    vecPtrIter_ = vecPtrs_.begin();
}

void Vector3DBuffer::freeLastSendBuffer()
{
    if(sendBufferIter_ != sendBuffers_.begin())
        --sendBufferIter_;
}

void Vector3DBuffer::freeLastRecvBuffer()
{
    if(recvBufferIter_ != recvBuffers_.begin())
    {
        std::vector<double> &recvBuffer = *(--recvBufferIter_);
        std::vector<Vector3D> &vec = **(--vecPtrIter_);

        for(int i = 0, j = 0, end = vec.size(); i < end; ++i, j += 3)
        {
            vec[i].x = recvBuffer[j];
            vec[i].y = recvBuffer[j + 1];
            vec[i].z = recvBuffer[j + 2];
        }

        *vecPtrIter_ = nullptr;
    }
}
