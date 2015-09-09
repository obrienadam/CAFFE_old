/**
 * @file    IndexMap.cpp
 * @author  Adam O'Brien <obrienadam89@gmail.com>
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * https://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 *
 * This file contains the method implementations for class IndexMap.
 */

#include "IndexMap.h"
#include "FvScheme.h"
#include "Output.h"
#include "Parallel.h"

IndexMap::IndexMap()
{

}

void IndexMap::initialize(int nCellsI, int nCellsJ, int nCellsK)
{
    nCellsIThisProc_ = nCellsI;
    nCellsJThisProc_ = nCellsJ;
    nCellsKThisProc_ = nCellsK;

    globalIndices_.resize(nCellsIThisProc_, nCellsJThisProc_, nCellsKThisProc_);
    cellStatuses_.resize(nCellsIThisProc_, nCellsJThisProc_, nCellsKThisProc_);

    cellStatuses_.assign(ACTIVE);
    generateLocalIndices();
}

int IndexMap::operator ()(int i, int j, int k, int varSetNo)
{
    if(i >= 0 && i < nCellsIThisProc_
            && j >= 0 && j < nCellsJThisProc_
            && k >= 0 && k < nCellsKThisProc_)
        return globalIndices_(i, j, k) + varSetNo*nActiveGlobal_;

    if(adjProcNoPtr_)
    {
        const std::array<int, 6> &adjProcNo = *adjProcNoPtr_;

        if(i >= nCellsIThisProc_ && adjProcNo[HexaFvmMesh::EAST] != Parallel::PROC_NULL)
            return adjGlobalIndices_[0](i - nCellsIThisProc_, j, k) + varSetNo*nActiveGlobal_;

        else if(i < 0 && adjProcNo[HexaFvmMesh::WEST] != Parallel::PROC_NULL)
            return adjGlobalIndices_[1](adjGlobalIndices_[1].sizeI() + i, j, k) + varSetNo*nActiveGlobal_;

        else if(j >= nCellsJThisProc_ && adjProcNo[HexaFvmMesh::NORTH] != Parallel::PROC_NULL)
            return adjGlobalIndices_[2](i, j - nCellsJThisProc_, k) + varSetNo*nActiveGlobal_;

        else if(j < 0 && adjProcNo[HexaFvmMesh::SOUTH] != Parallel::PROC_NULL)
            return adjGlobalIndices_[3](i, adjGlobalIndices_[3].sizeJ() + j, k) + varSetNo*nActiveGlobal_;

        else if(k >= nCellsKThisProc_ && adjProcNo[HexaFvmMesh::TOP] != Parallel::PROC_NULL)
            return adjGlobalIndices_[4](i, j, k - nCellsKThisProc_) + varSetNo*nActiveGlobal_;

        else if(k < 0 && adjProcNo[HexaFvmMesh::BOTTOM] != Parallel::PROC_NULL)
            return adjGlobalIndices_[5](i, j, adjGlobalIndices_[5].sizeK() + k) + varSetNo*nActiveGlobal_;
    }

    return -1;
}

int IndexMap::nActiveGlobal() const
{
    return nActiveGlobal_;
}

int IndexMap::nActiveLocal() const
{
    return nActiveLocalThisProc_;
}

int IndexMap::nActive() const
{
    return nActiveGlobal_;
}

bool IndexMap::isActive(int i, int j, int k)
{
    return cellStatuses_(i, j, k) == ACTIVE;
}

bool IndexMap::isGhost(int i, int j, int k)
{
    return cellStatuses_(i, j, k) == GHOST;
}

bool IndexMap::isInactive(int i, int j, int k)
{
    return cellStatuses_(i, j, k) == INACTIVE;
}

void IndexMap::setActive(int i, int j, int k)
{
    cellStatuses_(i, j, k) = ACTIVE;
}

void IndexMap::setGhost(int i, int j, int k)
{
    cellStatuses_(i, j, k) = GHOST;
}

void IndexMap::setInactive(int i, int j, int k)
{
    cellStatuses_(i, j, k) = INACTIVE;
}

void IndexMap::generateLocalIndices()
{
    nActiveLocalThisProc_ = 0;

    for(int i = 0; i < cellStatuses_.size(); ++i)
    {
        switch (cellStatuses_[i])
        {
        case ACTIVE: case GHOST:
            globalIndices_[i] = nActiveLocalThisProc_;
            ++nActiveLocalThisProc_;
            break;

        case INACTIVE:
            globalIndices_[i] = -1;
            break;
        }
    }

    nActiveGlobal_ = nActiveLocalThisProc_;
}

void IndexMap::generateGlobalIndices(std::shared_ptr< std::array<int, 6> > adjProcNoPtr)
{
    if(!adjProcNoPtr)
        return;

    int lowerGlobalIndex = 0;

    adjProcNoPtr_ = adjProcNoPtr;
    nActiveGlobal_ = 0;
    maxNLocal_ = 0;

    nActiveLocal_.resize(Parallel::nProcesses());
    nCellsILocal_.resize(Parallel::nProcesses());
    nCellsJLocal_.resize(Parallel::nProcesses());
    nCellsKLocal_.resize(Parallel::nProcesses());

    Parallel::allGather(nActiveLocalThisProc_, nActiveLocal_);
    Parallel::allGather(nCellsIThisProc_, nCellsILocal_);
    Parallel::allGather(nCellsJThisProc_, nCellsJLocal_);
    Parallel::allGather(nCellsKThisProc_, nCellsKLocal_);

    for(int i = 0; i < Parallel::nProcesses(); ++i)
    {
        if(i < Parallel::processNo())
            lowerGlobalIndex += nActiveLocal_[i];

        nActiveGlobal_ += nActiveLocal_[i];

        if(maxNLocal_ < nActiveLocal_[i])
            maxNLocal_ = nActiveLocal_[i];
    }

    globalIndices_.add(lowerGlobalIndex);

    std::array<int, 6> adjProcNo = *adjProcNoPtr_;

    for(int i = 0; i < 6; ++i)
    {
        if(adjProcNo[i] == Parallel::PROC_NULL)
            continue;

        adjGlobalIndices_[i].resize(nCellsILocal_[adjProcNo[i]], nCellsJLocal_[adjProcNo[i]], nCellsKLocal_[adjProcNo[i]]);

        Parallel::iSend(Parallel::processNo(), adjProcNo[i], i, globalIndices_);
        Parallel::iRecv(adjProcNo[i], Parallel::processNo(), i%2 == 0 ? i + 1 : i - 1, adjGlobalIndices_[i]);
    }

    Parallel::waitAll();
}
