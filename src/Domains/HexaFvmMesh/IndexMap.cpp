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
    :
      lowerGlobalIndex_(0),
      upperGlobalIndex_(0),
      gatheredNActiveLocal_(Parallel::nProcesses())
{

}

void IndexMap::initialize(int nCellsI, int nCellsJ, int nCellsK)
{
    nCellsI_ = nCellsI;
    nCellsJ_ = nCellsJ;
    nCellsK_ = nCellsK;
    localIndices_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    cellStatuses_.allocate(nCellsI_, nCellsJ_, nCellsK_);

    for(int i = 0; i < nActiveLocal_; ++i)
    {
        localIndices_(i) = i;
        cellStatuses_(i) = ACTIVE;
    }

    generateIndices();
}

int IndexMap::operator ()(int i, int j, int k, int varSetNo)
{
    if(i >= 0 && i < nCellsI_
            && j >= 0 && j < nCellsJ_
            && k >= 0 && k < nCellsK_)
    {
        return lowerGlobalIndex_ + localIndices_(i, j, k) + varSetNo*nActiveLocal_;
    }

    return -1;
}

bool IndexMap::isActive(int i, int j, int k)
{
    if(cellStatuses_(i, j, k) == ACTIVE)
        return true;
    else
        return false;
}

bool IndexMap::isGhost(int i, int j, int k)
{
    if(cellStatuses_(i, j, k) == GHOST)
        return true;
    else
        return false;
}

bool IndexMap::isInactive(int i, int j, int k)
{
    if(cellStatuses_(i, j, k) == INACTIVE)
        return true;
    else
        return false;
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

void IndexMap::generateIndices()
{
    nActiveLocal_ = 0;
    nActiveGlobal_ = 0;

    for(int i = 0; i < cellStatuses_.size(); ++i)
    {
        switch (cellStatuses_(i))
        {
        case ACTIVE: case GHOST:
            localIndices_(i) = nActiveLocal_;
            ++nActiveLocal_;
            break;

        case INACTIVE:
            localIndices_(i) = -1;
        }
    }

    lowerGlobalIndex_ = 0;
    Parallel::allGather(nActiveLocal_, gatheredNActiveLocal_);

    for(int i = 0; i < Parallel::nProcesses(); ++i)
    {
        if(i < Parallel::processNo())
            lowerGlobalIndex_ += gatheredNActiveLocal_[i];

        nActiveGlobal_ += gatheredNActiveLocal_[i];
    }

    upperGlobalIndex_ = lowerGlobalIndex_ + nActiveLocal_ - 1;
}
