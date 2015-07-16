/**
 * @file    IndexMap.cc
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

IndexMap::IndexMap()
{

}

void IndexMap::initialize(int nCellsI, int nCellsJ, int nCellsK)
{
    nCellsI_ = nCellsI;
    nCellsJ_ = nCellsJ;
    nCellsK_ = nCellsK;
    nActive_ = nCellsI_*nCellsJ_*nCellsK_;
    globalIndices_.allocate(nCellsI_, nCellsJ_, nCellsK_);
    cellStatuses_.allocate(nCellsI_, nCellsJ_, nCellsK_);

    for(int i = 0; i < nActive_; ++i)
    {
        globalIndices_(i) = i;
        cellStatuses_(i) = ACTIVE;
    }
}

int IndexMap::operator ()(int i, int j, int k, int varSetNo)
{
    if(i >= 0 && i < nCellsI_
            && j >= 0 && j < nCellsJ_
            && k >= 0 && k < nCellsK_)
    {
        return globalIndices_(i, j, k) + varSetNo*nActive_;
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

void IndexMap::generateGlobalIndices()
{
    nActive_ = 0;
    for(int i = 0; i < cellStatuses_.size(); ++i)
    {
        if(cellStatuses_(i) == ACTIVE || cellStatuses_(i) == GHOST)
        {
            globalIndices_(i) = nActive_;
            ++nActive_;
        }
        else
        {
            globalIndices_(i) = -1;
        }
    }
}
