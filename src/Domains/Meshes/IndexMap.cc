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

int IndexMap::operator ()(int i, int j, int k, int varSetNo)
{
    if(isActive_(i, j, k))
    {
        return indexField_(i, j, k) + varSetNo*nActive_;
    }
    else
    {
        Output::raiseException("IndexMap", "operator()", "attempted to access an index belonging to a cell that is either inactive or an interpolation cell.");
    }

    return -1;
}

void IndexMap::generateMap(Field<int> &cellStatus)
{
    int i, j, k;

    nActive_ = 0;

    nI_ = cellStatus.sizeI();
    nJ_ = cellStatus.sizeJ();
    nK_ = cellStatus.sizeK();

    indexField_.allocate(nI_, nJ_, nK_);
    isActive_.allocate(nI_, nJ_, nK_);

    for(k = 0; k < nK_; ++k)
    {
        for(j = 0; j < nJ_; ++j)
        {
            for(i = 0; i < nI_; ++i)
            {
                if(cellStatus(i, j, k) == ACTIVE)
                {
                    indexField_(i, j, k) = nActive_;
                    isActive_(i, j, k) = true;
                    ++nActive_;
                }
                else
                {
                    indexField_(i, j, k) = 0;
                    isActive_(i, j, k) = false;
                }
            }
        }
    }
}
