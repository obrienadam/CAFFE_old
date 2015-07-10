/**
 * @file    IndexMap.h
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
 * This file contains the interface for class IndexMap, which is used
 * to produce a map functor that converts between i,j,k indices of a mesh
 * to a global index to be used for assembling matrices.
 */

#ifndef INDEX_MAP_H
#define INDEX_MAP_H

#include "Field.h"

enum Ordering{I_INDEX, J_INDEX, K_INDEX};

class IndexMap
{
private:

    int nI_, nJ_, nK_, nActive_;

    Field<int> indexField_;

    Ordering index1_, index2_, index3_;

public:

    IndexMap(Ordering index1 = I_INDEX, Ordering index2 = J_INDEX, Ordering index3 = K_INDEX);

    int operator()(int i, int j, int k, int varSetNo);
    void generateMap(Field<int>& cellStatus);
    int nActive(){ return nActive_; }
};

#endif
