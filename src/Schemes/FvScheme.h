/**
 * @file    FvScheme.h
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
 * This file contains the interface for class FvScheme, which is an
 * abstract interface class used for deriving specific finite volume
 * spatial schemes.
 */

#ifndef FV_SCHEME_H
#define FV_SCHEME_H

#include "Field.h"

class FvScheme
{
public:

    enum GradientEvaluationMethod{LEAST_SQUARES, DIVERGENCE_THEOREM};
    enum InterpolationMethod{VOLUME_WEIGHTED, DISTANCE_WEIGHTED, NON_WEIGHTED};

    FvScheme();
    ~FvScheme();

    template <class T>
    static void interpolateInteriorFaces(InterpolationMethod method, Field<T>& field);

    template <class T, class GRAD_T>
    static void extrapolateInteriorFaces(GradientEvaluationMethod method, Field<T>& field, Field<GRAD_T>& gradField);
};

#include "FvSchemeI.h"

#endif
