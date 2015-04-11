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

#include <vector>
#include <string>

#include "Field.h"
#include "HexaFvmMesh.h"

class FvScheme
{
protected:

    std::string conservedFieldName_;
    HexaFvmMesh* meshPtr_;

public:

    FvScheme();

    virtual void initialize(HexaFvmMesh& mesh, std::string conservedFieldName = "phi");
    virtual int nConservedVariables() = 0;

    virtual void discretize(std::vector<double>& timeDerivatives_) = 0;
    virtual void updateSolution(std::vector<double>& timeDerivatives_) = 0;
    virtual double computeUpdateNorm(std::vector<double>& timeDerivatives_);


    /**
     * @brief This method is used for computing a weighted averaging coefficient based on a specified criteria.
     * @return An interpolation factor.
     */
    double getAlpha(int i, int j, int k, int direction);
};

#endif
