/**
 * @file    InitialConditions.h
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
 * This file contains the interface for class InitialConditions, which
 * is used for generating initial conditions on fields.
 */

#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H

#include <fstream>
#include <string>

#include <boost/property_tree/ptree.hpp>

#include "Field.h"

class InitialConditions
{

public:

    InitialConditions();

private:

    template <class T>
    void setInitialConditions(Field<T>& field);

    template <class T>
    void setSphere(Field<T>& field);

    template <class T>
    void setUniform(Field<T>& field);

    template <class T>
    void setBox(Field<T>& field);

    template <class T>
    void createUniform(T value, Field<T>& field);

    template <class T>
    void createSphere(double radius, Point3D center, T sphereInnerValue, Field<T>& field);

    template <class T>
    void createBox(double xLength, double yLength, double zLength, Point3D center, T boxInnerValue, Field<T>& field);

    template <class T>
    void smootheField();

    static double kCos(const Point3D &point);
    static double k

    boost::property_tree::ptree initialConditionsParameters_;
};

#include "InitialConditionsI.h"

#endif
