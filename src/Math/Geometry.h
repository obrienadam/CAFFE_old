#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "Point3D.h"

class Geometry
{


public:

    //- 2-Dimensional

    static double computeQuadrilateralArea(Point3D* points);

    //- 3-Dimensional

    static double computeHexahedronVolume(Point3D* points);

};

#endif
