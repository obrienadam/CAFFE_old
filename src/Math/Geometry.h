#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "Point3D.h"

class Geometry
{


public:

    //- 2-Dimensional

    static double computeQuadrilateralArea(Point3D* points);
    static Point3D computeQuadrilateralCentroid(Point3D* points);
    static Vector3D computeQuadrilateralNormal(Point3D* points);
    static bool checkQuadrilateralIsPlanar(Point3D* points);

    //- 3-Dimensional

    //- Hexahedron methods

    static double computeHexahedronVolume(Point3D* points);
    static Point3D computeHexahedronCentroid(Point3D* points);
    static bool checkHexahedronSurfacesIsPlanar(Point3D* points);

    //- Tetrahedron methods

    static double computeTetrahedronVolume(Point3D* points);
    static Point3D computeTetrahedronCentroid(Point3D* points);

};

#endif
