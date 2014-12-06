#include <math.h>

#include "Geometry.h"

double Geometry::computeQuadrilateralArea(Point3D* points)
{

    Vector3D r01, r02, r03;

    r01 = points[1] - points[0];
    r02 = points[2] - points[0];
    r03 = points[3] - points[0];

    return 0.5*(crossProduct(r02, r01).mag() + crossProduct(r03, r02).mag());

}

double Geometry::computeHexahedronVolume(Point3D* points)
{

    // This uses the method of calculating hexahedron volumes depicted in
    // Computational Fluid Mechanics and Heat Transfer, Third Edition, by
    // Richard H. Pletcher, John C. Tannehill and Gale Anderson

    Vector3D r16, r25, r57, r64, r27, r36;
    Vector3D s1265, s5674, s2376;

    r16 = points[6] - points[1];
    r25 = points[5] - points[2];
    s1265 = 0.5*crossProduct(r16, r25);

    r57 = points[7] - points[5];
    r64 = points[4] - points[6];
    s5674 = 0.5*crossProduct(r57, r64);

    r27 = points[7] - points[2];
    r36 = points[6] - points[3];
    s2376 = 0.5*crossProduct(r27, r36);

    return fabs(dotProduct(s1265 + s5674 + s2376, points[6] - points[0])/3.);

}
