#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <utility>

#include "Point3D.h"

class Geometry
{
private:

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
    static bool isInsideHexahedron(const Point3D& point, const Point3D* points);

    //- Tetrahedron methods
    static double computeTetrahedronVolume(Point3D* points);
    static Point3D computeTetrahedronCentroid(Point3D* points);

    //- Basic geometric methods
    virtual double surfaceArea() = 0;
    virtual double volume() = 0;

    //- Intersections
    virtual Point3D nearestIntersect(const Point3D& point) = 0;
    virtual std::pair<Point3D, Point3D> lineIntersect(const Point3D& pt1, const Point3D& pt2) = 0;

    //- Tests
    virtual bool isInside(const Point3D& point) = 0;

    //- Operators and movement
    virtual Geometry& operator+=(const Vector3D& translationVec) = 0;
    virtual Geometry& operator-=(const Vector3D& translationVec) = 0;
    virtual Geometry& rotate(const double radians) = 0;
    virtual Geometry& scale(const double a) = 0;
};

#endif
