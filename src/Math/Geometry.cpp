#include <math.h>

#include "Geometry.h"

double Geometry::computeQuadrilateralArea(Point3D* points)
{
    Vector3D r01 = points[1] - points[0];
    Vector3D r02 = points[2] - points[0];
    Vector3D r03 = points[3] - points[0];

    return 0.5*(cross(r02, r01).mag() + cross(r03, r02).mag());
}

Point3D Geometry::computeQuadrilateralCentroid(Point3D* points)
{
    Point3D diag1Midpoint = points[0] + 0.5*(points[2] - points[0]);
    Point3D diag2Midpoint = points[1] + 0.5*(points[3] - points[1]);

    return diag1Midpoint + 0.5*(diag2Midpoint - diag1Midpoint);
}

Vector3D Geometry::computeQuadrilateralNormal(Point3D* points)
{
    return 0.5*(cross(points[1] - points[0], points[2] - points[0]) + cross(points[1] - points[0], points[2] - points[0]));
}

bool Geometry::checkQuadrilateralIsPlanar(Point3D *points)
{
    // Check to see if a quadrilateral is planary by comparing the normals formed by
    // the triangles sharing the diagonal

    const double TOLER(1e-12);

    Vector3D normal1 = cross(points[1] - points[0], points[2] - points[0]).unitVector();
    Vector3D normal2 = cross(points[2] - points[0], points[3] - points[0]).unitVector();

    return (normal1 - normal2).mag() <= TOLER;
}

double Geometry::computeHexahedronVolume(Point3D* points)
{
    // This uses the method of calculating hexahedron volumes depicted in
    // Computational Fluid Mechanics and Heat Transfer, Third Edition, by
    // Richard H. Pletcher, John C. Tannehill and Gale Anderson

    Vector3D r16 = points[6] - points[1];
    Vector3D r25 = points[5] - points[2];
    Vector3D s1265 = 0.5*cross(r16, r25);

    Vector3D r57 = points[7] - points[5];
    Vector3D r64 = points[4] - points[6];
    Vector3D s5674 = 0.5*cross(r57, r64);

    Vector3D r27 = points[7] - points[2];
    Vector3D r36 = points[6] - points[3];
    Vector3D s2376 = 0.5*cross(r27, r36);

    return fabs(dot(s1265 + s5674 + s2376, points[6] - points[0])/3.);
}

Point3D Geometry::computeHexahedronCentroid(Point3D *points)
{
    // The strategy is to break the hexahedron into tetrahedrons. The centroid is then
    // defined as the volume weighted average of the tetrahedron centroids

    Point3D tet[5][4], tetCentroids[5], centroid(0., 0., 0.);
    double tetVolumes[5], sumVolumes(0.);

    // Break the hexahedron into 5 tetrahedrons

    tet[0][0] = points[4];
    tet[0][1] = points[3];
    tet[0][2] = points[6];
    tet[0][3] = points[7];

    tet[1][0] = points[4];
    tet[1][1] = points[1];
    tet[1][2] = points[6];
    tet[1][3] = points[5];

    tet[2][0] = points[4];
    tet[2][1] = points[1];
    tet[2][2] = points[6];
    tet[2][3] = points[3];

    tet[3][0] = points[4];
    tet[3][1] = points[0];
    tet[3][2] = points[1];
    tet[3][3] = points[3];

    tet[4][0] = points[3];
    tet[4][1] = points[1];
    tet[4][2] = points[2];
    tet[4][3] = points[6];

    for(int i = 0; i < 5; ++i)
    {
        tetCentroids[i] = computeTetrahedronCentroid(tet[i]);
        tetVolumes[i] = computeTetrahedronVolume(tet[i]);
        sumVolumes += tetVolumes[i];
    }

    for(int i = 0; i < 5; ++i)
    {
        centroid += tetVolumes[i]*tetCentroids[i];
    }

    return centroid/sumVolumes;
}

bool Geometry::checkHexahedronSurfacesIsPlanar(Point3D *points)
{
    Point3D tmpVertices[4];

    // Check east face

    tmpVertices[0] = points[1];
    tmpVertices[1] = points[2];
    tmpVertices[2] = points[6];
    tmpVertices[3] = points[5];

    if(!checkQuadrilateralIsPlanar(tmpVertices))
        return false;

    // Check west face

    tmpVertices[0] = points[0];
    tmpVertices[1] = points[4];
    tmpVertices[2] = points[7];
    tmpVertices[3] = points[3];

    if(!checkQuadrilateralIsPlanar(tmpVertices))
        return false;

    // Check north face

    tmpVertices[0] = points[2];
    tmpVertices[1] = points[3];
    tmpVertices[2] = points[7];
    tmpVertices[3] = points[6];

    if(!checkQuadrilateralIsPlanar(tmpVertices))
        return false;

    // Check south face

    tmpVertices[0] = points[0];
    tmpVertices[1] = points[1];
    tmpVertices[2] = points[5];
    tmpVertices[3] = points[4];

    if(!checkQuadrilateralIsPlanar(tmpVertices))
        return false;

    // Check top face

    tmpVertices[0] = points[4];
    tmpVertices[1] = points[5];
    tmpVertices[2] = points[6];
    tmpVertices[3] = points[7];

    if(!checkQuadrilateralIsPlanar(tmpVertices))
        return false;

    // Check bottom face

    tmpVertices[0] = points[0];
    tmpVertices[1] = points[3];
    tmpVertices[2] = points[2];
    tmpVertices[3] = points[1];

    if(!checkQuadrilateralIsPlanar(tmpVertices))
        return false;

    return true;
}

bool Geometry::isInsideHexahedron(const Point3D &point, const Point3D *points)
{
    Point3D tmpVertices[4], faceXc;
    Vector3D faceNorm;

    // Check east face
    tmpVertices[0] = points[1];
    tmpVertices[1] = points[2];
    tmpVertices[2] = points[6];
    tmpVertices[3] = points[5];

    faceXc = computeQuadrilateralCentroid(tmpVertices);
    faceNorm = computeQuadrilateralNormal(tmpVertices);

    if(dot(faceXc - point, faceNorm) < 0.)
        return false;

    // Check west face
    tmpVertices[0] = points[0];
    tmpVertices[1] = points[4];
    tmpVertices[2] = points[7];
    tmpVertices[3] = points[3];

    faceXc = computeQuadrilateralCentroid(tmpVertices);
    faceNorm = computeQuadrilateralNormal(tmpVertices);

    if(dot(faceXc - point, faceNorm) < 0.)
        return false;

    // Check north face
    tmpVertices[0] = points[2];
    tmpVertices[1] = points[3];
    tmpVertices[2] = points[7];
    tmpVertices[3] = points[6];

    faceXc = computeQuadrilateralCentroid(tmpVertices);
    faceNorm = computeQuadrilateralNormal(tmpVertices);

    if(dot(faceXc - point, faceNorm) < 0.)
        return false;

    // Check south face
    tmpVertices[0] = points[0];
    tmpVertices[1] = points[1];
    tmpVertices[2] = points[5];
    tmpVertices[3] = points[4];

    faceXc = computeQuadrilateralCentroid(tmpVertices);
    faceNorm = computeQuadrilateralNormal(tmpVertices);

    if(dot(faceXc - point, faceNorm) < 0.)
        return false;

    // Check top face
    tmpVertices[0] = points[4];
    tmpVertices[1] = points[5];
    tmpVertices[2] = points[6];
    tmpVertices[3] = points[7];

    faceXc = computeQuadrilateralCentroid(tmpVertices);
    faceNorm = computeQuadrilateralNormal(tmpVertices);

    if(dot(faceXc - point, faceNorm) < 0.)
        return false;

    // Check bottom face
    tmpVertices[0] = points[0];
    tmpVertices[1] = points[3];
    tmpVertices[2] = points[2];
    tmpVertices[3] = points[1];

    faceXc = computeQuadrilateralCentroid(tmpVertices);
    faceNorm = computeQuadrilateralNormal(tmpVertices);

    if(dot(faceXc - point, faceNorm) < 0.)
        return false;

    return true;
}

double Geometry::computeTetrahedronVolume(Point3D *points)
{
    Vector3D a, b, c;

    a = points[1] - points[0];
    b = points[2] - points[0];
    c = points[3] - points[0];

    return fabs(dot(a, cross(b, c)))/6.;
}

Point3D Geometry::computeTetrahedronCentroid(Point3D *points)
{
    Point3D sum(0., 0., 0.);

    for(int i = 0; i < 4; ++i)
        sum += points[i];

    return 0.25*sum;
}
