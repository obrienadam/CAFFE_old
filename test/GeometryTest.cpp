#include <iostream>
#include <math.h>

#include "../src/Math/Point3D.h"
#include "../src/Math/Geometry.h"
#include "../src/Math/Sphere.h"

int main()
{
  using namespace std;

  Point3D cubeVertices[8];
  double dx(1), dy(1), dz(1);

  cubeVertices[0] = Point3D(0., 0., 0.);
  cubeVertices[1] = Point3D(dx, 0., 0.);
  cubeVertices[2] = Point3D(dx, dy, 0.);
  cubeVertices[3] = Point3D(0., dy, 0.);
  cubeVertices[4] = Point3D(0., 0., dz);
  cubeVertices[5] = Point3D(dx, 0., dz);
  cubeVertices[6] = Point3D(dx, dy, dz);
  cubeVertices[7] = Point3D(0., dy, dz);

  cout << "Cube Volume: " << fabs(dx*dy*dz) << endl
       << "Computed: " << Geometry::computeHexahedronVolume(cubeVertices) << endl;

  cout << "XY-plane Area: " << fabs(dx*dy) << endl
       << "Computed: " << Geometry::computeQuadrilateralArea(cubeVertices) << endl;

  cout << "Testing the sphere nearest intersect method...\n";

  Sphere sphere(1., Point3D(0., 0., 0.));
  Point3D testPt(1.2, 1.1, 1.1);
  
  cout << "Sphere center     : " << sphere.center << endl
       << "Sphere radius     : " << sphere.radius << endl
       << "Test point        : " << testPt << endl
       << "Nearest intersect : " << sphere.nearestIntersect(testPt) << endl;

  cout << "Testing the sphere line segment intersect method...\n";

  Point3D pt1(0.5, 0.5, 0.5), pt2(2.3, 1.2, 1.5);

  cout << "Sphere center     : " << sphere.center << endl
       << "Sphere radius     : " << sphere.radius << endl
       << "Test point 1      : " << pt1 << endl
       << "Test point 2      : " << pt2 << endl
       << "Intersect point 1 : " << sphere.lineIntersect(pt1, pt2).first << endl
       << "Intersect point 2 : " << sphere.lineIntersect(pt1, pt2).second << endl;

  testPt = Point3D(0.8, 0.5, 0.5);

  cout << "Testing if a point " << testPt << " is within a hexahedron, for a cube defined with parameters " << dx << " " << dy << " " << dz << endl;

  cout << Geometry::isInsideHexahedron(testPt, cubeVertices) << endl;

  return 0;
}
