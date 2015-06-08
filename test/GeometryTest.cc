#include <iostream>
#include <math.h>

#include "../src/Math/Point3D.h"
#include "../src/Math/Geometry.h"

int main()
{
  using namespace std;

  Point3D cubeVertices[8];
  double dx(-0.2), dy(-0.3), dz(-0.423);

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
  
  return 0;
}
