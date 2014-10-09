#include <iostream>
#include <cstdlib>

#include "../src/DataStructures/SmartPointer3D.h"

typedef SmartPointer3D<double> double3D;

int main()
{

  using namespace std;

  double3D ran3D(17, 25, 19);
  
  double3D::iterator begin = ran3D.begin();
  double3D::iterator end = ran3D.end();
  double3D::iterator itr;
  
  int i = 1;

  for(itr = begin; itr != end; ++itr)
    {
      *itr = i;
      cout << *itr << endl;
      ++i;
    }
  
return 0;
}
