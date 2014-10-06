#include <iostream>

#include "HexaFdmMesh.h"
#include "RunControl.h"

using namespace std;

int main(int argc, const char* argv[])
{
  RunControl runControl(argc, argv);
  HexaFdmMesh mesh(30, 30, 30);

  runControl.displayStartMessage();

  while(runControl.continueRun())
    {

    }

  return 0;
}
