#include <iostream>

#include "RunControl.h"
#include "HexaFdmMesh.h"
#include "AdvectionDiffusionField.h"

int main(int argc, const char* argv[])
{
  using namespace std;

  RunControl runControl(argc, argv);
  HexaFdmMesh mesh(30, 30, 30);

  runControl.displayStartMessage();

  while(runControl.continueRun())
    {

    }

  return 0;
}
