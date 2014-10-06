#include <iostream>

#include "RunControl.h"
#include "HexaFdmMesh.h"
#include "AdvectionDiffusionField.h"

using namespace std;

int main(int argc, const char* argv[])
{
  RunControl runControl(argc, argv);
  HexaFdmMesh mesh(30, 30, 30);
  AdvectionDiffusionField phiField(30, 30, 30);

  runControl.displayStartMessage();

  while(runControl.continueRun())
    {

    }

  return 0;
}
