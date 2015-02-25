#include <sstream>

#include "Version.h"
#include "Output.h"

void Version::display()
{
    using namespace std;

    ostringstream versionInfo;

    versionInfo << "Computational Algorithms for Fluid Equations (CAFFE) 0.0.0\n"
                << "This is free software; see the source for copying conditions. There is NO\n"
                << "warranty; not even for merchantability or fitness for a particular purpose.";


    Output::printToScreen(versionInfo.str());
}
