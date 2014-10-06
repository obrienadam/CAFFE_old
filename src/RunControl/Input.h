#ifndef INPUT_H
#define INPUT_H

#include <string>
#include <fstream>

class Input
{

 protected:

  std::ifstream inFile_;

 public:

  Input();
  ~Input();

  void openInputFile(std::string filename);

};

#endif
