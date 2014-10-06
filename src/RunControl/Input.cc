#include "Input.h"

Input::Input()
{

}

Input::~Input()
{
  if(inFile_.is_open())
    inFile_.close();
}
