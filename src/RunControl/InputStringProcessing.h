#ifndef INPUTSTRINGPROCESSING_H
#define INPUTSTRINGPROCESSING_H

#include <string>
#include <vector>

class InputStringProcessing
{
private:

public:

    //- Processes an input line. Removes comments and optionally removes whitespace

    static std::string processBuffer(std::string& buffer, bool removeAllWhitespace = true);

    //- Read in the next element. Assume that the line represents doubles in the format (n1 n2 n3 ... N)

    static double getNextElement(std::string &buffer);

    //- Partition the line, using a specified delimiting character

    static std::vector<std::string> partition(std::string buffer, std::string delimiter);
};

#endif
