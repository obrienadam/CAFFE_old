#ifndef INPUT_H
#define INPUT_H

#include <string>
#include <fstream>
#include <map>

class Input
{
private:

    std::string filename_;
    std::ifstream fin_;

    void processBuffer(std::string& buffer);

public:

    std::map<std::string, int> inputInts;
    std::map<std::string, double> inputDoubles;
    std::map<std::string, std::string> inputStrings;

    Input();
    Input(std::string filename);
    ~Input();

    void openInputFile(std::string filename);

    //- Show everything in the input for debugging purposes

    void print();
};

#endif
