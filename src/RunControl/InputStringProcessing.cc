#include <boost/algorithm/string.hpp>

#include "InputStringProcessing.h"

std::string InputStringProcessing::processBuffer(std::string &buffer, bool removeAllWhitespace)
{
    using namespace boost::algorithm;

    // Trim the leading/lagging whitespace

    trim(buffer);

    // Remove all whitespace from buffer if desired (default)

    if(removeAllWhitespace)
        erase_all(buffer, " ");

    // Check if it is a comment line, if so discard the input

    if(buffer[0] == '#')
    {
        buffer.clear();
    }

    // Remove any comments on the line

    buffer = buffer.substr(0, buffer.find("#"));
    return buffer;
}

double InputStringProcessing::getNextElement(std::string &buffer)
{
    using namespace boost::algorithm;

    double element;

    // Ensure the buffer is trimmed

    trim_left_if(buffer, is_any_of("( "));

    // Extract a double element from the string

    element = stod(buffer.substr(0, buffer.find_first_of(" )")));

    // Remove the extracted element from the string

    buffer = buffer.substr(buffer.find_first_of(" )"), buffer.back());
    return element;
}

std::vector<std::string> InputStringProcessing::partition(std::string buffer, std::string delimiter)
{
    using namespace std;
    using namespace boost::algorithm;

    vector<string> partitionedString;

    split(partitionedString, buffer, is_any_of(delimiter));

    return partitionedString;
}
