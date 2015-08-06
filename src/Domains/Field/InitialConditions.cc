#include "InitialConditions.h"

#include <boost/property_tree/info_parser.hpp>

InitialConditions::InitialConditions()
{
    using namespace boost::property_tree;

    read_info("case/initialConditions.info", initialConditionsParameters_);
}

template <>
void InitialConditions::setInitialConditions(Field<double> &field)
{
    std::string icType = initialConditionsParameters_.get<std::string>(field.name + ".type");

    if(icType == "uniform")
    {
        createUniform(initialConditionsParameters_.get<double>(field.name + ".value"), field);
    }
    else if(icType == "sphere")
    {
        createSphere(initialConditionsParameters_.get<double>(field.name + ".radius"),
                     std::stov(initialConditionsParameters_.get<std::string>(field.name + ".center")),
                     initialConditionsParameters_.get<double>(field.name + ".value"),
                     field);
    }
    else if(icType == "box")
    {
        createBox(initialConditionsParameters_.get<double>(field.name + ".xLength"),
                  initialConditionsParameters_.get<double>(field.name + ".yLength"),
                  initialConditionsParameters_.get<double>(field.name + ".yLength"),
                  std::stov(initialConditionsParameters_.get<std::string>(field.name + ".center")),
                  initialConditionsParameters_.get<double>(field.name + ".value"),
                  field);
    }
    else
    {
        Output::raiseException("InitialConditions", "setInitialConditions", "Unrecognized initial condition type \"" + icType + "\".");
    }
}

template<>
void InitialConditions::setInitialConditions(Field<Vector3D> &field)
{
    using namespace std;

    string icType = initialConditionsParameters_.get<string>(field.name + ".type");

    if(icType == "uniform")
    {
        createUniform(stov(initialConditionsParameters_.get<string>(field.name + ".value")), field);
    }
    else if(icType == "sphere")
    {
        createSphere(initialConditionsParameters_.get<double>(field.name + ".radius"),
                     stov(initialConditionsParameters_.get<string>(field.name + ".center")),
                     stov(initialConditionsParameters_.get<string>(field.name + ".value")),
                     field);
    }
    else if(icType == "box")
    {
        createBox(initialConditionsParameters_.get<double>(field.name + ".xLength"),
                  initialConditionsParameters_.get<double>(field.name + ".yLength"),
                  initialConditionsParameters_.get<double>(field.name + ".yLength"),
                  stov(initialConditionsParameters_.get<string>(field.name + ".center")),
                  stov(initialConditionsParameters_.get<string>(field.name + ".value")),
                  field);
    }
    else
    {
        Output::raiseException("InitialConditions", "setInitialConditions", "Unrecognized initial condition type \"" + icType + "\".");
    }
}
