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

template<>
void InitialConditions::writeRestart(const Field<double> &field)
{
    int i, j, k;
    std::ofstream fout((field.name + "Restart.rst").c_str());
    const HexaFvmMesh &mesh = field.getMesh();

    for(k = 0; k < mesh.nCellsK(); ++k)
    {
        for(j = 0; j < mesh.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh.nCellsI(); ++i)
            {
                fout << field(i, j, k) << " ";
            }

            fout << std::endl;
        }
        fout << std::endl;
    }
}

template<>
void InitialConditions::writeRestart(const Field<Vector3D> &field)
{
    int i, j, k;
    std::ofstream fout((field.name + "Restart.rst").c_str());
    const HexaFvmMesh &mesh = field.getMesh();

    for(k = 0; k < mesh.nCellsK(); ++k)
    {
        for(j = 0; j < mesh.nCellsJ(); ++j)
        {
            for(i = 0; i < mesh.nCellsI(); ++i)
            {
                fout << field(i, j, k).x << " " << field(i, j, k).y << " " << field(i, j, k).z << " ";
            }

            fout << std::endl;
        }
        fout << std::endl;
    }
}
