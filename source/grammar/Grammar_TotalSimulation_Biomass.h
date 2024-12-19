#ifndef OpenSMOKE_Grammar_TotalSimulation_analysis_H
#define OpenSMOKE_Grammar_TotalSimulation_analysis_H

#include <string>
#include <dictionary/OpenSMOKE_Dictionary.h>
#include <dictionary/OpenSMOKE_DictionaryGrammar.h>
#include <dictionary/OpenSMOKE_DictionaryKeyWord.h>

namespace BioSMOKE
{
class Grammar_TotalSimulation_analysis : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
{
  protected:
    virtual void DefineRules()
    {

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Type", OpenSMOKE::SINGLE_STRING,
                                                          "Type of run: Isothermal |  Non-Isothermal", true));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ConstantVolume", OpenSMOKE::SINGLE_BOOL, 
                                                          "Wethere to run a simualtion at constant volume or constant porosity", true));

        AddKeyWord(
            OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Porosity", OpenSMOKE::SINGLE_DOUBLE, "Void fraction  (-)", true));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Nlayer", OpenSMOKE::SINGLE_INT,
                                                          "Number of discretization intervals (-)", true));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InitialRadius", OpenSMOKE::SINGLE_MEASURE,
                                                          "Initial sphere radius (es. 0.1m)", true));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TotalSimulationTime", OpenSMOKE::SINGLE_MEASURE,
                                                          "Total time to simulate(i.e. 100 s)", true));

        AddKeyWord(
            OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Darcy_number", OpenSMOKE::SINGLE_DOUBLE, "Darcy number", true));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
            "@hext", OpenSMOKE::SINGLE_MEASURE, "Convective heat exchange coefficient (i.e. 20 W/m2/K)", true));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Lambda", OpenSMOKE::SINGLE_MEASURE,
                                                          "Thermal conductivity for the solid", true));
    }
};

inline void Get_TotalSimulation_analysisFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary &dictionary,
                                                       std::string &analysis, bool &energy_balance, bool &volume_loss,
                                                       double &final_time, double &porosity, int &number_of_layers,
                                                       double &initial_radius, double &Da_number, double &hext,
                                                       double &lambda_solid, std::vector<std::string> &output_species)
{
    Grammar_TotalSimulation_analysis grammar_TotalSimulation_analysis;
    dictionary.SetGrammar(grammar_TotalSimulation_analysis);

    analysis = "Total_Analysis";

    // Type of simulation
    {
        std::string value;
        if (dictionary.CheckOption("@Type") == true)
        {
            dictionary.ReadString("@Type", value);
            if (value == "Isothermal")
                energy_balance = false;
            else if (value == "Non_Isothermal")
                energy_balance = true;
            else
            {
                OpenSMOKE::FatalErrorMessage("Unknown simulation type: " + value +
                                             "  Conditions:Isothermal or Non-Isothermal");
            }
        }
    }

    // Initial radius
    {
        if (dictionary.CheckOption("@InitialRadius") == true)
        {
            double value;
            std::string units;
            dictionary.ReadMeasure("@InitialRadius", value, units);
            if (units == "m")
                initial_radius = value;
            else if (units == "cm")
                initial_radius = value / 100.;
            else if (units == "mm")
                initial_radius = value / 1000.;
            else
                OpenSMOKE::FatalErrorMessage("Unknown initial radius units of measurment");
        }
    }

    // Volume
    {
        std::string value;
        if (dictionary.CheckOption("@ConstantVolume") == true)
        {
            dictionary.ReadString("@ConstantVolume", value);
            if (value == "true")
                volume_loss = false;
            else if (value == "false")
                volume_loss = true;
            else
            {
                OpenSMOKE::FatalErrorMessage("Unknown value for ConstantVolume: " + value + "  Conditions: true or false");
            }
        }
    }

    // Total simulation time
    {
        if (dictionary.CheckOption("@TotalSimulationTime") == true)
        {
            double value;
            std::string units;
            dictionary.ReadMeasure("@TotalSimulationTime", value, units);
            if (units == "s")
                final_time = value;
            else if (units == "ms")
                final_time = value / 1000.;
            else if (units == "min")
                final_time = value * 60.;
            else if (units == "h")
                final_time = value * 3600.;
            else
                OpenSMOKE::FatalErrorMessage("Unknown total simulation time units");
        }
    }

    // porosity
    {
        if (dictionary.CheckOption("@Porosity") == true)
            dictionary.ReadDouble("@Porosity", porosity);
    }

    // Number of layers
    {
        if (dictionary.CheckOption("@Nlayer") == true)
            dictionary.ReadInt("@Nlayer", number_of_layers);
    }

    // Darcy number
    {
        if (dictionary.CheckOption("@Darcy_number") == true)
            dictionary.ReadDouble("@Darcy_number", Da_number);
    }

    // Heat exchange coefficient
    {
        if (dictionary.CheckOption("@hext") == true)
        {
            double value;
            std::string units;
            dictionary.ReadMeasure("@hext", value, units);

            if (units == "W/m2/K")
                hext = value;
            else
                OpenSMOKE::FatalErrorMessage("Unknown heat exchange coefficent units");
        }
    }

    // Solid thermal conductivity
    {
        if (dictionary.CheckOption("@lambda") == true)
        {
            double value;
            std::string units;
            dictionary.ReadMeasure("@lambda", value, units);

            if (units == "W/m/K")
                lambda_solid = value;
            else
                OpenSMOKE::FatalErrorMessage("Unknown thermal conductivity units");
        }
    }
}
} // namespace BioSMOKE

#endif /* OpenSMOKE_Grammar_TotalSimulation_analysis_H */
