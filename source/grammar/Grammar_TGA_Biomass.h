#ifndef OpenSMOKE_Grammar_TGA_analysis_H
#define OpenSMOKE_Grammar_TGA_analysis_H

#include <string>
#include <dictionary/OpenSMOKE_Dictionary.h>
#include <dictionary/OpenSMOKE_DictionaryGrammar.h>

namespace BioSMOKE
{
class Grammar_TGA_analysis : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
{
  protected:
    virtual void DefineRules()
    {

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@HeatingRate", OpenSMOKE::SINGLE_MEASURE,
                                                          "Heat rates of particle (i.e. 1 °C/min)", false));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TotalSimulationTime", OpenSMOKE::SINGLE_MEASURE,
                                                          "Total time to simulate(i.e. 100 s)", false));

        // AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
        //     "@OutputSpecies", OpenSMOKE::VECTOR_STRING, "List of species which will be written on ASCII file",
        //     false));
    }
};

inline void Get_TGAanalysisFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary &dictionary, std::string &analysis,
                                          double &Heat_Rates, double &final_time,
                                          std::vector<std::string> &output_species)
{
    Grammar_TGA_analysis grammar_TGA_analysis;
    dictionary.SetGrammar(grammar_TGA_analysis);

    analysis = "TGA";

    // Heating rate
    {
        if (dictionary.CheckOption("@HeatingRate") == true)
        {
            double value;
            std::string units;
            dictionary.ReadMeasure("@HeatingRate", value, units);

            if (units == "K/s")
                Heat_Rates = value;
            else if (units == "K/min")
                Heat_Rates = value / 60.;
            else if (units == "K/h")
                Heat_Rates = value / 3600.;
            else
                OpenSMOKE::FatalErrorMessage("Unknown heating rate unit of measurement");
        }
    }

    // Read residence time
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
                OpenSMOKE::FatalErrorMessage("Unknown time unit of measurement");
        }
    }
}

} // namespace BioSMOKE

#endif /* OpenSMOKE_Grammar_TGA_analysis_H */
