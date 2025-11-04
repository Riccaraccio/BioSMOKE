#pragma once

namespace BioSMOKE
{
class Grammar_BioSMOKE : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
{
  protected:
    virtual void DefineRules()
    {
        // clang-format off
        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsFolder",
                                                          OpenSMOKE::SINGLE_PATH,
                                                          "Name of the folder containing the kinetic scheme (XML Version)",
                                                          true,
                                                          "@KineticsPreProcessor",
                                                          "none",
                                                          "none"));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsPreProcessor",
                                                          OpenSMOKE::SINGLE_DICTIONARY,
                                                          "Name of the dictionary containing the list of kinetic files to be interpreted",
                                                          true,
                                                          "@KineticsFolder",
                                                          "none",
                                                          "none"));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Type",
                                                          OpenSMOKE::SINGLE_STRING,
                                                          "Type of run: only TG | Complete analysis",
                                                          true));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TGA",
                                                          OpenSMOKE::SINGLE_DICTIONARY,
                                                          "Name of the dictionary containing the command for the TGA simulation",
                                                          false));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Particle",
                                                          OpenSMOKE::SINGLE_DICTIONARY,
                                                          "Name of the dictionary containing the command for the complete simulation",
                                                          false));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InletGasStatus",
                                                          OpenSMOKE::SINGLE_DICTIONARY,
                                                          "Name of the dictionary defining the inlet gas composition, temperature and pressure",
                                                          true));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InletSolidStatus",
                                                          OpenSMOKE::SINGLE_DICTIONARY,
                                                          "Name of the dictionary defining the inlet gas composition, temperature and pressure",
                                                          true));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OdeParameters",
                                                          OpenSMOKE::SINGLE_DICTIONARY,
                                                          "Dictionary containing the numerical parameters for solving the stiff ODE system",
                                                          false));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TemperatureProfile",
                                                          OpenSMOKE::SINGLE_DICTIONARY,
                                                          "Dictionary containing the numerical parameters for porous medium",
                                                          false));
        
        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Options",
                                                          OpenSMOKE::SINGLE_DICTIONARY,
                                                          "Dictionary containing additional options for BioSMOKE",
                                                          false));

        // TODO
        // AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SensitivityAnalysis", 
				// 												OpenSMOKE::SINGLE_DICTIONARY, 
				// 												"Dictionary containing additional options for solving the sensitivity analysis", 
				// 												false) );

        // clang-format on
    }
};
} // namespace BioSMOKE
