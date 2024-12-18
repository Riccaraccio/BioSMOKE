#ifndef OpenSMOKE_Grammar_TGA_analysis_H
#define	OpenSMOKE_Grammar_TGA_analysis_H

#include <string>
#include "boost/filesystem.hpp"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

namespace OpenSMOKE
{
    class Grammar_TGA_analysis : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
    {
    protected:
        virtual void DefineRules()
	{
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Heat_Rates", 
                                                                OpenSMOKE::SINGLE_MEASURE, 
                                                                "Heat rates of particle (i.e. 1 °C/min)", 
                                                                false) );
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SimulationTime", 
                                                                OpenSMOKE::SINGLE_MEASURE, 
                                                                "Total time to simulate(i.e. 100 s)", 
                                                                false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OutputSpecies",
															   OpenSMOKE::VECTOR_STRING,
															  "List of species which will be written on ASCII file",
															  false));
        }
    };
    
    void Get_TGAanalysisFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, std::string& analysis, 
                                       double& Heat_Rates, double& final_time, std::vector<std::string>& output_species )
    {
        Grammar_TGA_analysis grammar_TGA_analysis;
        dictionary.SetGrammar(grammar_TGA_analysis);

        
        analysis = "TGA";
           
        // Heat rates
        {
            if (dictionary.CheckOption("@Heat_Rates") == true)
            {
                double value;
                std::string units;
                dictionary.ReadMeasure("@Heat_Rates", value, units);

                if (units == "K/s")			Heat_Rates = value;
                else if (units == "K/min")		Heat_Rates = value/60.;
                else if (units ==  "K/h")               Heat_Rates = value/3600.;
                else OpenSMOKE::FatalErrorMessage("Unknown heat rates units");
            }
        }
        
        // Read residence time
        {
            if (dictionary.CheckOption("@SimulationTime") == true)
            {
                double value;
                std::string units;
                dictionary.ReadMeasure("@SimulationTime", value, units);
                if (units == "s")	  final_time = value;
                else if (units == "ms")   final_time = value/1000.;
                else if (units == "min")  final_time = value*60.;
                else if (units == "h")    final_time = value*3600.;
                else OpenSMOKE::FatalErrorMessage("Unknown time units");

            }

			if (dictionary.CheckOption("@OutputSpecies") == true)
			{
				std::vector<std::string> names;
				dictionary.ReadOption("@OutputSpecies",names);

				for (unsigned int j = 0; j < names.size(); j++)
				{
					output_species.push_back(names[j]);
				}
			}
            else{
                output_species.push_back("all");
            }
				
        }
    }

}

#endif	/* OpenSMOKE_Grammar_TGA_analysis_H */

