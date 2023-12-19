#ifndef OpenSMOKE_Grammar_TotalSimulation_analysis_H
#define	OpenSMOKE_Grammar_TotalSimulation_analysis_H

#include <string>
#include "boost/filesystem.hpp"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

namespace OpenSMOKE
{
    class Grammar_TotalSimulation_analysis : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
    {
    protected:
        virtual void DefineRules()
	{
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Type", 
                                                                OpenSMOKE::SINGLE_STRING, 
                                                                "Type of run: Isothermal |  Adiabatic", 
                                                                true) );
            
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Volume", 
                                                                OpenSMOKE::SINGLE_STRING, 
                                                                "Type of run: ConstantVolume |  VariableVolume", 
                                                                true) );
            
                                                                	
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Epsi", 
                                                                OpenSMOKE::SINGLE_DOUBLE, 
                                                                "Void fraction  (-)", 
                                                                true ) );
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Tau", 
                                                                OpenSMOKE::SINGLE_DOUBLE, 
                                                                "Tortuosity factor  (-)", 
                                                                false ) );
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Nlayer", 
                                                                OpenSMOKE::SINGLE_INT, 
                                                                "Number of interval of sphere (-)", 
                                                                true ) );
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Radius", 
                                                                OpenSMOKE::SINGLE_MEASURE, 
                                                                "Sphere raggio (es. 0.1m)", 
                                                                true ) );
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SimulationTime", 
                                                                OpenSMOKE::SINGLE_MEASURE, 
                                                                "Total time to simulate(i.e. 100 s)", 
                                                                true));
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Darcy_number", 
                                                                OpenSMOKE::SINGLE_DOUBLE, 
                                                                "Darcy number", 
                                                                true) );

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OutputSpecies",
				OpenSMOKE::VECTOR_STRING,
				"List of species which will be written on ASCII file",
				false));
        }
    };
    
    void Get_TotalSimulation_analysisFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, std::string& analysis, bool& energyBalance,   
                                                    bool& volumeLoss, double& final_time, double& epsi, double& tau, int& ntr, double& raggioTot,double& Da, std::vector<std::string>& output_species)
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
                if (value == "Isothermal")	energyBalance = false;
                else if (value == "Non_Isothermal")	energyBalance = true;
                else 
                {
                    OpenSMOKE::FatalErrorMessage("Unknown simulation type: " + value + "  Conditions:Isothermal or Non_Isothermal");

                }
            }
        }
        
        //RaggioTot
        {
            if (dictionary.CheckOption("@Radius") == true)
            {
                double value;
                std::string units;
                dictionary.ReadMeasure("@Radius", value, units);
                if (units == "m")	      raggioTot = value;
                else if (units == "cm")   raggioTot = value/100.;
                else if (units == "mm")   raggioTot = value/1000.;
                else OpenSMOKE::FatalErrorMessage("Unknown time units");
            }
        }
        
        
        // Volume Condition
        {
            std::string value;
            if (dictionary.CheckOption("@Volume") == true)
            {
                dictionary.ReadString("@Volume", value);
                if (value == "ConstantVolume")	volumeLoss = false;
                else if (value == "VariableVolume")	volumeLoss = true;
                else 
                {
                    OpenSMOKE::FatalErrorMessage("Unknown simulation type: " + value + "  Conditions:ConstantVolume or VariableVolume");

                }
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
        }

        // epsi
        {
            if (dictionary.CheckOption("@Epsi") == true)
                dictionary.ReadDouble("@Epsi", epsi);       
        }

        // tau
        {
            if (dictionary.CheckOption("@Tau") == true)
                dictionary.ReadDouble("@Tau", tau);       
        }

        // intervall
        {
            if (dictionary.CheckOption("@Nlayer") == true)
                dictionary.ReadInt("@Nlayer", ntr);       
        }
        
        // Darcy
        {
            if (dictionary.CheckOption("@Darcy_number") == true)
                dictionary.ReadDouble("@Darcy_number", Da);       
        }

		if (dictionary.CheckOption("@OutputSpecies") == true)
		{
			std::vector<std::string> names;
			dictionary.ReadOption("@OutputSpecies", names);

			for (unsigned int j = 0; j < names.size(); j++)
			{
				output_species.push_back(names[j]);
			}
		}
    }
}

#endif	/* OpenSMOKE_Grammar_TotalSimulation_analysis_H */

