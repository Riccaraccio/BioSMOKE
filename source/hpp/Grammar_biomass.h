#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OpenSMOKE
{
    class Grammar_Biomass : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
    {
    protected:
        
        virtual void DefineRules()
        {
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsFolder", 
                                                                OpenSMOKE::SINGLE_PATH, 
                                                                "Name of the folder containing the kinetic scheme (XML Version)", 
                                                                true,
                                                                "@KineticsPreProcessor",
                                                                "none",
                                                                "none") );	
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@KineticsPreProcessor", 
                                                                OpenSMOKE::SINGLE_DICTIONARY, 
                                                                "Name of the dictionary containing the list of kinetic files to be interpreted", 
                                                                true,
                                                                "@KineticsFolder",
                                                                "none",
                                                                "none") );
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Analysis", 
                                                                OpenSMOKE::SINGLE_STRING, 
                                                                "Type of run: only TG | Complete analysis", 
                                                                true) );
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TGA", 
                                                                OpenSMOKE::SINGLE_DICTIONARY, 
                                                                "Name of the dictionary containing the command for the TGA simulation", 
                                                                true) );
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Total_Analysis", 
                                                                OpenSMOKE::SINGLE_DICTIONARY, 
                                                                "Name of the dictionary containing the command for the complete simulation", 
                                                                true) );
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InletGasStatus", 
                                                                OpenSMOKE::SINGLE_DICTIONARY, 
                                                                "Name of the dictionary defining the inlet gas composition, temperature and pressure", 
                                                                false) );
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InletSolidStatus", 
                                                                OpenSMOKE::SINGLE_DICTIONARY, 
                                                                "Name of the dictionary defining the inlet gas composition, temperature and pressure", 
                                                                false) );
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OdeParameters", 
                                                                OpenSMOKE::SINGLE_DICTIONARY, 
                                                                "Dictionary containing the numerical parameters for solving the stiff ODE system", 
                                                                false) );	
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TemperatureProfile",
																OpenSMOKE::SINGLE_DICTIONARY,
																"Dictionary containing the numerical parameters for porous medium",
																false));
        }
    };
}
        