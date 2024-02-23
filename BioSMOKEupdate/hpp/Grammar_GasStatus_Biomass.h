#ifndef OpenSMOKE_Grammar_Gas_H
#define	OpenSMOKE_Grammar_Gas_H

#include <string>
#include "boost/filesystem.hpp"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

namespace OpenSMOKE
{
    class Grammar_Gas : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
    {
    protected:
        virtual void DefineRules()
	{
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Temperature", 
                                                                OpenSMOKE::SINGLE_MEASURE, 
                                                                "Temperature of the mixture (i.e. 500 K)", 
                                                                false) );	

            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Pressure", 
                                                                OpenSMOKE::SINGLE_MEASURE, 
                                                                "Pressure of the mixture (i.e. 1 atm)", 
                                                                false) );
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Density", 
                                                                OpenSMOKE::SINGLE_MEASURE, 
                                                                "Density of the mixture (i.e. 1 g/cm3)", 
                                                                false) );	
			
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MoleFractions", 
                                                                OpenSMOKE::VECTOR_STRING_DOUBLE, 
                                                                "Mole fractions of the mixture (i.e. CH4 0.60 H2 0.40)", 
                                                                true,
                                                                "@MassFractions @Moles @Masses" ,
                                                                "none",
                                                                "none") );

            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MassFractions", 
                                                                OpenSMOKE::VECTOR_STRING_DOUBLE, 
                                                                "Mass fractions of the mixture (i.e. CH4 0.60 H2 0.40)", 
                                                                true,
                                                                "@MoleFractions @Moles @Masses" ,
                                                                "none",
                                                                "none") );

            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Moles", 
                                                                OpenSMOKE::VECTOR_STRING_DOUBLE, 
                                                                "Moles (relative) of the mixture (i.e. CH4 2 H2 1)", 
                                                                true,
                                                                "@MoleFractions @MassFractions @Masses",
                                                                "none",
                                                                "none") );

            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Masses", 
                                                                OpenSMOKE::VECTOR_STRING_DOUBLE, 
                                                                "Masses (relative) of the mixture (i.e. CH4 2 H2 1)", 
                                                                true,
                                                                "@MoleFractions @MassFractions @Moles",
                                                                "none",
                                                                "none") );
 

        }
    };
    
    void GetGasStatusFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMapXML,
				    double& T, double& P_Pa, double& rho, OpenSMOKE::OpenSMOKEVectorDouble& omega, OpenSMOKE::OpenSMOKEVectorDouble& omegaSup)
    {
        Grammar_Gas grammar_gas_status;
        dictionary.SetGrammar(grammar_gas_status);

        unsigned int state_variables = 0;
        bool temperature_assigned = false;
        bool pressure_assigned = false;
        bool density_assigned = false;

        // Temperature
        {
            if (dictionary.CheckOption("@Temperature") == true)
            {
                double value;
                std::string units;
                dictionary.ReadMeasure("@Temperature", value, units);

                if (units == "K")			T = value;
                else if (units == "C")		T = value + 273.15;
                else OpenSMOKE::FatalErrorMessage("Unknown temperature units");

                state_variables++;
                temperature_assigned = true;
            }
        }

        // Pressure
        {
            if (dictionary.CheckOption("@Pressure") == true)
            {
                double value;
                std::string units;
                dictionary.ReadMeasure("@Pressure", value, units);

                if (units == "Pa")			P_Pa = value;
                else if (units == "bar")	P_Pa = value*1.e5;
                else if (units == "atm")	P_Pa = value*101325.;
                else OpenSMOKE::FatalErrorMessage("Unknown pressure units");

                state_variables++;
                pressure_assigned = true;
            }
        }

        // Density
        {
            if (dictionary.CheckOption("@Density") == true)
            {
                double value;
                std::string units;
                dictionary.ReadMeasure("@Density", value, units);

                if (units == "kg/m3")		rho = value;
                else if (units == "g/cm3")	rho = value*1.e3;
                else OpenSMOKE::FatalErrorMessage("Unknown density units");

                state_variables++;
                density_assigned = true;
            }
        }

        if (state_variables != 2)
                OpenSMOKE::FatalErrorMessage("The status of a gas mixture requires any 2 (and only 2) among: @Temperature, @Pressure and @Density");

        // Composition Internal
        {
            {
                std::vector<std::string> names;
                std::vector<double> values;

                if (dictionary.CheckOption("@MoleFractions") == true)
                {
                    dictionary.ReadOption("@MoleFractions", names, values);

                    const double sum =std::accumulate(values.begin(),values.end(),0.);
                    if (sum<(1.-1e-6) || sum>(1.+1e-6))
                            OpenSMOKE::FatalErrorMessage("The mole fractions must sum to 1.");

                    OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsMapXML.NumberOfSpecies());
                    for(unsigned int i=1;i<thermodynamicsMapXML.NumberOfSpecies();i++)
                           x[thermodynamicsMapXML.IndexOfSpecies(names[i-1])] = values[i-1]/sum;
                                                                    
                    ChangeDimensions(thermodynamicsMapXML.NumberOfSpecies(), &omega, true);
                    double MW;
                    thermodynamicsMapXML.MassFractions_From_MoleFractions(omega.GetHandle(),MW,x.GetHandle());
                    
                }
                else if (dictionary.CheckOption("@MassFractions") == true)
                {
                    dictionary.ReadOption("@MassFractions", names, values);
                    
                    const double sum =std::accumulate(values.begin(),values.end(),0.);
                    if (sum<(1.-1e-6) || sum>(1.+1e-6))
                            OpenSMOKE::FatalErrorMessage("The mass fractions must sum to 1.");
                   
                    ChangeDimensions(thermodynamicsMapXML.NumberOfSpecies(), &omega, true);
                    
                    for(unsigned int i=0;i<names.size();i++)
                        omega[thermodynamicsMapXML.IndexOfSpecies(names[i])] = values[i]/sum;
                    
                }
                else if (dictionary.CheckOption("@Moles") == true)
                {
                    dictionary.ReadOption("@Moles", names, values);
                    const double sum =std::accumulate(values.begin(),values.end(),0.);
                    OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsMapXML.NumberOfSpecies());
                    
                    for(unsigned int i=0;i<names.size();i++)
                        x[thermodynamicsMapXML.IndexOfSpecies(names[i])] = values[i]/sum;
                    
                    ChangeDimensions(thermodynamicsMapXML.NumberOfSpecies(), &omega, true);
                    double MW;
                    thermodynamicsMapXML.MassFractions_From_MoleFractions(omega.GetHandle(),MW,x.GetHandle());
                    
                }
                else  
                {
                    dictionary.ReadOption("@Masses", names, values);
                    const double sum =std::accumulate(values.begin(),values.end(),0.);
                    ChangeDimensions(thermodynamicsMapXML.NumberOfSpecies(), &omega, true);
                    for(unsigned int i=0;i<names.size();i++)
                            omega[thermodynamicsMapXML.IndexOfSpecies(names[i])] = values[i]/sum;
                    
                }
            }
        }
        
        // Composition Surface
        {
            {
                std::vector<std::string> names;
                std::vector<double> values;

                if (dictionary.CheckOption("@MoleFractionsSup") == true)
                {
                    dictionary.ReadOption("@MoleFractionsSup", names, values);

                    const double sum =std::accumulate(values.begin(),values.end(),0.);
                    if (sum<(1.-1e-6) || sum>(1.+1e-6))
                            OpenSMOKE::FatalErrorMessage("The mole fractions must sum to 1.");

                    OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsMapXML.NumberOfSpecies());
                    for(unsigned int i=0;i<names.size();i++)
                            x[thermodynamicsMapXML.IndexOfSpecies(names[i])] = values[i]/sum;
                    ChangeDimensions(thermodynamicsMapXML.NumberOfSpecies(), &omegaSup, true);
                    double MW;
                    thermodynamicsMapXML.MassFractions_From_MoleFractions(omegaSup.GetHandle(),MW,x.GetHandle());
                }
                else if (dictionary.CheckOption("@MassFractionsSup") == true)
                {
                    dictionary.ReadOption("@MassFractionsSup", names, values);

                    const double sum =std::accumulate(values.begin(),values.end(),0.);
                    if (sum<(1.-1e-6) || sum>(1.+1e-6))
                            OpenSMOKE::FatalErrorMessage("The mass fractions must sum to 1.");

                    ChangeDimensions(thermodynamicsMapXML.NumberOfSpecies(), &omegaSup, true);
                    for(unsigned int i=0;i<names.size();i++)
                            omegaSup[thermodynamicsMapXML.IndexOfSpecies(names[i])] = values[i]/sum;
                }
                else if (dictionary.CheckOption("@MolesSup") == true)
                {
                    dictionary.ReadOption("@MolesSup", names, values);
                    const double sum =std::accumulate(values.begin(),values.end(),0.);
                    OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsMapXML.NumberOfSpecies());
                    for(unsigned int i=0;i<names.size();i++)
                            x[thermodynamicsMapXML.IndexOfSpecies(names[i])] = values[i]/sum;
                    ChangeDimensions(thermodynamicsMapXML.NumberOfSpecies(), &omegaSup, true);
                    double MW;
                    thermodynamicsMapXML.MassFractions_From_MoleFractions(omegaSup.GetHandle(),MW,x.GetHandle());
                }
                else  
                {
                    dictionary.ReadOption("@MassesSup", names, values);
                    const double sum =std::accumulate(values.begin(),values.end(),0.);
                    ChangeDimensions(thermodynamicsMapXML.NumberOfSpecies(), &omegaSup, true);
                    for(unsigned int i=0;i<names.size();i++)
                            omegaSup[thermodynamicsMapXML.IndexOfSpecies(names[i])] = values[i]/sum;
                }
            }
        }
                
        if (density_assigned == true)
        {
            double MW; 
            MW = thermodynamicsMapXML.MolecularWeight_From_MassFractions(omega.GetHandle());
            if (temperature_assigned == true)	P_Pa	= rho*PhysicalConstants::R_J_kmol*T/MW;
            if (pressure_assigned == true)	T	= P_Pa*MW/PhysicalConstants::R_J_kmol/rho;
        }
        else
        {
            double MW;
            MW = thermodynamicsMapXML.MolecularWeight_From_MassFractions(omega.GetHandle());
            
            rho = (P_Pa*MW)/(PhysicalConstants::R_J_kmol*T);
        }
    }
}

#endif	/* OpenSMOKE_Grammar_Gas_H */

