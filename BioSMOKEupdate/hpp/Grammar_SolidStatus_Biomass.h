#ifndef OpenSMOKE_Grammar_SolidStatusGG_H
#define	OpenSMOKE_Grammar_SolidStatusGG_H

#include <string>
#include "boost/filesystem.hpp"
#include "dictionary/OpenSMOKE_Dictionary.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"

namespace OpenSMOKE
{
    class Grammar_SolidStatusGG : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
    {
    protected:
        virtual void DefineRules()
        {
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Temperature", 
                                                                OpenSMOKE::SINGLE_MEASURE, 
                                                                "Temperature of the mixture (i.e. 500 K)", 
                                                                true) );
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@TemperatureSup", 
                                                                OpenSMOKE::SINGLE_MEASURE, 
                                                                "Bulk gas temperature, employed for convective heat transfer evaluation (i.e. 500 K)", 
                                                                true) );
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@hext", 
                                                                OpenSMOKE::SINGLE_MEASURE, 
                                                                "Convective heat exchange coefficient (i.e. 20 W/m2/K)", 
                                                                true) );
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Pressure", 
                                                                OpenSMOKE::SINGLE_MEASURE, 
                                                                "Pressure of the mixture (i.e. 1 atm)", 
                                                                true) );	
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Density", 
                                                                OpenSMOKE::SINGLE_MEASURE,
                                                                "Skeletal Density of solid mixture (i.e. ASH 400 Kg/m3 CELL 850 Kg/m3)", 
                                                                true) );	

            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MoleFractions", 
                                                                OpenSMOKE::VECTOR_STRING_DOUBLE, 
                                                                "Mole fractions of the solid mixture (i.e. ASH 0.60 C 0.40)", 
                                                                true,
                                                                "@MassFractions @Moles @Masses",
                                                                "none",
                                                                "none") );

            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MassFractions", 
                                                                OpenSMOKE::VECTOR_STRING_DOUBLE, 
                                                                "Mass fractions of the solid mixture (i.e. ASH 0.60 CELL 0.40)", 
                                                                true,
                                                                "@MoleFractions @Moles @Masses",
                                                                "none",
                                                                "none") );

            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Moles", 
                                                                OpenSMOKE::VECTOR_STRING_DOUBLE, 
                                                                "Moles (relative) of the mixture (i.e. ASH 2 CELL 1)", 
                                                                true,
                                                                "@MoleFractions @MassFractions @Masses",
                                                                "none",
                                                                "none") );

            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Masses", 
                                                                OpenSMOKE::VECTOR_STRING_DOUBLE, 
                                                                "Masses (relative) of the solid mixture (i.e. ASH 2 CELL 1)", 
                                                                true,
                                                                "@MoleFractions @MassFractions @Moles ",
                                                                "none",
                                                                "none") );
            
            /*AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MoleFractionsSup", 
                                                                OpenSMOKE::VECTOR_STRING_DOUBLE, 
                                                                "Mole fractions of the solid mixture at surface (i.e. ASH 0.60 C 0.40)", 
                                                                true,
                                                                "@MassFractionsSup @MolesSup @MassesSup",
                                                                "none",
                                                                "none") );

            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MassFractionsSup", 
                                                                OpenSMOKE::VECTOR_STRING_DOUBLE, 
                                                                "Mass fractions of the solid mixture at surface (i.e. ASH 0.60 CELL 0.40)", 
                                                                true,
                                                                "@MoleFractionsSup @MolesSup @MassesSup",
                                                                "none",
                                                                "none") );

            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MolesSup", 
                                                                OpenSMOKE::VECTOR_STRING_DOUBLE, 
                                                                "Moles (relative) of the mixture at surface (i.e. ASH 2 CELL 1)", 
                                                                true,
                                                                "@MoleFractionsSup @MassFractionsSup @MassesSup",
                                                                "none",
                                                                "none") );

            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MassesSup", 
                                                                OpenSMOKE::VECTOR_STRING_DOUBLE, 
                                                                "Masses (relative) of the solid mixture at surface (i.e. ASH 2 CELL 1)", 
                                                                true,
                                                                "@MoleFractionsSup @MassFractionsSup @MolesSup ",
                                                                "none",
                                                                "none") );*/
            
            /*AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@lambda_A", 
                                                                OpenSMOKE::VECTOR_STRING_DOUBLE, 
                                                                "Thermal conductivity for the solid [A term]--> lamda a+b*T", 
                                                                true ));
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@lambda_B", 
                                                                OpenSMOKE::VECTOR_STRING_DOUBLE, 
                                                                "Thermal conductivity for the solid [B term]--> lamda a+b*T", 
                                                                true ));*/
            
            AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@lambda", 
                                                                OpenSMOKE::SINGLE_MEASURE, 
                                                                "Thermal conductivity for the solid", 
                                                                true) );

			
            
        }
    };
                                                                                  
    void GetSolidStatusFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN& thermodynamicsSolidMapXML,
                                      double& T, double& P_Pa, double& Ts, double& rho,OpenSMOKE::OpenSMOKEVectorDouble& omega, 
                                      OpenSMOKE::OpenSMOKEVectorDouble& omegaSup,OpenSMOKE::OpenSMOKEVectorDouble& lambda_A,OpenSMOKE::OpenSMOKEVectorDouble& lambda_B,double& hext, double& lambda_S)
    
    
    {
        Grammar_SolidStatusGG grammar_solid_status;
        dictionary.SetGrammar(grammar_solid_status);

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
            }
        }
        
        // Temperature at surface
        {
            if (dictionary.CheckOption("@TemperatureSup") == true)
            {
                double value;
                std::string units;
                dictionary.ReadMeasure("@TemperatureSup", value, units);

                if (units == "K")			Tsup_solid = value;
                else if (units == "C")		Tsup_solid = value + 273.15;
                else OpenSMOKE::FatalErrorMessage("Unknown temperature units");
            }
        }
        
        // Heat exchange coefficient
        {
            if (dictionary.CheckOption("@hext") == true)
            {
                double value;
                std::string units;
                dictionary.ReadMeasure("@hext", value, units);

                if (units == "W/m2/K")			hext = value;
                else OpenSMOKE::FatalErrorMessage("Unknown temperature units");
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

                
            }
        }

        
        // Density Solid
        
        {
            if (dictionary.CheckOption("@Density") == true)
            {
                double value;
                std::string units;
                dictionary.ReadMeasure("@Density", value, units);

                if (units == "Kg/m3")	        rho = value;
                else OpenSMOKE::FatalErrorMessage("Unknown density units");

                
            }
        }
        
        
         // Heat exchange coefficient
        {
            if (dictionary.CheckOption("@lambda") == true)
            {
                double value;
                std::string units;
                dictionary.ReadMeasure("@lambda", value, units);

                if (units == "W/m/K")			lambdaS = value;
                else OpenSMOKE::FatalErrorMessage("Unknown temperature units");
            }
        }  
        
        {
            std::vector<std::string> names;
            std::vector<double> values;
            
            /*if (dictionary.CheckOption("@lambda_A") == true)
            {
                dictionary.ReadOption("@lambda_A", names, values);
                
                for(unsigned int i=0;i<thermodynamicsSolidMapXML.number_of_materials();i++)
		{          
                    ChangeDimensions(thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size(), &lambda_A, true);
                    
                    for(unsigned int j=0;j<names.size();j++)
                    {
                        int index = 0.;
                        for(unsigned int k=0;k<thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size();k++)
                        {
                            if (names[j] == thermodynamicsSolidMapXML.matrix_names_solid_species()[i][k])
                            {
                                index = k;
                                break;
                            }
                            
                        }
                        lambda_A[index+1] = values[j];
                    }
                }
            }*/
        }
        
        // Thermal conductivity -- B term
        
        {
            /*std::vector<std::string> names;
            std::vector<double> values;
            
            if (dictionary.CheckOption("@lambda_B") == true)
            {
                dictionary.ReadOption("@lambda_B", names, values);
                
                for(unsigned int i=0;i<thermodynamicsSolidMapXML.number_of_materials();i++)
		{          
                    ChangeDimensions(thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size(), &lambda_B, true);
                    for(unsigned int j=0;j<names.size();j++)
                    {
                        int index = 0.;
                        for(unsigned int k=0;k<thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size();k++)
                        {
                            if (names[j] == thermodynamicsSolidMapXML.matrix_names_solid_species()[i][k])
                            {
                                index = k;
                                break;
                            }
                            
                        }
                        lambda_B[index+1] = values[j];
                    }
                }
            }*/
        }

            
        // Composition internal
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
                    
                    for(unsigned int i=0;i<thermodynamicsSolidMapXML.number_of_materials();i++)
                    {
                        OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size());
                        for(unsigned int j=0;j<names.size();j++)
                        {
                            int index = 0.;
                            for(unsigned int k=0;k<thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size();k++)
                            {
                                if (names[j] == thermodynamicsSolidMapXML.matrix_names_solid_species()[i][k])
                                {
                                    index = k;
                                    break;
                                }
                            }
                            x[index+1] = values[j]/sum;
                        }
                        ChangeDimensions(thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size(), &omega, true);
                        double MW;
                        thermodynamicsSolidMapXML.MassFractions_From_MoleFractions(omega.GetHandle(),MW,x.GetHandle());
                        
                    }
                }
                
                else if (dictionary.CheckOption("@MassFractions") == true)
                {
                    dictionary.ReadOption("@MassFractions", names, values);
                    const double sum =std::accumulate(values.begin(),values.end(),0.);
                    if (sum<(1.-1e-6) || sum>(1.+1e-6))
                            OpenSMOKE::FatalErrorMessage("The mass fractions must sum to 1.");
                    
                    for(unsigned int i=0;i<thermodynamicsSolidMapXML.number_of_materials();i++)
                    {
                        ChangeDimensions(thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size(), &omega, true);
                        for(unsigned int j=0;j<names.size();j++)
                        {
                            int index = 0.;
                            for(unsigned int k=0;k<thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size();k++)
                            {
                                if (names[j] == thermodynamicsSolidMapXML.matrix_names_solid_species()[i][k])
                                {
                                    index = k;
                                    break;
                                }

                            }
                            omega[index+1] = values[j]/sum;
                        }
                    }
                }
                
                else if (dictionary.CheckOption("@Moles") == true)
                {
                    dictionary.ReadOption("@Moles", names, values);
                    const double sum =std::accumulate(values.begin(),values.end(),0.);
                    
                    for(unsigned int i=0;i<thermodynamicsSolidMapXML.number_of_materials();i++)
                    {
                        OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size());
                        for(unsigned int j=0;j<names.size();j++)
                        {
                            int index = 0.;
                            for(unsigned int k=0;k<thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size();k++)
                            {
                                if (names[j] == thermodynamicsSolidMapXML.matrix_names_solid_species()[i][k])
                                {
                                    index = k;
                                    break;
                                }
                            }
                            x[index+1] = values[j]/sum;
                        }
                        
                        ChangeDimensions(thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size(), &omega, true);
                        double MW;
                        thermodynamicsSolidMapXML.MassFractions_From_MoleFractions(omega.GetHandle(),MW,x.GetHandle());
                    }
                    
                }
                else  
                {
                    dictionary.ReadOption("@Masses", names, values);
                    const double sum =std::accumulate(values.begin(),values.end(),0.);
                    
                    for(unsigned int i=0;i<thermodynamicsSolidMapXML.number_of_materials();i++)
                    {
                        ChangeDimensions(thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size(), &omega, true);
                        for(unsigned int j=0;j<names.size();j++)
                        {
                            int index = 0.;
                            for(unsigned int k=0;k<thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size();k++)
                            {
                                if (names[j] == thermodynamicsSolidMapXML.matrix_names_solid_species()[i][k])
                                {
                                    index = k;
                                    break;
                                }

                            }
                            omega[index+1] = values[j]/sum;
                        }
                    }
                }
            }

            
        }
        
        // Composition at surface
        /*{
            {
                std::vector<std::string> names;
                std::vector<double> values;

                if (dictionary.CheckOption("@MoleFractionsSup") == true)
                {
                    dictionary.ReadOption("@MoleFractionsSup", names, values);
                    const double sum =std::accumulate(values.begin(),values.end(),0.);
                    if (sum<(1.-1e-6) || sum>(1.+1e-6))
                            OpenSMOKE::FatalErrorMessage("The mole fractions must sum to 1.");
                    
                    for(unsigned int i=0;i<thermodynamicsSolidMapXML.number_of_materials();i++)
                    {
                        OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size());
                        for(unsigned int j=0;j<names.size();j++)
                        {
                            int index = 0.;
                            for(unsigned int k=0;k<thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size();k++)
                            {
                                if (names[j] == thermodynamicsSolidMapXML.matrix_names_solid_species()[i][k])
                                {
                                    index = k;
                                    break;
                                }
                            }
                            x[index+1] = values[j]/sum;
                        }
                        
                        ChangeDimensions(thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size(), &omegaSup, true);
                        double MW;
                        thermodynamicsSolidMapXML.MassFractions_From_MoleFractions(omegaSup,MW,x);
                    }
                }
                
                else if (dictionary.CheckOption("@MassFractionsSup") == true)
                {
                    dictionary.ReadOption("@MassFractionsSup", names, values);
                    const double sum =std::accumulate(values.begin(),values.end(),0.);
                    if (sum<(1.-1e-6) || sum>(1.+1e-6))
                            OpenSMOKE::FatalErrorMessage("The mass fractions must sum to 1.");
                    
                    for(unsigned int i=0;i<thermodynamicsSolidMapXML.number_of_materials();i++)
                    {
                        ChangeDimensions(thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size(), &omegaSup, true);
                        for(unsigned int j=0;j<names.size();j++)
                        {
                            int index = 0.;
                            for(unsigned int k=0;k<thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size();k++)
                            {
                                if (names[j] == thermodynamicsSolidMapXML.matrix_names_solid_species()[i][k])
                                {
                                    index = k;
                                    break;
                                }

                            }
                            omegaSup[index+1] = values[j]/sum;
                        }
                    }
                    
                }
                
                else if (dictionary.CheckOption("@MolesSup") == true)
                {
                    dictionary.ReadOption("@MolesSup", names, values);
                    const double sum =std::accumulate(values.begin(),values.end(),0.);
                    
                    for(unsigned int i=0;i<thermodynamicsSolidMapXML.number_of_materials();i++)
                    {
                        OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size());
                        for(unsigned int j=0;j<names.size();j++)
                        {
                            int index = 0.;
                            for(unsigned int k=0;k<thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size();k++)
                            {
                                if (names[j] == thermodynamicsSolidMapXML.matrix_names_solid_species()[i][k])
                                {
                                    index = k;
                                    break;
                                }
                            }
                            x[index+1] = values[j]/sum;
                        }
                        
                        ChangeDimensions(thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size(), &omegaSup, true);
                        double MW;
                        thermodynamicsSolidMapXML.MassFractions_From_MoleFractions(omegaSup,MW,x);
                        
                    }
                    
                }
                else  
                {
                    dictionary.ReadOption("@MassesSup", names, values);
                    const double sum =std::accumulate(values.begin(),values.end(),0.);
                    
                    for(unsigned int i=0;i<thermodynamicsSolidMapXML.number_of_materials();i++)
                    {
                        ChangeDimensions(thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size(), &omegaSup, true);
                        for(unsigned int j=0;j<names.size();j++)
                        {
                            int index = 0.;
                            for(unsigned int k=0;k<thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size();k++)
                            {
                                if (names[j] == thermodynamicsSolidMapXML.matrix_names_solid_species()[i][k])
                                {
                                    index = k;
                                    break;
                                }

                            }
                            omegaSup[index+1] = values[j]/sum;
                        }
                    }
                    
                }
            }
        }*/


    }
    
    void GetSolidStatusFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary, double& T, double& P_Pa, 
                                      double& Ts,OpenSMOKE::OpenSMOKEVectorDouble& rho,OpenSMOKE::OpenSMOKEVectorDouble& omega,OpenSMOKE::OpenSMOKEVectorDouble& omegaSup)
    {
        Grammar_SolidStatusGG grammar_solid_status;
        dictionary.SetGrammar(grammar_solid_status);

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

                
            }
        }
        
        // Temperature at surface
        {
            if (dictionary.CheckOption("@TemperatureSup") == true)
            {
                double value;
                std::string units;
                dictionary.ReadMeasure("@TemperatureSup", value, units);

                if (units == "K")			Ts = value;
                else if (units == "C")		Ts = value + 273.15;
                else OpenSMOKE::FatalErrorMessage("Unknown temperature units");
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

                
            }
        }

        // Density
        
        {
            std::vector<std::string> names;
            std::vector<double> values;
            
            if (dictionary.CheckOption("@Density") == true)
            {
                dictionary.ReadOption("@Density", names, values);
                ChangeDimensions(names.size()+1, &rho, true);
                for(unsigned int i=0;i<names.size();i++)
                    rho[i+1] = values[i];
            }
        }

            
        // Composition
        {
            {
                std::vector<std::string> names;
                std::vector<double> values;

                if (dictionary.CheckOption("@MassFractions") == true)
                {
                    dictionary.ReadOption("@MassFractions", names, values);

                    const double sum =std::accumulate(values.begin(),values.end(),0.);
                    if (sum<(1.-1e-6) || sum>(1.+1e-6))
                            OpenSMOKE::FatalErrorMessage("The mass fractions must sum to 1.");

                    ChangeDimensions(names.size()+1, &omega, true);
                    for(unsigned int i=0;i<names.size();i++)
                            omega[i+1] = values[i]/sum;
                }
                
                else 
                {
                    dictionary.ReadOption("@Masses", names, values);
                    const double sum =std::accumulate(values.begin(),values.end(),0.);
                    ChangeDimensions(names.size()+1, &omega, true);
                    for(unsigned int i=0;i<names.size();i++)
                            omega[i+1] = values[i]/sum;
                }
            }
        }
        
        // Composition at surface
        
        {
            {
                std::vector<std::string> names;
                std::vector<double> values;

                if (dictionary.CheckOption("@MassFractionsSup") == true)
                {
                    dictionary.ReadOption("@MassFractionsSup", names, values);

                    const double sum =std::accumulate(values.begin(),values.end(),0.);
                    if (sum<(1.-1e-6) || sum>(1.+1e-6))
                            OpenSMOKE::FatalErrorMessage("The mass fractions must sum to 1.");

                    ChangeDimensions(names.size()+1, &omegaSup, true);
                    for(unsigned int i=0;i<names.size();i++)
                            omegaSup[i+1] = values[i]/sum;
                }
                
                else 
                {
                    dictionary.ReadOption("@MassesSup", names, values);
                    const double sum =std::accumulate(values.begin(),values.end(),0.);
                    ChangeDimensions(names.size()+1, &omegaSup, true);
                    for(unsigned int i=0;i<names.size();i++)
                            omegaSup[i+1] = values[i]/sum;
                }
            }
        }
    }
}

#endif	/* OpenSMOKE_Grammar_SolidStatus_Options_H */

