#ifndef OpenSMOKE_Grammar_SolidStatus_H
#define OpenSMOKE_Grammar_SolidStatus_H

#include <string>
#include <boost/filesystem.hpp>
#include <dictionary/OpenSMOKE_Dictionary.h>
#include <dictionary/OpenSMOKE_DictionaryGrammar.h>

namespace BioSMOKE
{
class Grammar_SolidStatus : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
{
  protected:
    virtual void DefineRules()
    {
        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Temperature", OpenSMOKE::SINGLE_MEASURE,
                                                          "Temperature of the mixture (i.e. 500 K)", true));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Pressure", OpenSMOKE::SINGLE_MEASURE,
                                                          "Pressure of the mixture (i.e. 1 atm)", true));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
            "@Density", OpenSMOKE::SINGLE_MEASURE,
            "Skeletal Density of solid mixture (i.e. ASH 400 Kg/m3 CELL 850 Kg/m3)", true));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MoleFractions", OpenSMOKE::VECTOR_STRING_DOUBLE,
                                                          "Mole fractions of the solid mixture (i.e. ASH 0.60 C 0.40)",
                                                          true, "@MassFractions @Moles @Masses", "none", "none"));

        AddKeyWord(
            OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MassFractions", OpenSMOKE::VECTOR_STRING_DOUBLE,
                                                   "Mass fractions of the solid mixture (i.e. ASH 0.60 CELL 0.40)",
                                                   true, "@MoleFractions @Moles @Masses", "none", "none"));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Moles", OpenSMOKE::VECTOR_STRING_DOUBLE,
                                                          "Moles (relative) of the mixture (i.e. ASH 2 CELL 1)", true,
                                                          "@MoleFractions @MassFractions @Masses", "none", "none"));

        AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord(
            "@Masses", OpenSMOKE::VECTOR_STRING_DOUBLE, "Masses (relative) of the solid mixture (i.e. ASH 2 CELL 1)",
            true, "@MoleFractions @MassFractions @Moles ", "none", "none"));
    }
};

void GetSolidStatusFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary &dictionary,
                                  OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN &thermodynamicsSolidMapXML, double &T,
                                  double &P_Pa, double &rho, OpenSMOKE::OpenSMOKEVectorDouble &omega)
{
    Grammar_SolidStatus grammar_solid_status;
    dictionary.SetGrammar(grammar_solid_status);

    // Temperature
    {
        if (dictionary.CheckOption("@Temperature") == true)
        {
            double value;
            std::string units;
            dictionary.ReadMeasure("@Temperature", value, units);

            if (units == "K")
                T = value;
            else if (units == "C")
                T = value + 273.15;
            else
                OpenSMOKE::FatalErrorMessage("Unknown temperature units");
        }
    }

    // Pressure
    {
        if (dictionary.CheckOption("@Pressure") == true)
        {
            double value;
            std::string units;
            dictionary.ReadMeasure("@Pressure", value, units);

            if (units == "Pa")
                P_Pa = value;
            else if (units == "bar")
                P_Pa = value * 1.e5;
            else if (units == "atm")
                P_Pa = value * 101325.;
            else
                OpenSMOKE::FatalErrorMessage("Unknown pressure units");
        }
    }

    // Density Solid

    {
        if (dictionary.CheckOption("@Density") == true)
        {
            double value;
            std::string units;
            dictionary.ReadMeasure("@Density", value, units);

            if (units == "Kg/m3")
                rho = value;
            else
                OpenSMOKE::FatalErrorMessage("Unknown density units");
        }
    }

    

    // Composition internal
    {
        {
            std::vector<std::string> names;
            std::vector<double> values;

            if (dictionary.CheckOption("@MoleFractions") == true)
            {
                dictionary.ReadOption("@MoleFractions", names, values);
                const double sum = std::accumulate(values.begin(), values.end(), 0.);
                if (sum < (1. - 1e-6) || sum > (1. + 1e-6))
                    OpenSMOKE::FatalErrorMessage("The mole fractions must sum to 1.");

                for (unsigned int i = 0; i < thermodynamicsSolidMapXML.number_of_materials(); i++)
                {
                    OpenSMOKE::OpenSMOKEVectorDouble x(
                        thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size());
                    for (unsigned int j = 0; j < names.size(); j++)
                    {
                        int index = 0.;
                        for (unsigned int k = 0; k < thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size();
                             k++)
                        {
                            if (names[j] == thermodynamicsSolidMapXML.matrix_names_solid_species()[i][k])
                            {
                                index = k;
                                break;
                            }
                        }
                        x[index + 1] = values[j] / sum;
                    }
                    ChangeDimensions(thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size(), &omega, true);
                    double MW;
                    thermodynamicsSolidMapXML.MassFractions_From_MoleFractions(omega.GetHandle(), MW, x.GetHandle());
                }
            }

            else if (dictionary.CheckOption("@MassFractions") == true)
            {
                dictionary.ReadOption("@MassFractions", names, values);
                const double sum = std::accumulate(values.begin(), values.end(), 0.);
                if (sum < (1. - 1e-6) || sum > (1. + 1e-6))
                    OpenSMOKE::FatalErrorMessage("The mass fractions must sum to 1.");

                for (unsigned int i = 0; i < thermodynamicsSolidMapXML.number_of_materials(); i++)
                {
                    ChangeDimensions(thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size(), &omega, true);
                    for (unsigned int j = 0; j < names.size(); j++)
                    {
                        int index = 0.;
                        for (unsigned int k = 0; k < thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size();
                             k++)
                        {
                            if (names[j] == thermodynamicsSolidMapXML.matrix_names_solid_species()[i][k])
                            {
                                index = k;
                                break;
                            }
                        }
                        omega[index + 1] = values[j] / sum;
                    }
                }
            }

            else if (dictionary.CheckOption("@Moles") == true)
            {
                dictionary.ReadOption("@Moles", names, values);
                const double sum = std::accumulate(values.begin(), values.end(), 0.);

                for (unsigned int i = 0; i < thermodynamicsSolidMapXML.number_of_materials(); i++)
                {
                    OpenSMOKE::OpenSMOKEVectorDouble x(
                        thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size());
                    for (unsigned int j = 0; j < names.size(); j++)
                    {
                        int index = 0.;
                        for (unsigned int k = 0; k < thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size();
                             k++)
                        {
                            if (names[j] == thermodynamicsSolidMapXML.matrix_names_solid_species()[i][k])
                            {
                                index = k;
                                break;
                            }
                        }
                        x[index + 1] = values[j] / sum;
                    }

                    ChangeDimensions(thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size(), &omega, true);
                    double MW;
                    thermodynamicsSolidMapXML.MassFractions_From_MoleFractions(omega.GetHandle(), MW, x.GetHandle());
                }
            }
            else
            {
                dictionary.ReadOption("@Masses", names, values);
                const double sum = std::accumulate(values.begin(), values.end(), 0.);

                for (unsigned int i = 0; i < thermodynamicsSolidMapXML.number_of_materials(); i++)
                {
                    ChangeDimensions(thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size(), &omega, true);
                    for (unsigned int j = 0; j < names.size(); j++)
                    {
                        int index = 0.;
                        for (unsigned int k = 0; k < thermodynamicsSolidMapXML.matrix_names_solid_species()[i].size();
                             k++)
                        {
                            if (names[j] == thermodynamicsSolidMapXML.matrix_names_solid_species()[i][k])
                            {
                                index = k;
                                break;
                            }
                        }
                        omega[index + 1] = values[j] / sum;
                    }
                }
            }

            // Check if species is included in the solid mechanism
            for (int j = 0; j < names.size(); j++)
            {
                if (thermodynamicsSolidMapXML.IndexOfSpecies(names[j]) < 0)
                {
                    OpenSMOKE::FatalErrorMessage("Check the provided species");
                }
            }

            // Check if any of the composition is less than zero
            for (int j = 1; j <= omega.Size(); j++)
            {
                if (omega[j] < 0)
                {
                    OpenSMOKE::FatalErrorMessage("Check the provided composition");
                }
            }
        }
    }
}
} // namespace BioSMOKE
#endif /* OpenSMOKE_Grammar_SolidStatus_Options_H */
