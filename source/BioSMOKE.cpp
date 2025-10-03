/* ----------------------------------------------------------------------------------- *\
|                                                                                       |
|                  ____  _      ____  __  __  ___  _  _______                           |
|                 | __ )(_) ___/ ___||  \/  |/ _ \| |/ / ____|_ __  _ __                |
|                 |  _ \| |/ _ \___ \| |\/| | | | | ' /|  _| | '_ \| '_ \               |
|                 | |_) | | (_) |__) | |  | | |_| | . \| |___| |_) | |_) |              |
|                 |____/|_|\___/____/|_|  |_|\___/|_|\_\_____| .__/| .__/               |
|                                                            |_|   |_|                  |
|                                                                                       |
| ------------------------------------------------------------------------------------- |
|  See license and copyright at the end of this file.                                   |
| ------------------------------------------------------------------------------------- |
|                                                                                       |
|           Authors: Riccardo Caraccio <riccardo.caraccio@polimi.it>                    |
|                    Timoteo Dinelli   <timoteo.dinelli@polimi.it>                      |
|                    Andrea Locaspi    <andrea.locaspi@polimi.it>                       |
|                                                                                       |
|               CRECK Modeling Lab <https://www.creckmodeling.polimi.it>                |
|               Department of Chemistry, Materials and Chemical Engineering             |
|               Politecnico di Milano, P.zza Leonardo da Vinci 32, 20133 Milano         |
|                                                                                       |
\* ----------------------------------------------------------------------------------- */
#include <memory>
#include <string>

// ==================================================
// OpenSMOKEpp Headers
// ==================================================
#include <OpenSMOKEpp>
// ==================================================
// Maps
#include <maps/Maps_CHEMKIN>
#include <maps/ThermodynamicsMap_Solid_CHEMKIN.h>
// ==================================================
// Dictionaries
#include <dictionary/OpenSMOKE_DictionaryManager.h>
#include <dictionary/OpenSMOKE_DictionaryGrammar.h>
#include <dictionary/OpenSMOKE_DictionaryKeyWord.h>
// ==================================================
// Ideal Reactors utilities
#include <idealreactors/utilities/Utilities>
#include <idealreactors/plugflow/PlugFlowReactor_Profile.h>

// ==================================================
// Internal Headers
// ==================================================
#include "grammar/Grammar_BioSMOKE.h"
#include "grammar/Grammar_SolidStatus.h"
#include "grammar/Grammar_TGA_Biomass.h"
#include "grammar/Grammar_TotalSimulation_Biomass.h"

#include "TGAnalysis.h"

int main(int argc, char **argv)
{
    OpenSMOKE::OpenSMOKE_logo("BioSMOKEpp", "");

    std::string input_file_name_ = "input.dic";
    std::string main_dictionary_name_ = "BioSMOKE";

    // Program options from command line
    {
        namespace po = boost::program_options;
        po::options_description description("Options for the BioSMOKE numerical solver");
        description.add_options()("help", "print help messages")("np", po::value<unsigned int>(),
                                                                 "number of threads (default 1")(
            "input", po::value<std::string>(),
            "name of the file containing the main dictionary (default \"input.dic\")")(
            "dictionary", po::value<std::string>(), "name of the main dictionary to be used (default \"BioSMOKE\")");

        po::variables_map vm;
        try
        {
            po::store(po::parse_command_line(argc, argv, description), vm); // can throw

            if (vm.count("help"))
            {
                std::cout << "Basic Command Line Parameters" << std::endl;
                std::cout << description << std::endl;
                return OPENSMOKE_SUCCESSFULL_EXIT;
            }

            if (vm.count("input"))
            {
                input_file_name_ = vm["input"].as<std::string>();
            }

            if (vm.count("dictionary"))
            {
                main_dictionary_name_ = vm["dictionary"].as<std::string>();
            }

            po::notify(vm); // throws on error, so do after help in case  there are any problems
        }
        catch (po::error &e)
        {
            std::cerr << "Fatal error: " << e.what() << std::endl << std::endl;
            std::cerr << description << std::endl;
            return OPENSMOKE_FATAL_ERROR_EXIT;
        }
    }

    // Defines the grammar rules
    BioSMOKE::Grammar_BioSMOKE grammar_biosmoke;

    // Define the dictionaries
    OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
    dictionaries.ReadDictionariesFromFile(input_file_name_);
    dictionaries(main_dictionary_name_).SetGrammar(grammar_biosmoke);

    // Read kinetic model
    boost::filesystem::path kinetic_folder;
    if (dictionaries(main_dictionary_name_).CheckOption("@KineticsFolder") == true)
    {
        dictionaries(main_dictionary_name_).ReadPath("@KineticsFolder", kinetic_folder);
    }

    std::shared_ptr<OpenSMOKE::KineticsMap_CHEMKIN> kineticsMapXML;
    std::shared_ptr<OpenSMOKE::ThermodynamicsMap_CHEMKIN> thermodynamicsMapXML;
    std::shared_ptr<OpenSMOKE::TransportPropertiesMap_CHEMKIN> transportMapXML;
    std::shared_ptr<OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN> thermodynamicSolidMapXML;
    std::shared_ptr<OpenSMOKE::KineticsMap_Solid_CHEMKIN> kineticsSolidMapXML;

    // Read the kinetic scheme in XML format
    {
        if (!boost::filesystem::exists(kinetic_folder / "kinetics.xml"))
        {
            OpenSMOKE::FatalErrorMessage("The kinetic mechanism does not exist! Please check");
        }

        boost::property_tree::ptree ptree;
        boost::property_tree::read_xml((kinetic_folder / "kinetics.xml").string(), ptree);

        double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
        thermodynamicsMapXML = std::make_shared<OpenSMOKE::ThermodynamicsMap_CHEMKIN>(ptree);
        kineticsMapXML = std::make_shared<OpenSMOKE::KineticsMap_CHEMKIN>(*thermodynamicsMapXML, ptree);
        transportMapXML = std::make_shared<OpenSMOKE::TransportPropertiesMap_CHEMKIN>(ptree);
        double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
        std::cout << "Time to read Gas Phase XML file: " << tEnd - tStart << std::endl;

        // Read the solid-phase reactions
        std::cout << "Reading the kinetic scheme of the solid phase in XML format" << std::endl;

        if (!boost::filesystem::exists(kinetic_folder / "kinetics.solid.xml"))
        {
            OpenSMOKE::FatalErrorMessage("The solid-phase kinetic mechanism does not exist! Please check");
        }
        else
        {
            boost::property_tree::ptree ptree;
            boost::property_tree::read_xml((kinetic_folder / "kinetics.solid.xml").string(), ptree);

            double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
            thermodynamicSolidMapXML = std::make_shared<OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN>(ptree);
            kineticsSolidMapXML =
                std::make_shared<OpenSMOKE::KineticsMap_Solid_CHEMKIN>(*thermodynamicSolidMapXML, ptree, 1);
            double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
            std::cout << "Time to read Solid Phase XML file: " << tEnd - tStart << std::endl;
        }
    }

    std::string analysis_type;
    std::vector<std::string> output_species;
    BioSMOKE::Analysis_Type type;

    // TGA Analysis
    double heating_rate, final_time;

    // Total Analysis
    bool energy_balance, volume_loss;
    double porosity, initial_radius, Da_number, ext_heat_transf_coeff, lambda_solid;
    int number_of_layers;

    if (dictionaries(main_dictionary_name_).CheckOption("@Type") == true)
    {

        dictionaries(main_dictionary_name_).ReadString("@Type", analysis_type);

        if (analysis_type == "TG_Analysis")
        {
            std::string name_of_solid_status_subdictionary;
            if (dictionaries(main_dictionary_name_).CheckOption("@TGA"), true)
            {
                dictionaries(main_dictionary_name_).ReadDictionary("@TGA", name_of_solid_status_subdictionary);
            }

            BioSMOKE::Get_TGAanalysisFromDictionary(dictionaries(name_of_solid_status_subdictionary), heating_rate,
                                                    final_time, output_species);

            // TODO
            // if (output_species_[0] == "all")
            // {
            //     output_species_ = kineticsSolidMapXML->NamesOfSpecies();
            // }
            type = BioSMOKE::THERMOGRAVIMETRIC_ANALYSIS;
        }
        else if (analysis_type == "Total_Analysis")
        {
            std::string name_of_solid_status_subdictionary;
            if (dictionaries(main_dictionary_name_).CheckOption("@Total_Analysis") == true)
            {
                dictionaries(main_dictionary_name_)
                    .ReadDictionary("@Total_Analysis", name_of_solid_status_subdictionary);
            }

            BioSMOKE::Get_TotalSimulation_analysisFromDictionary(
                dictionaries(name_of_solid_status_subdictionary), energy_balance, volume_loss, final_time, porosity,
                number_of_layers, initial_radius, Da_number, ext_heat_transf_coeff, lambda_solid, output_species);

            type = BioSMOKE::ONE_DIMENSIONAL_SPHERICAL_PARTICLE;
        }
        else
        {
            std::cout << "Unknown simulation type: " << analysis_type << std::endl;
            OpenSMOKE::FatalErrorMessage("Avaiable options: TG_Analysis or Total_Analysis");
        }
    }

    std::string name_of_profile_subdictionary;
    bool is_temperature_profile = false;
    OpenSMOKE::PlugFlowReactor_Profile *temperature_profile; // not really clean
    if (dictionaries(main_dictionary_name_).CheckOption("@TemperatureProfile") == true)
    {
        dictionaries(main_dictionary_name_).ReadDictionary("@TemperatureProfile", name_of_profile_subdictionary);

        OpenSMOKE::OpenSMOKEVectorDouble x, y;
        std::string x_variable, y_variable;

        GetXYProfileFromDictionary(dictionaries(name_of_profile_subdictionary), x, y, x_variable, y_variable);
        is_temperature_profile = true;
        temperature_profile = new OpenSMOKE::PlugFlowReactor_Profile(x, y, x_variable);
    }

    double T_gas, P_Pa_gas;
    OpenSMOKE::OpenSMOKEVectorDouble omega0_gas;
    // Read initial conditions
    {
        std::string name_of_gas_status_subdictionary;
        if (dictionaries(main_dictionary_name_).CheckOption("@InletGasStatus") == true)
        {
            dictionaries(main_dictionary_name_).ReadDictionary("@InletGasStatus", name_of_gas_status_subdictionary);
        }
        else
        {
            OpenSMOKE::FatalErrorMessage("Initial gas status not defined");
        }

        GetGasStatusFromDictionary(dictionaries(name_of_gas_status_subdictionary), *thermodynamicsMapXML, T_gas,
                                   P_Pa_gas, omega0_gas);

        // TODO
        // if (output_species_[0] == "all")
        // {
        //     output_species_ = kineticsSolidMapXML->NamesOfSpecies();
        // }
    }

    double T_solid, P_Pa_solid, rho_solid;
    OpenSMOKE::OpenSMOKEVectorDouble omega0_solid;
    {
        std::string name_of_solid_status_subdictionary;
        if (dictionaries(main_dictionary_name_).CheckOption("@InletSolidStatus") == true)
        {
            dictionaries(main_dictionary_name_).ReadDictionary("@InletSolidStatus", name_of_solid_status_subdictionary);
        }
        else
        {
            OpenSMOKE::FatalErrorMessage("Initial solid status not defined");
        }

        BioSMOKE::GetSolidStatusFromDictionary(dictionaries(name_of_solid_status_subdictionary),
                                               *thermodynamicSolidMapXML, T_solid, P_Pa_solid, rho_solid, omega0_solid);

        if (is_temperature_profile == true)
        {
            T_solid = temperature_profile->Get(0.);
        }
    }

    std::shared_ptr<OpenSMOKE::ODE_Parameters> ode_parameters; // TODO

    // Convert OpenSMOKE::OpenSMOKEVectorDouble to std::vector<double>
    std::vector<double> omega0gas(omega0_gas.GetHandle(), omega0_gas.GetHandle() + omega0_gas.Size());
    std::vector<double> omega0solid(omega0_solid.GetHandle(), omega0_solid.GetHandle() + omega0_solid.Size());

    if (type == BioSMOKE::THERMOGRAVIMETRIC_ANALYSIS)
    {
        // clang-format off
        BioSMOKE::TGAnalysis tga_analysis(  *thermodynamicsMapXML,
                                            *kineticsMapXML,
                                            *transportMapXML,
                                            *thermodynamicSolidMapXML,
                                            *kineticsSolidMapXML,
                                            *ode_parameters,
                                            T_solid,
                                            P_Pa_solid,
                                            rho_solid,
                                            omega0gas,
                                            omega0solid,
                                            heating_rate);
        tga_analysis.Solve(0., final_time);
        // clang-format on
    }

    return OPENSMOKE_SUCCESSFULL_EXIT;
}
