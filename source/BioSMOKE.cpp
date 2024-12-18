#include <memory>
#include <string>

// OpenSMOKE
#include <OpenSMOKEpp>
// CHEMKIN maps
#include <maps/Maps_CHEMKIN>
#include <maps/ThermodynamicsMap_Solid_CHEMKIN.h>

// Reactor utilities
#include <idealreactors/utilities/Utilities>
#include <surfacereactors/utilities/Utilities>
#include <utilities/ropa/OnTheFlyROPA.h>

// Numerical parameters
#include <math/native-dae-solvers/parameters/DaeSolver_Parameters.h>
#include <math/native-nls-solvers/NonLinearSolver.h>
#include <math/native-nls-solvers/KernelSparse.h>
#include <math/native-ode-solvers/MultiValueSolver>
#include <math/external-ode-solvers/ODE_Parameters.h>

// Internal headers
#include "grammar/Grammar_BioSMOKE.h"
#include "grammar/Grammar_SolidStatus.h"

int main(int argc, char **argv)
{
    OpenSMOKE::OpenSMOKE_logo("BioSMOKEpp", "Riccardo Caraccio (riccardo.caraccio@polimi.it)");

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
            "dictionary", po::value<std::string>(),
            "name of the main dictionary to be used (default \"BatchReactor\")");

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
                input_file_name_ = vm["input"].as<std::string>();

            if (vm.count("dictionary"))
                main_dictionary_name_ = vm["dictionary"].as<std::string>();

            po::notify(vm); // throws on error, so do after help in case  there are any problems
        }
        catch (po::error &e)
        {
            std::cerr << "Fatal error: " << e.what() << std::endl << std::endl;
            std::cerr << description << std::endl;
            return OPENSMOKE_FATAL_ERROR_EXIT;
        }
    }

    // Define the dictionaries
    OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
    BioSMOKE::Grammar_BioSMOKE grammar_biosmoke;
    dictionaries.ReadDictionariesFromFile(input_file_name_);
    dictionaries(main_dictionary_name_).SetGrammar(grammar_biosmoke);

    // Read kinetic model
    boost::filesystem::path kinetic_folder;
    if (dictionaries(main_dictionary_name_).CheckOption("@KineticsFolder") == true)
        dictionaries(main_dictionary_name_).ReadPath("@KineticsFolder", kinetic_folder);

    std::shared_ptr<OpenSMOKE::KineticsMap_CHEMKIN> kinetics_map_XML;
    std::shared_ptr<OpenSMOKE::ThermodynamicsMap_CHEMKIN> thermodynamicsMapXML;
    std::shared_ptr<OpenSMOKE::TransportPropertiesMap_CHEMKIN> transportMapXML;
    std::shared_ptr<OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN> thermodynamic_solid_map_XML;
    std::shared_ptr<OpenSMOKE::KineticsMap_Solid_CHEMKIN> kinetics_solid_map_XML;

    // Read the kinetic scheme in XML format
    {
        if (!boost::filesystem::exists(kinetic_folder / "kinetics.xml"))
            OpenSMOKE::FatalErrorMessage("The kinetic mechanism does not exist! Please check");

        boost::property_tree::ptree ptree;
        boost::property_tree::read_xml((kinetic_folder / "kinetics.xml").string(), ptree);

        double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
        thermodynamicsMapXML = std::make_shared<OpenSMOKE::ThermodynamicsMap_CHEMKIN>(ptree);
        kinetics_map_XML = std::make_shared<OpenSMOKE::KineticsMap_CHEMKIN>(*thermodynamicsMapXML, ptree);
        transportMapXML = std::make_shared<OpenSMOKE::TransportPropertiesMap_CHEMKIN>(ptree);
        double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
        std::cout << "Time to read Gas Phase XML file: " << tEnd - tStart << std::endl;

        // Read the solid-phase reactions
        std::cout << "Reading the kinetic scheme (solid phase) in XML format" << std::endl;

        if (!boost::filesystem::exists(kinetic_folder / "kinetics.solid.xml"))
        {
            std::cout << "No reactions in the solid phase" << std::endl;
        }
        else
        {
            boost::property_tree::ptree ptree;
            boost::property_tree::read_xml((kinetic_folder / "kinetics.solid.xml").string(), ptree);

            double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
            thermodynamic_solid_map_XML = std::make_shared<OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN>(ptree);
            kinetics_solid_map_XML =
                std::make_shared<OpenSMOKE::KineticsMap_Solid_CHEMKIN>(*thermodynamic_solid_map_XML, ptree, 1);
            double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
            std::cout << "Time to read Solid Phase XML file: " << tEnd - tStart << std::endl;
        }
    }

    std::string analysis_type;
    if (dictionaries(main_dictionary_name_).CheckOption("@Type") == true)
    {

        dictionaries(main_dictionary_name_).ReadString("@Type", analysis_type);

        if (analysis_type == "TG_Analysis")
        {
            std::string name_of_solid_status_subdictionary;
            if (dictionaries(main_dictionary_name_).CheckOption("@TGA"), true)
                dictionaries(main_dictionary_name_).ReadDictionary("@TGA", name_of_solid_status_subdictionary);

            // TODO
            // Get_TGAanalysisFromDictionary(dictionaries(name_of_solid_status_subdictionary), analysis, Heat_Rates,
            //                               final_time, output_species_);
            //
            // if (output_species_[0] == "all")
            // {
            //     output_species_ = kinetics_solid_map_XML->NamesOfSpecies();
            // }
        }
        else if (analysis_type == "Total_Analysis")
        {
            std::string name_of_solid_status_subdictionary;
            if (dictionaries(main_dictionary_name_).CheckOption("@Total_Analysis") == true)
                dictionaries(main_dictionary_name_)
                    .ReadDictionary("@Total_Analysis", name_of_solid_status_subdictionary);

            // TODO
            // Get_TotalSimulation_analysisFromDictionary(dictionaries(name_of_solid_status_subdictionary), analysis,
            //                                            energyBalance, volumeLoss, final_time, epsi, tau, intervalli,
            //                                            raggioTot, Da, output_species_);
        }
        else
        {
            std::cout << "Unknown";
        }
    }

    std::string name_of_gas_status_subdictionary;
    if (dictionaries(main_dictionary_name_).CheckOption("@TemperatureProfile") == true)
    {
        dictionaries(main_dictionary_name_).ReadDictionary("@TemperatureProfile", name_of_gas_status_subdictionary);

        std::vector<double> x, y;
        std::string x_variable, y_variable;

        // TODO
        // GetXYProfileFromDictionary(dictionaries(name_of_gas_status_subdictionary), x, y, x_variable, y_variable);
        // temperature_profile = true;
        // vector_profile = new OpenSMOKE::BiomassTemperature_Profile(x, y);
        // temperature_profile_ = vector_profile;
        // final_time = x[x.Size()];
    }

    double T_gas, P_Pa_gas;
    OpenSMOKE::OpenSMOKEVectorDouble omega;
    // Read initial conditions
    {
        std::string name_of_gas_status_subdictionary;
        if (dictionaries(main_dictionary_name_).CheckOption("@InitialGasStatus") == true)
            dictionaries(main_dictionary_name_).ReadDictionary("@InitialGasStatus", name_of_gas_status_subdictionary);

        GetGasStatusFromDictionary(dictionaries(name_of_gas_status_subdictionary), *thermodynamicsMapXML, T_gas,
                                   P_Pa_gas, omega);
    }

    double T_solid, P_solid, rho_solid, hext, lambda_solid;
    OpenSMOKE::OpenSMOKEVectorDouble omega_solid;
    {
        std::string name_of_solid_status_subdictionary;
        if (dictionaries(main_dictionary_name_).CheckOption("@InitialSolidStatus") == true)
            dictionaries(main_dictionary_name_)
                .ReadDictionary("@InitialSolidStatus", name_of_solid_status_subdictionary);

        BioSMOKE::GetSolidStatusFromDictionary(dictionaries(name_of_solid_status_subdictionary),
                                               *thermodynamic_solid_map_XML, T_solid, P_solid, rho_solid, omega_solid,
                                               hext, lambda_solid);
        // TODO
        // if (temperature_profile == true)
        //     Tsolid = temperature_profile_->Get(0.);
    }

    return 0;
}
