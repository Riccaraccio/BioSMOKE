/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|   License                                                               |
|                                                                         |
|   Copyright(C) 2016-2012  Alberto Cuoci and Alessandro Stagni           |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
#include <sys/stat.h>
#include <sys/types.h>

// Mattia 2017 Conversione LINUX-->WIN
#ifdef _WIN32
#include <stdlib.h>
#endif

// #ifdef linux
// #include <dirent.h>
// #endif

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"
#include "Eigen/Dense"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// Reactor utilities
#include "idealreactors/utilities/Utilities"
#include "surfacereactors/utilities/Utilities"
#include "utilities/ropa/OnTheFlyROPA.h"

// Numerical parameters
#include "math/native-dae-solvers/parameters/DaeSolver_Parameters.h"
#include "math/native-nls-solvers/NonLinearSolver.h"
#include "math/native-nls-solvers/KernelSparse.h"
#include "math/native-ode-solvers/MultiValueSolver"
#include "OpenSMOKE_ODESystem_Interface.h"
#include "math/external-ode-solvers/ODE_Parameters.h"

#include "BiomassTemperature_Profile.h"
#include "memoryAllocation.H"

OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML;
OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML;
OpenSMOKE::TransportPropertiesMap_CHEMKIN* transportMapXML;
OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN* thermodynamicsSolidMapXML;
std::vector<OpenSMOKE::KineticsMap_Solid_CHEMKIN*> kineticsSolidMapXML;

#include "hpp/Grammar_biomass.h"
#include "hpp/Grammar_GasStatus_Biomass.h"
#include "hpp/Grammar_SolidStatus_Biomass.h"
#include "hpp/Grammar_TGA_Biomass.h"
#include "hpp/Grammar_TotalSimulation_Biomass.h"

// License
// #include "license/OpenSMOKELicenseUtilities.hpp"

using namespace std;    
#include "createOutput.H"
#include "AuxiliaryFunctions.h"
#include "AuxiliaryFunctions.hpp"


int main(int argc, char** argv)
{
    std::cout.setf(std::ios::scientific);

    boost::filesystem::path executable_file = boost::filesystem::system_complete(argv[0]);
    boost::filesystem::path executable_folder = executable_file.parent_path();

    OpenSMOKE::OpenSMOKE_logo("bioSMOKE 2024", "Riccardo Caraccio (riccardo.caraccio@polimi.it)");

    unsigned int max_number_allowed_species = 100000;
    //OpenSMOKE::OpenSMOKE_CheckLicense(executable_folder, "bioSMOKE", max_number_allowed_species);

    std::string input_file_name_ = "input.dic";
    std::string main_dictionary_name_ = "Biomass";

    // Program options from command line
    {
        namespace po = boost::program_options;
        po::options_description description("Options for the OpenSMOKE_Biomass");
        description.add_options()
            ("help", "print help messages")
            ("input", po::value<std::string>(), "name of the file containing the main dictionary (default \"input.dic\")")
            ("dictionary", po::value<std::string>(), "name of the main dictionary to be used (default \"Biomass\")");

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
        catch (po::error& e)
        {
            std::cerr << "Fatal error: " << e.what() << std::endl << std::endl;
            std::cerr << description << std::endl;
            return OPENSMOKE_FATAL_ERROR_EXIT;
        }
    }

    ///////////////////////////////////////////////////////////////////////
    //                       READING INPUT FILE                          //
    ///////////////////////////////////////////////////////////////////////

    // Defines the grammar rules
    OpenSMOKE::Grammar_Biomass grammar_biomass;
    OpenSMOKE::BiomassTemperature_Profile* vector_profile = NULL;

    // Define the dictionaries
    OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
    dictionaries.ReadDictionariesFromFile(input_file_name_);
    dictionaries(main_dictionary_name_).SetGrammar(grammar_biomass);

    // Kinetic scheme
    boost::filesystem::path path_kinetics_output;
    if (dictionaries(main_dictionary_name_).CheckOption("@KineticsFolder") == true)
        dictionaries(main_dictionary_name_).ReadPath("@KineticsFolder", path_kinetics_output);


    // Read the kinetic scheme in XML format
    {
        {
            if (!boost::filesystem::exists(path_kinetics_output / "kinetics.xml"))
                OpenSMOKE::FatalErrorMessage("The kinetic mechanism does not exist! Please check");

            boost::property_tree::ptree ptree;
            boost::property_tree::read_xml((path_kinetics_output / "kinetics.xml").string(), ptree);

            double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
            thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree);
            kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML, ptree);
            transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(ptree);
            double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
            std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;

            if (thermodynamicsMapXML->NumberOfSpecies() > max_number_allowed_species)
            {
                std::stringstream tag; tag << max_number_allowed_species;
                OpenSMOKE::FatalErrorMessage("The OpenSMOKE++ license you are using is limited to " + tag.str() + " species");
            }
        }

        // Read the solid-phase reactions
        {
            std::cout << "Reading the kinetic scheme (solid phase) in XML format" << std::endl;

            if (!boost::filesystem::exists(path_kinetics_output / "kinetics.solid.xml"))
            {
                std::cout << "No reactions in the solid phase" << std::endl;
            }
            else
            {
                boost::property_tree::ptree ptree;
                boost::property_tree::read_xml((path_kinetics_output / "kinetics.solid.xml").string(), ptree);

                double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

                thermodynamicsSolidMapXML = new OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN(ptree);
                kineticsSolidMapXML.resize(thermodynamicsSolidMapXML->number_of_materials());
                for (unsigned int k = 0; k < thermodynamicsSolidMapXML->number_of_materials(); k++)
                    kineticsSolidMapXML[k] = new OpenSMOKE::KineticsMap_Solid_CHEMKIN(*thermodynamicsSolidMapXML, ptree, k + 1);
                //kineticsSolidMapXML = new OpenSMOKE::KineticsMap_Solid_CHEMKIN(*thermodynamicsSolidMapXML, ptree, 1);

                double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
                std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;
            }
        }
    }

    std::string typeAnalysis;
    if (dictionaries(main_dictionary_name_).CheckOption("@Analysis") == true)
    {

        dictionaries(main_dictionary_name_).ReadString("@Analysis", typeAnalysis);

        if (typeAnalysis == "TG_Analysis")
        {
            std::string name_of_solid_status_subdictionary;
            if (dictionaries(main_dictionary_name_).CheckOption("@TGA"), true)
                dictionaries(main_dictionary_name_).ReadDictionary("@TGA", name_of_solid_status_subdictionary);

            Get_TGAanalysisFromDictionary(dictionaries(name_of_solid_status_subdictionary), analysis, Heat_Rates, final_time, output_species_);

            if (output_species_[0] == "all")
            {
                output_species_ = kineticsSolidMapXML[0]->NamesOfSpecies();
            }
        }
        else if (typeAnalysis == "Total_Analysis")
        {
            std::string name_of_solid_status_subdictionary;
            if (dictionaries(main_dictionary_name_).CheckOption("@Total_Analysis") == true)
                dictionaries(main_dictionary_name_).ReadDictionary("@Total_Analysis", name_of_solid_status_subdictionary);

            Get_TotalSimulation_analysisFromDictionary(dictionaries(name_of_solid_status_subdictionary), analysis, energyBalance,
                volumeLoss, final_time, epsi, tau, intervalli, raggioTot, Da, output_species_);
        }

        std::string name_of_gas_status_subdictionary;
        if (dictionaries(main_dictionary_name_).CheckOption("@TemperatureProfile") == true)
        {
            dictionaries(main_dictionary_name_).ReadDictionary("@TemperatureProfile", name_of_gas_status_subdictionary);

            OpenSMOKE::OpenSMOKEVectorDouble x, y;
            std::string x_variable, y_variable;
            GetXYProfileFromDictionary(dictionaries(name_of_gas_status_subdictionary), x, y, x_variable, y_variable);

            temperature_profile = true;
            vector_profile = new OpenSMOKE::BiomassTemperature_Profile(x, y);
            temperature_profile_ = vector_profile;
            final_time = x[x.Size()];
        }

    }
    else OpenSMOKE::FatalErrorMessage("Unknown simulation type: " + typeAnalysis + "  Conditions:TG or Total_Analysis");

    // Read initial conditions
    {
        std::string name_of_gas_status_subdictionary;
        if (dictionaries(main_dictionary_name_).CheckOption("@InletGasStatus") == true)
            dictionaries(main_dictionary_name_).ReadDictionary("@InletGasStatus", name_of_gas_status_subdictionary);

        GetGasStatusFromDictionary(dictionaries(name_of_gas_status_subdictionary), *thermodynamicsMapXML, Tbulk, Pbulk, rhoGas, omegaIn_gas);
    }

    {
        std::string name_of_solid_status_subdictionary;
        if (dictionaries(main_dictionary_name_).CheckOption("@InletSolidStatus") == true)
            dictionaries(main_dictionary_name_).ReadDictionary("@InletSolidStatus", name_of_solid_status_subdictionary);

        GetSolidStatusFromDictionary(dictionaries(name_of_solid_status_subdictionary), *thermodynamicsSolidMapXML, Tsolid, Psolid,
            rhoSolid, omegaIn_solid, hext, lambdaS);

        if (temperature_profile == true)
            Tsolid = temperature_profile_->Get(0.);

    }

    //Calcolo il numero di specie solide, gas ed il numero di reazioni totali;

    std::cout << "Number of materials:           " << thermodynamicsSolidMapXML->number_of_materials() << std::endl;
    solid_species = thermodynamicsSolidMapXML->number_of_solid_species();
    std::cout << "Number of solid-phase species: " << solid_species << std::endl;
    gas_species = thermodynamicsSolidMapXML->number_of_gas_species();
    std::cout << "Number of gas-phase species: " << gas_species << std::endl;
    std::cout << "Number of intervals: " << intervalli << std::endl;
    Neq = (gas_species + solid_species + 1);
    std::cout << "Number of equations for each interval: " << Neq << std::endl;
    Ntot_eq = Neq * intervalli;
    std::cout << "Number of total equations: " << Ntot_eq << std::endl;

    for (unsigned int k = 0; k < thermodynamicsSolidMapXML->number_of_materials(); k++)
    {
        solid_reactions = kineticsSolidMapXML[k]->NumberOfReactions();
        std::cout << "Number of solid reactions: " << solid_reactions << std::endl;
    }
    gas_reactions = kineticsMapXML->NumberOfReactions();
    std::cout << "Number of gas reactions: " << gas_reactions << std::endl;


    std::cout << "*Density read from input file [kg/m3]\t"  << rhoSolid << std::endl; ;
    std::cout << "N.B It's not apparent density " << std::endl;
 
    ChangeDimensions(solid_species + gas_species, &MW_tot, true);
    ChangeDimensions(solid_species, &MW_solid, true);
    ChangeDimensions(intervalli, &P, true);
    ChangeDimensions(gas_species, &xGasIn, true);
    ChangeDimensions(gas_species, &cGas_TGA, true);
    ChangeDimensions(gas_species, &Kc, true);
    ChangeDimensions(gas_species, &Diff, true);

    // Calculate the properties for the entering gas
    thermodynamicsMapXML->SetPressure(Pbulk);       transportMapXML->SetPressure(Pbulk);
    thermodynamicsMapXML->SetTemperature(Tbulk);    transportMapXML->SetTemperature(Tbulk);

    MW_gas_in = thermodynamicsMapXML->MolecularWeight_From_MassFractions(omegaIn_gas.GetHandle());
    thermodynamicsMapXML->MoleFractions_From_MassFractions(xGasIn.GetHandle(), MW_gas_in, omegaIn_gas.GetHandle()); 
    ExtGas_cp = thermodynamicsMapXML->cpMolar_Mixture_From_MoleFractions(xGasIn.GetHandle()) / MW_gas_in;
    ExtGas_lambda = transportMapXML->ThermalConductivity(xGasIn.GetHandle());
    transportMapXML->MassDiffusionCoefficients(Diff.GetHandle(), xGasIn.GetHandle());

    // Mass transfer coefficents calculated trought the Colburn analogy
    for (int i = 1; i <= gas_species; i++) {
        Kc[i] = hext * pow((pow(Diff[i], 2) / (pow(ExtGas_lambda, 2) * rhoGas * ExtGas_cp)), 1.0 / 3.0);
    }

    for (int i = 1; i <= P.Size(); i++)
        P[i] = Psolid;

    double rhoGas_internal = Psolid * MW_gas_in / PhysicalConstants::R_J_kmol / Tsolid;

    // MW_tot = [MW_solid, MW_gas]
    for (unsigned int j = 0; j < solid_species; j++)
        MW_tot[j + 1] = thermodynamicsSolidMapXML->MW(gas_species + j);

    for (unsigned int j = 0; j < thermodynamicsMapXML->NumberOfSpecies(); j++)
        MW_tot[j + 1 + solid_species] = thermodynamicsMapXML->MW(j);

    for (unsigned int j = 0; j < solid_species; j++)
        MW_solid[j + 1] = thermodynamicsSolidMapXML->MW(gas_species + j);


    ChangeDimensions(solid_reactions, solid_species + gas_species, &nu, true);

    for (unsigned int k = 0; k < thermodynamicsSolidMapXML->number_of_materials(); k++)
    {
        kineticsSolidMapXML[k]->stoichiometry().BuildStoichiometricMatrix();

        for (int z = 0; z < kineticsSolidMapXML[k]->stoichiometry().stoichiometric_matrix_products().outerSize(); ++z)
        {
            if (z <= (gas_species - 1))
            {
                for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsSolidMapXML[k]->stoichiometry().stoichiometric_matrix_products(), z); it; ++it)
                {
                    nu[it.index() + 1][z + 1 + solid_species] = it.value();
                }
            }
            else
            {
                for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsSolidMapXML[k]->stoichiometry().stoichiometric_matrix_products(), z); it; ++it)
                {
                    nu[it.index() + 1][z + 1 - gas_species] = it.value();
                }

            }
        }

        for (int z = 0; z < kineticsSolidMapXML[k]->stoichiometry().stoichiometric_matrix_reactants().outerSize(); ++z)
        {
            if (z <= (gas_species - 1))
            {
                for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsSolidMapXML[k]->stoichiometry().stoichiometric_matrix_reactants(), z); it; ++it)
                {
                    nu[it.index() + 1][z + 1 + solid_species] = -it.value();
                }
            }
            else
            {
                for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsSolidMapXML[k]->stoichiometry().stoichiometric_matrix_reactants(), z); it; ++it)
                {
                    nu[it.index() + 1][z + 1 - gas_species] = -it.value();
                }

            }
        }
    }

    ////////////////////////////////////////////////////////////////////////
    //                    PREPARE OUTPUT FILEs                            //
    ////////////////////////////////////////////////////////////////////////
    //  
    // Mattia 2017 Conversione LINUX-->WIN

    #ifdef linux
        DIR* dp;
        dp = opendir("Output");
        if (dp == NULL)
        {
            const int dir_err = mkdir("Output", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            if (-1 == dir_err)
            {
                printf("Error creating directory!n");
                exit(1);
            }
        }
    #endif

    #ifdef _WIN32
        if (!boost::filesystem::exists("Output"))
            OpenSMOKE::CreateDirectory("Output");

        /*boost::filesystem::remove_all("Output");
        std::string mkdir_command = "mkdir ";
        std::string outputFolderName = "Output";
        std::string message2 = mkdir_command + outputFolderName;
        system(message2.c_str());*/
    #endif

    if (typeAnalysis == "TG_Analysis")
    {
        {
            TG.open("Output/TG.ris", std::ios::out);
            TG.setf(std::ios::scientific);

            unsigned int counter = 1;
            OpenSMOKE::PrintTagOnASCIILabel(25, TG, "t[s]", counter);
            OpenSMOKE::PrintTagOnASCIILabel(25, TG, "t[min]", counter);
            OpenSMOKE::PrintTagOnASCIILabel(25, TG, "T[°C]", counter);
            OpenSMOKE::PrintTagOnASCIILabel(25, TG, "T[K]", counter);
            OpenSMOKE::PrintTagOnASCIILabel(25, TG, "TG", counter);

            TG << std::endl;
        }

        {
            finalYield.open("Output/finalYield.ris", std::ios::out);
            finalYield.setf(std::ios::scientific);
        }

        {
            finalYieldSolid.open("Output/finalYieldSolid.ris", std::ios::out);
            finalYieldSolid.setf(std::ios::scientific);
        }

        {
            finalYieldGas.open("Output/finalYieldGas.ris", std::ios::out);
            finalYieldGas.setf(std::ios::scientific);
        }

        {
            yieldScheletal.open("Output/yieldSpeciesSelective.ris", std::ios::out);
            yieldScheletal.setf(std::ios::scientific);
        }

        {
            yieldSpecies.open("Output/yieldSpecies.ris", std::ios::out);
            yieldSpecies.setf(std::ios::scientific);

            unsigned int counter = 1;
            OpenSMOKE::PrintTagOnASCIILabel(25, yieldSpecies, "t[s]", counter);
            OpenSMOKE::PrintTagOnASCIILabel(25, yieldSpecies, "t[min]", counter);
            OpenSMOKE::PrintTagOnASCIILabel(25, yieldSpecies, "T[°C]", counter);
            OpenSMOKE::PrintTagOnASCIILabel(25, yieldSpecies, "T[K]", counter);

            unsigned int NS_ = thermodynamicsSolidMapXML->NumberOfSpecies();
            unsigned int NSolid_ = thermodynamicsSolidMapXML->number_of_solid_species();
            unsigned int NGas_ = thermodynamicsSolidMapXML->number_of_gas_species();

            std::vector<unsigned int> widths_of_formation_rates_species_;
            widths_of_formation_rates_species_.resize(NS_);

            for (unsigned int i = 0; i < NSolid_; i++)
                widths_of_formation_rates_species_[i] = OpenSMOKE::CalculateSpeciesFieldWidth(thermodynamicsSolidMapXML->NamesOfSpecies()[i + NGas_], NSolid_);

            for (unsigned int i = NSolid_; i < NS_; i++)
                widths_of_formation_rates_species_[i] = OpenSMOKE::CalculateSpeciesFieldWidth(thermodynamicsMapXML->NamesOfSpecies()[i - NSolid_], NGas_);

            for (int i = 0; i < output_species_.size(); i++)
            {
                OpenSMOKE::PrintTagOnASCIILabel(25, yieldSpecies, output_species_[i], counter);
            }

            yieldSpecies << std::endl;

        }
    }
    if (typeAnalysis == "Total_Analysis")
    {
        {
            temp.open("Output/temperature.ris", std::ios::out);
            temp.setf(std::ios::scientific);

            unsigned int counter = 1;
            OpenSMOKE::PrintTagOnASCIILabel(25, temp, "t[s]", counter);

            for (int i = 1; i <= intervalli; i++)
            {
                string shell = boost::lexical_cast<string>(i);
                OpenSMOKE::PrintTagOnASCIILabel(25, temp, "T[K]_shell_" + shell, counter);
            }

            temp << std::endl;
        }

        {
            press.open("Output/pressure.ris", std::ios::out);
            press.setf(std::ios::scientific);

            unsigned int counter = 1;
            OpenSMOKE::PrintTagOnASCIILabel(25, press, "t[s]", counter);

            for (int i = 1; i <= intervalli; i++)
            {
                string shell = boost::lexical_cast<string>(i);
                OpenSMOKE::PrintTagOnASCIILabel(25, press, "P[Pa]_shell_" + shell, counter);
            }

            press << std::endl;
        }

        {
            fmass.open("Output/mass_shells.ris", std::ios::out);
            fmass.setf(std::ios::scientific);

            unsigned int counter = 1;
            OpenSMOKE::PrintTagOnASCIILabel(25, fmass, "t[s]", counter);

            for (int i = 1; i <= intervalli; i++)
            {
                string shell = boost::lexical_cast<string>(i);
                OpenSMOKE::PrintTagOnASCIILabel(25, fmass, "m[kg]_shell_" + shell, counter);
            }

            fmass << std::endl;
        }

        {
            por.open("Output/porosity.ris", std::ios::out);
            por.setf(std::ios::scientific);

            unsigned int counter = 1;
            OpenSMOKE::PrintTagOnASCIILabel(25, por, "t[s]", counter);

            for (int i = 1; i <= intervalli; i++)
            {
                string shell = boost::lexical_cast<string>(i);
                OpenSMOKE::PrintTagOnASCIILabel(25, por, "epsi_shell_" + shell, counter);
            }

            por << std::endl;

        }

        {
            volume.open("Output/volume.ris", std::ios::out);
            volume.setf(std::ios::scientific);

            unsigned int counter = 1;
            OpenSMOKE::PrintTagOnASCIILabel(25, volume, "t[s]", counter);
            OpenSMOKE::PrintTagOnASCIILabel(25, volume, "V(t)/Vstart", counter);
            OpenSMOKE::PrintTagOnASCIILabel(25, volume, "D(t)/Dstart", counter);

            for (int i = 1; i <= intervalli; i++)
            {
                string shell = boost::lexical_cast<string>(i);
                OpenSMOKE::PrintTagOnASCIILabel(25, volume, "V[m3]_shell_" + shell, counter);
            }

            volume << std::endl;
        }

        {
            TG.open("Output/TG.ris", std::ios::out);
            TG.setf(std::ios::scientific);

            unsigned int counter = 1;
            OpenSMOKE::PrintTagOnASCIILabel(25, TG, "t[s]", counter);
            OpenSMOKE::PrintTagOnASCIILabel(25, TG, "m(t)/mStart", counter);

            TG << std::endl;

        }

        {
            yieldSpecies.open("Output/yieldSpeciesfromR.ris", std::ios::out);
            yieldSpecies.setf(std::ios::scientific);

            unsigned int counter = 1;
            OpenSMOKE::PrintTagOnASCIILabel(25, yieldSpecies, "t[s]", counter);

            unsigned int NS_ = thermodynamicsSolidMapXML->NumberOfSpecies();
            unsigned int NSolid_ = thermodynamicsSolidMapXML->number_of_solid_species();
            unsigned int NGas_ = thermodynamicsSolidMapXML->number_of_gas_species();

            std::vector<unsigned int> widths_of_formation_rates_species_;
            widths_of_formation_rates_species_.resize(NS_);

            for (unsigned int i = 0; i < NSolid_; i++)
                widths_of_formation_rates_species_[i] = OpenSMOKE::CalculateSpeciesFieldWidth(thermodynamicsSolidMapXML->NamesOfSpecies()[i + NGas_], NSolid_);

            for (unsigned int i = NSolid_; i < NS_; i++)
                widths_of_formation_rates_species_[i] = OpenSMOKE::CalculateSpeciesFieldWidth(thermodynamicsMapXML->NamesOfSpecies()[i - NSolid_], NGas_);

            OpenSMOKE::PrintTagOnASCIILabel(25, yieldSpecies, "total mass[-]", counter);

            for (int i = 0; i < NS_; i++)
            {
                if (i < NSolid_)
                    OpenSMOKE::PrintTagOnASCIILabel(25, yieldSpecies, thermodynamicsSolidMapXML->NamesOfSpecies()[i + NGas_] + "[kg/kg_s0]", counter);
                else
                    OpenSMOKE::PrintTagOnASCIILabel(25, yieldSpecies, thermodynamicsMapXML->NamesOfSpecies()[i - NSolid_] + "[kg/kg_s0]", counter);
            }

            yieldSpecies << std::endl;

        }

        {
            fmassTotal.open("Output/massSpeciesTotal.ris", std::ios::out);
            fmassTotal.setf(std::ios::scientific);

            unsigned int counter = 1;
            OpenSMOKE::PrintTagOnASCIILabel(25, fmassTotal, "t[s]", counter);

            unsigned int NS_ = thermodynamicsSolidMapXML->NumberOfSpecies();
            unsigned int NSolid_ = thermodynamicsSolidMapXML->number_of_solid_species();
            unsigned int NGas_ = thermodynamicsSolidMapXML->number_of_gas_species();

            std::vector<unsigned int> widths_of_formation_rates_species_;
            widths_of_formation_rates_species_.resize(NS_);

            for (unsigned int i = 0; i < NSolid_; i++)
                widths_of_formation_rates_species_[i] = OpenSMOKE::CalculateSpeciesFieldWidth(thermodynamicsSolidMapXML->NamesOfSpecies()[i + NGas_], NSolid_);

            for (unsigned int i = NSolid_; i < NS_; i++)
                widths_of_formation_rates_species_[i] = OpenSMOKE::CalculateSpeciesFieldWidth(thermodynamicsMapXML->NamesOfSpecies()[i - NSolid_], NGas_);

            OpenSMOKE::PrintTagOnASCIILabel(25, fmassTotal, "total mass[kg]", counter);

            for (int i = 0; i < NS_; i++)
            {
                if (i < NSolid_)
                    OpenSMOKE::PrintTagOnASCIILabel(25, fmassTotal, thermodynamicsSolidMapXML->NamesOfSpecies()[i + NGas_] + "[kg]", counter);
                else
                    OpenSMOKE::PrintTagOnASCIILabel(25, fmassTotal, thermodynamicsMapXML->NamesOfSpecies()[i - NSolid_] + "[kg]", counter);
            }

            fmassTotal << std::endl;

        }
        /*{

            const boost::filesystem::path file_name = ("Output/massTotal.Cyl");
            fmassTotal.open(file_name.c_str(), std::ios::out);
            fmassTotal.setf(std::ios::scientific);

            unsigned int counter = 1;
            OpenSMOKE::PrintTagOnASCIILabel(25, fmassTotal, "t[s]", counter);

            unsigned int NS_ = thermodynamicsMapXML->NumberOfSpecies();
            std::vector<unsigned int> widths_of_formation_rates_species_;
            widths_of_formation_rates_species_.resize(NS_);
            for (unsigned int i = 0; i < NS_; i++)
                widths_of_formation_rates_species_[i] = OpenSMOKE::CalculateSpeciesFieldWidth(thermodynamicsMapXML->NamesOfSpecies()[i], NS_);

            for (int i = 0; i < NS_; i++)
            {
                OpenSMOKE::PrintTagOnASCIILabel(25, fmassTotal, thermodynamicsMapXML->NamesOfSpecies()[i] + "[kg]", counter);
            }

            fmassTotal << std::endl;

        }*/

        {

            const boost::filesystem::path file_name = ("Output/gasFlowRate.Cyl");
            gasFlowRate.open(file_name.c_str(), std::ios::out);
            gasFlowRate.setf(std::ios::scientific);

            unsigned int counter = 1;
            OpenSMOKE::PrintTagOnASCIILabel(25, gasFlowRate, "t[s]", counter);

            unsigned int NS_ = thermodynamicsMapXML->NumberOfSpecies();
            std::vector<unsigned int> widths_of_formation_rates_species_;
            widths_of_formation_rates_species_.resize(NS_);
            for (unsigned int i = 0; i < NS_; i++)
                widths_of_formation_rates_species_[i] = OpenSMOKE::CalculateSpeciesFieldWidth(thermodynamicsMapXML->NamesOfSpecies()[i], NS_);

            for (int i = 0; i < NS_; i++)
            {
                OpenSMOKE::PrintTagOnASCIILabel(25, gasFlowRate, thermodynamicsMapXML->NamesOfSpecies()[i] + "[kg/s]", counter);
            }

            gasFlowRate << std::endl;

        }

        {
            const boost::filesystem::path file_name = ("Output/gasFlowRateSelective.Cyl");
            gasFlowRateSelective.open(file_name.c_str(), std::ios::out);
            gasFlowRateSelective.setf(std::ios::scientific);

            unsigned int counter = 1;
            OpenSMOKE::PrintTagOnASCIILabel(25, gasFlowRateSelective, "t[s]", counter);
            OpenSMOKE::PrintTagOnASCIILabel(25, gasFlowRateSelective, "LightGas[kg/s]", counter);
            OpenSMOKE::PrintTagOnASCIILabel(25, gasFlowRateSelective, "Tar[kg/s]", counter);


            gasFlowRateSelective << std::endl;
        }

        {
            species.resize(output_species_.size());
            for (int i = 0; i < output_species_.size(); i++)
            {

                const boost::filesystem::path file_name = ("Output/" + output_species_[i] + ".Cyl");
                species[i].open(file_name.c_str(), std::ios::out);
                species[i].setf(std::ios::scientific);

                unsigned int counter = 1;
                OpenSMOKE::PrintTagOnASCIILabel(25, species[i], "t[s]", counter);
                //OpenSMOKE::PrintTagOnASCIILabel(25, temp, "t[min]", counter);

                for (int j = 1; j <= intervalli; j++)
                {
                    string shell = boost::lexical_cast<string>(j);
                    OpenSMOKE::PrintTagOnASCIILabel(25, species[i], "xShell_" + shell + "[wt]", counter);
                }

                species[i] << std::endl;

            }

        }

        {
            QreactTot.open("Output/QreactTot.ris", std::ios::out);
            QreactTot.setf(std::ios::scientific);

            unsigned int counter = 1;
            OpenSMOKE::PrintTagOnASCIILabel(25, QreactTot, "t[s]", counter);
            OpenSMOKE::PrintTagOnASCIILabel(25, QreactTot, "QreactTot[W]", counter);


            QreactTot << std::endl;
        }

        {
            QreactShell.open("Output/QreactShell.ris", std::ios::out);
            QreactShell.setf(std::ios::scientific);

            unsigned int counter = 1;
            OpenSMOKE::PrintTagOnASCIILabel(25, QreactShell, "t[s]", counter);


            for (int i = 1; i <= intervalli; i++)
            {
                string shell = boost::lexical_cast<string>(i);
                OpenSMOKE::PrintTagOnASCIILabel(25, QreactShell, "Qreact_shell_" + shell + "[W]", counter);
            }

            QreactShell << std::endl;

        }
    }


    ////////////////////////////////////////////////////////////////////////
    //                           INTEGRATION                              //
    ////////////////////////////////////////////////////////////////////////

    ChangeDimensions(intervalli, &epsi_var, true);
    ChangeDimensions(intervalli, &massTotsolid_initial_parz, true);

    if (typeAnalysis == "TG_Analysis")
    {
        // Costruisco i calcoli imponendo volume unitario // x0 = [mi_solid, mi_gas, T]
        Eigen::VectorXd x0(Neq), x;

        for (int j = 0; j < Neq; j++)
        {
            if (j < solid_species)
                x0[j] = omegaIn_solid[j + 1] * rhoSolid *1;
            else if (j >= solid_species && j < solid_species + gas_species)
                x0[j] = omegaIn_gas[j + 1 - solid_species] * rhoGas_internal *1;
            else
                x0[j] = Tsolid;
        }

        // Compute initial mass
        massTotsolid_initial = 0.;
        for (int j = 0; j < (solid_species); j++)
            massTotsolid_initial += x0[j];

        // Create the solver
        typedef OdeSMOKE::KernelDense<OpenSMOKE::ODESystem_Interface> denseOde;
        typedef OdeSMOKE::MethodGear<denseOde> methodGear;
        OdeSMOKE::MultiValueSolver<methodGear> ode_solver;

        // Default ODE parameters
        OpenSMOKE::ODE_Parameters* ode_parameters_ = new OpenSMOKE::ODE_Parameters;
        ode_parameters_->SetAbsoluteTolerance(1e-14);
        ode_parameters_->SetRelativeTolerance(1e-9);

        // Set the ODE system of equations function and print function
        ode_solver.SetSystemOfEquations(&OdeTGA);
        ode_solver.SetPrintFunction(&MyStepPrintTGA);

        // Set initial conditions
        ode_solver.SetInitialConditions(0, x0);

        // Set linear algebra options
        ode_solver.SetLinearAlgebraSolver(ode_parameters_->linear_algebra());
        ode_solver.SetFullPivoting(ode_parameters_->full_pivoting());

        // Set relative and absolute tolerances
        ode_solver.SetAbsoluteTolerances(ode_parameters_->absolute_tolerance());
        ode_solver.SetRelativeTolerances(ode_parameters_->relative_tolerance());

        // Solve the system
        double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
        OdeSMOKE::OdeStatus status = ode_solver.Solve(final_time);
        double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();


        // Check the solution
        if (status > 0)
        {
            ode_solver.Solution(x);
        }


    }
    else if (typeAnalysis == "Total_Analysis")
    {
        // SET THE DIMENSIONS OF THE FOLLOWING VECTORS:
        // - raggio -> settato sul numero degli intervalli, definisce il raggio di ogni
        //             sfera concentrica
        // - r -> settato sul numero degli intervalli, definisce la coordinata di ogni punto
        //        della griglia imposto al centro di ogni guscio
        // - V -> settato sul numero degli intervalli, definisce il volume di ogni guscio 
        //        da non confondere con il volume di ogni sfera concentrica
        // - S -> settato sul numero degli intervalli, definisce la superficie esterna
        //        di ogni guscio che in questo caso coincide con la superficie di ogni sfera
        //        concentrica
        // N.B. Tutti i gusci sono settati a isomassa e isovolumeS

        ChangeDimensions(intervalli, &raggio, true);
        ChangeDimensions(intervalli, &Vrecord, true);
        ChangeDimensions(intervalli, &r, true);
        ChangeDimensions(intervalli, &V, true);
        ChangeDimensions(intervalli, &S, true);
        ChangeDimensions(intervalli, &Vstart, true);
        ChangeDimensions(intervalli, &Rstart, true);

        raggio[1] = raggioTot / pow(intervalli, 1. / 3.);

        if (intervalli > 1)
        {
            for (int i = 2; i <= intervalli; i++)
            {
                raggio[i] = pow(i, 1. / 3.) * raggio[1];
            }
        }

        r[1] = raggio[1] / 2.;
        if (intervalli > 1)
        {
            for (int k = 2; k <= intervalli; k++)
            {
                r[k] = (raggio[k] + raggio[k - 1]) / 2.;
            }
        }

        V[1] = 4. / 3. * PhysicalConstants::pi * pow(raggio[1], 3.);
        S[1] = 4. * PhysicalConstants::pi * pow(raggio[1], 2.);

        for (int k = 2; k <= intervalli; k++)
        {
            V[k] = 4. / 3. * PhysicalConstants::pi * (pow(raggio[k], 3.) - pow(raggio[k - 1], 3.));
            S[k] = 4. * PhysicalConstants::pi * pow(raggio[k], 2.);
        }

        for (int i = 1; i <= intervalli; i++)
        {
            Vstart[i] = V[i];
            Rstart[i] = raggio[i];
        }

        Vtot_iniz = V.SumElements();

        // Update value of epsi
        for (int i = 1; i <= intervalli; i++)
            epsi_var[i] = epsi;

        // Initial Conditions
        // Il vettore delle x0=m_0_j è così definito per l'i-esimo intervallo: 
        // n°specie solide + n°species gas + temperatura

        Eigen::VectorXd x0(Ntot_eq), x;

        for (int i = 1; i <= intervalli; i++)
        {
            for (int j = 0; j < Neq; j++)
            {
                if (j < solid_species)
                    x0[j + Neq * (i - 1)] = omegaIn_solid[j + 1] * rhoSolid * V[i] * (1 - epsi_var[i]);
                else if (j >= solid_species && j < solid_species + gas_species)
                    x0[j + Neq * (i - 1)] = omegaIn_gas[j + 1 - solid_species] * rhoGas_internal * V[i] * epsi_var[i];
                else
                    x0[j + Neq * (i - 1)] = Tsolid;
            }
        }

        // Compute initial mass
        massTotsolid_initial = 0.;
        for (int i = 1; i <= intervalli; i++)
        {
            for (int j = 0; j < (solid_species); j++)
                massTotsolid_initial += x0[j + Neq * (i - 1)];
        }

        massTotsolid_initial_parz = 0.;
        for (int i = 1; i <= intervalli; i++)
        {
            for (int j = 0; j < (solid_species); j++)
                massTotsolid_initial_parz[i] += x0[j + Neq * (i - 1)];
        }

        ChangeDimensions(gas_species, &massGasSpecies_outIntegrated, true);
        ChangeDimensions(gas_species, &massGasSpeciesFlowrate_old, true);
        ChangeDimensions(intervalli, gas_species + solid_species, &FormationRate_old, true);
        ChangeDimensions(gas_species, &massGasSpeciesTotal_fromR, true);
        ChangeDimensions(solid_species, &massSolidSpeciesTotal_fromR, true);

        massGasSpeciesFlowrate_old = 0.;
        massGasSpecies_outIntegrated = 0.;
        FormationRate_old = 0.;
        massGasSpeciesTotal_fromR = 0.;

        for (int i = 1; i <= intervalli; i++)
        {
            for (int j = 0; j < gas_species + solid_species; j++)
            {
                if (j < solid_species) 
                    massSolidSpeciesTotal_fromR[j + 1] += x0[j + Neq * (i - 1)];

                else
                    massGasSpeciesTotal_fromR[j + 1 - solid_species] += x0[j + Neq * (i - 1)];
			}
        }

     
        ChangeDimensionsFunction();
        // Create the solver
        typedef OdeSMOKE::KernelDense<OpenSMOKE::ODESystem_Interface> denseOde;
        typedef OdeSMOKE::MethodGear<denseOde> methodGear;
        OdeSMOKE::MultiValueSolver<methodGear> ode_solver;

        // Default ODE parameters
        OpenSMOKE::ODE_Parameters* ode_parameters_ = new OpenSMOKE::ODE_Parameters;
        ode_parameters_->SetAbsoluteTolerance(1e-14);
        ode_parameters_->SetRelativeTolerance(1e-9);


        // Set the ODE system of equations function and print function
        ode_solver.SetSystemOfEquations(&OdeTotal);
        ode_solver.SetPrintFunction(&MyStepPrintTotal);

        // Set initial conditions
        ode_solver.SetInitialConditions(0, x0);


        // Set linear algebra options
        ode_solver.SetLinearAlgebraSolver(ode_parameters_->linear_algebra());
        ode_solver.SetFullPivoting(ode_parameters_->full_pivoting());

        // Set relative and absolute tolerances
        ode_solver.SetAbsoluteTolerances(ode_parameters_->absolute_tolerance());
        ode_solver.SetRelativeTolerances(ode_parameters_->relative_tolerance());

        // Set minimum and maximum values
        Eigen::VectorXd yMin(Ntot_eq);
        Eigen::VectorXd yMax(Ntot_eq);

        for (int i = 1; i <= intervalli; i++)
        {
            for (int j = 0; j < Neq; j++)
            {
                if (j < (solid_species + gas_species)) 
                {
                    yMin[j + Neq * (i - 1)] = 0.;
                    yMax[j + Neq * (i - 1)] = 100.;
                }
                else
                {
                    yMin[j + Neq * (i - 1)] = 10.;
                    yMax[j + Neq * (i - 1)] = 10000.;
                }
                    
            }
        }

        ode_solver.SetMinimumValues(yMin);
        ode_solver.SetMaximumValues(yMax);
        
        
        // Solve the system
        double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
        OdeSMOKE::OdeStatus status = ode_solver.Solve(final_time);
        double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

        // Check the solution
        if (status > 0)
        {
            ode_solver.Solution(x);
        }

    }        
    else
        OpenSMOKE::FatalErrorMessage("Unknown simulation type");
    // _fcloseall;

    cout << "-----------------------------------------------------------------------------" << endl;
    cout << "//                               ALL DONE                                  //" << endl;
    cout << "-----------------------------------------------------------------------------" << endl;

    return 0;
}
