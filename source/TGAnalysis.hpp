#include "TGAnalysis.h"

#include "BioSMOKE_OdeInterfaces.h"
#include <math/native-ode-solvers/MultiValueSolver>
#include <math/OpenSMOKEFunctions.h>

namespace BioSMOKE
{
// clang-format off

TGAnalysis::TGAnalysis(OpenSMOKE::ThermodynamicsMap_CHEMKIN &thermodynamicsMap,
                       OpenSMOKE::KineticsMap_CHEMKIN &kineticsMap,
                       OpenSMOKE::TransportPropertiesMap_CHEMKIN &transportMap,
                       OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN &thermodynamicsSolidMap,
                       OpenSMOKE::KineticsMap_Solid_CHEMKIN &kineticsSolidMap,
                       OpenSMOKE::ODE_Parameters &ode_parameters,
                       BioSMOKE::BioSMOKE_Options &biosmoke_options,
                       const double T0,
                       const double P0,
                       const double rho0_solid,
                       const std::vector<double> &omega0_gas,
                       const std::vector<double> &omega0_solid,
                       const double heating_rate) :
            BaseSolver(thermodynamicsMap, kineticsMap, transportMap, thermodynamicsSolidMap,
                       kineticsSolidMap, ode_parameters, biosmoke_options)
{ // clang-format on

    iteration_ = 0;
    counter_file_video_ = 0;
    counter_file_ASCII_ = 0;

    V0_solid_ = 1;
    T0_solid_ = T0;
    P0_solid_ = P0;
    T0_gas_ = T0;
    P0_gas_ = P0;
    omega0_solid_ = omega0_solid;
    omega0_gas_ = omega0_gas;
    heating_rate_ = heating_rate;
    rho0_solid_ = rho0_solid;

    NGS_ = thermodynamicsMap_.NumberOfSpecies();
    NSS_ = thermodynamicsSolidMap_.number_of_solid_species();
    NC_ = NGS_ + NSS_;
    NE_ = NC_ + 1;

    y0_.resize(NE_);
    yf_.resize(NE_);
    omega_solid_.resize(NSS_);
    omega_gas_.resize(NGS_);
    x0_solid_.resize(NSS_);
    x_solid_.resize(NSS_);
    x0_gas_.resize(NGS_);
    x_gas_.resize(NGS_);
    mass_solid_.resize(NSS_);
    mass_gas_.resize(NGS_);

    thermodynamicsMap_.SetTemperature(T0_solid_);
    thermodynamicsMap_.SetPressure(P0_solid_);
    kineticsMap_.SetTemperature(T0_solid_);
    kineticsMap_.SetPressure(P0_solid_);

    thermodynamicsMap_.MoleFractions_From_MassFractions(x0_gas_.data(), MW0_gas_, omega0_gas_.data());
    rho0_gas_ = P0_gas_ * MW0_gas_ / (PhysicalConstants::R_J_kmol * T0_gas_);
    MW_gas_ = MW0_gas_;
    rho_gas_ = rho0_gas_;
    mass0_tot_gas_ = 0.;

    thermodynamicsSolidMap_.SetTemperature(T0_solid_);
    thermodynamicsSolidMap_.SetPressure(P0_solid_);
    kineticsSolidMap_.SetTemperature(T0_solid_);
    kineticsSolidMap_.SetPressure(P0_solid_);
    thermodynamicsSolidMap_.SolidMoleFractions_From_SolidMassFractions(x0_solid_.data(), MW0_solid_,
                                                                       omega0_solid_.data());
    MW_solid_ = MW0_solid_;
    V_solid_ = V0_solid_;
    rho_solid_ = rho0_solid;
    mass0_tot_solid_ = V0_solid_ * rho0_solid_;
    mass_tot_solid_ = mass0_tot_solid_;
    omega_solid_ = omega0_solid_;

    T_ = T0_solid_;
    P_ = P0_solid_;

    OpenAllFiles();
}

void TGAnalysis::PrepareASCIIFile(const boost::filesystem::path output_file_ascii)
{
    PrepareASCIIFile(fASCII_, output_file_ascii);
}

void TGAnalysis::PrepareASCIIFile(std::ofstream &fOutput, const boost::filesystem::path output_file_ascii)
{
    indices_of_output_species_.resize(biosmoke_options_.output_species().size());
    for (unsigned int i = 0; i < biosmoke_options_.output_species().size(); i++)
        indices_of_output_species_[i] =
            thermodynamicsSolidMap_.IndexOfSpecies(biosmoke_options_.output_species()[i]) - 1;

    if (indices_of_output_species_.size() != 0)
    {
        widths_of_output_species_.resize(biosmoke_options_.output_species().size());
        for (unsigned int i = 0; i < biosmoke_options_.output_species().size(); i++)
            widths_of_output_species_[i] =
                OpenSMOKE::CalculateSpeciesFieldWidth(biosmoke_options_.output_species()[i], NC_);
    }
    else
    {
        widths_of_output_species_.resize(NC_);
        for (unsigned int i = 0; i < NC_; i++)
            widths_of_output_species_[i] =
                OpenSMOKE::CalculateSpeciesFieldWidth(thermodynamicsSolidMap_.NamesOfSpecies()[i], NC_);
    }

    fOutput.open(output_file_ascii.c_str(), std::ios::out);

    unsigned int counter = 1;
    fOutput.setf(std::ios::scientific);
    OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "t[s]", counter);
    OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "T[K]", counter);
    OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "Ms/Ms0[-]", counter);
    OpenSMOKE::PrintTagOnASCIILabel(20, fOutput, "Mg/Ms0[-]", counter);

    if (indices_of_output_species_.size() != 0)
    {
        for (unsigned int i = 0; i < indices_of_output_species_.size(); i++)
            OpenSMOKE::PrintTagOnASCIILabel(
                widths_of_output_species_[i], fOutput,
                thermodynamicsSolidMap_.NamesOfSpecies()[indices_of_output_species_[i]] + "_M", counter);
        // for (unsigned int i = 0; i < indices_of_output_species_.size(); i++)
        //     OpenSMOKE::PrintTagOnASCIILabel(
        //         widths_of_output_species_[i], fOutput,
        //         thermodynamicsMap_.NamesOfSpecies()[indices_of_output_species_[i] - 1] + "_x", counter);
    }
    else
    {
        for (unsigned int i = 0; i < NC_; i++)
            OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput,
                                            thermodynamicsSolidMap_.NamesOfSpecies()[i] + "_M", counter);
        // for (unsigned int i = 0; i < NC_; i++)
        //     OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput,
        //                                     thermodynamicsSolidMap_.NamesOfSpecies()[i] + "_x", counter);
    }

    fOutput << std::endl;
}

void TGAnalysis::PrepareXMLFile(const boost::filesystem::path output_file_xml)
{
    fXML_.open(output_file_xml.c_str(), std::ios::out);
    OpenSMOKE::SetXMLFile(fXML_);
    fXML_ << "<Type> ThermogravimetricAnalysis </Type>" << std::endl;
    unsigned int counter = 2;
    fXML_ << "<additional>" << std::endl;
    fXML_ << 6 << std::endl;
    fXML_ << "time [s] " << counter++ << std::endl;
    fXML_ << "temperature [K] " << counter++ << std::endl;
    fXML_ << "pressure [Pa] " << counter++ << std::endl;
    fXML_ << "mol-weight [kg/kmol] " << counter++ << std::endl;
    fXML_ << "density [kg/m3] " << counter++ << std::endl;
    fXML_ << "heat-release [W/m3] " << counter++ << std::endl;
    fXML_ << "</additional>" << std::endl;

    fXML_ << "<t-p-mw>" << std::endl;
    fXML_ << 1 << " " << 2 << " " << 3 << std::endl;
    fXML_ << "</t-p-mw>" << std::endl;

    fXML_ << "<mass-fractions>" << std::endl;
    fXML_ << thermodynamicsSolidMap_.NumberOfSpecies() << std::endl;
    for (unsigned int j = 0; j < NC_; j++)
        fXML_ << thermodynamicsSolidMap_.NamesOfSpecies()[j] << " " << thermodynamicsSolidMap_.MW(j) << " " << counter++
              << std::endl;
    fXML_ << "</mass-fractions>" << std::endl;
    fXML_ << "<profiles>" << std::endl;
}

void TGAnalysis::CloseXMLFile()
{
    fXML_ << "</profiles>" << std::endl;
    fXML_ << "<profiles-size> " << std::endl;
    fXML_ << counter_file_XML_ << " " << 1 + (NC_ + 1) << std::endl;
    fXML_ << "</profiles-size> " << std::endl;

    // if (on_the_fly_post_processing_.is_active() == true)
    // {
    // 	fXML_ << "<formation-rates>" << std::endl;
    // 	fXML_ << "<!--units: kg/m3/s-->" << std::endl;
    // 	fXML_ << fXML_formation_rates_.str();
    // 	fXML_ << "</formation-rates>" << std::endl;

    // 	fXML_ << "<reaction-rates>" << std::endl;
    // 	fXML_ << "<!--units: kmol/m3/s-->" << std::endl;
    // 	fXML_ << fXML_reaction_rates_.str();
    // 	fXML_ << "</reaction-rates>" << std::endl;
    // }

    fXML_ << "</opensmoke>" << std::endl;
}

void TGAnalysis::OpenAllFiles()
{
    if (biosmoke_options_.verbose_output() == true)
    {
        if (!boost::filesystem::exists(biosmoke_options_.output_path()))
            OpenSMOKE::CreateDirectory(biosmoke_options_.output_path());

        if (biosmoke_options_.verbose_ascii_file() == true)
            PrepareASCIIFile(biosmoke_options_.output_path() / "Output.out");

        if (biosmoke_options_.verbose_xml_file() == true)
            PrepareXMLFile(biosmoke_options_.output_path() / "Output.xml");
    }
}

void TGAnalysis::PrintFinalStatus(std::ostream &fOutput, const double t)
{
    fOutput.setf(std::ios::scientific);
    fOutput << std::setw(20) << std::left << t;
    fOutput << std::setw(20) << std::left << T_;
    fOutput << std::setw(20) << std::left << mass_tot_solid_ / mass0_tot_solid_;
    fOutput << std::setw(20) << std::left << mass_tot_gas_ / mass0_tot_solid_;

    if (indices_of_output_species_.size() != 0)
    {
        for (unsigned int i = 0; i < indices_of_output_species_.size(); i++)
            if (indices_of_output_species_[i] < NGS_)
                fOutput << std::setw(widths_of_output_species_[i]) << std::left
                        << mass_gas_[indices_of_output_species_[i]] / (mass_tot_solid_ + mass_tot_gas_);
            else
                fOutput << std::setw(widths_of_output_species_[i]) << std::left
                        << mass_solid_[indices_of_output_species_[i] - NGS_] / (mass_tot_solid_ + mass_tot_gas_);
    }
    else
    {
        for (unsigned int i = 0; i < NC_; i++)
            if (i < NGS_)
                fOutput << std::setw(widths_of_output_species_[i]) << std::left
                        << mass_gas_[i] / (mass_tot_solid_ + mass_tot_gas_);
            else
                fOutput << std::setw(widths_of_output_species_[i]) << std::left
                        << mass_solid_[i - NGS_] / (mass_tot_solid_ + mass_tot_gas_);
    }
    fOutput << std::endl;
}

void TGAnalysis::CloseAllFiles()
{
    if (biosmoke_options_.verbose_output() == true)
    {
        if (biosmoke_options_.verbose_ascii_file() == true)
            fASCII_.close();

        if (biosmoke_options_.verbose_xml_file() == true)
            CloseXMLFile();
    }

    if (biosmoke_options_.sensitivity_analysis() == true)
        CloseSensitivityXMLFiles();
}

int TGAnalysis::Equations(const double t, const std::vector<double> &y, std::vector<double> &dy)
{
    // recover unknowns: mass_gas <> mass_solid <> T
    for (unsigned int i = 0; i < NE_; i++)
    {
        if (i < NGS_)
            mass_gas_[i] = y[i];
        else if (i < NGS_ + NSS_)
            mass_solid_[i - NGS_] = y[i];
        else
            T_ = y[i];
    }

    if (is_temperature_profile_ == true)
    {
        T_ = biosmoke_profile_->Get(t);
    }

    // set maps conditions
    thermodynamicsSolidMap_.SetTemperature(T_);
    thermodynamicsSolidMap_.SetPressure(P_);
    kineticsSolidMap_.SetTemperature(T_);
    kineticsSolidMap_.SetPressure(P_);

    // calculate total masses
    mass_tot_solid_ = std::accumulate(mass_solid_.begin(), mass_solid_.end(), 0.0);
    mass_tot_gas_ = std::accumulate(mass_gas_.begin(), mass_gas_.end(), 0.0);

    // caluclate solid mass fractions
    for (unsigned int i = 0; i < NSS_; i++)
        omega_solid_[i] = mass_solid_[i] / mass_tot_solid_;

    MW_solid_ = thermodynamicsSolidMap_.SolidMolecularWeight_From_SolidMassFractions(omega_solid_.data());

    // calculate gas mass fractions
    if (mass_tot_gas_ > 0)
        for (unsigned int i = 0; i < NGS_; i++)
            omega_gas_[i] = mass_gas_[i] / mass_tot_gas_;

    // calculate solid concentrations
    std::vector<double> cSolid_(NSS_, 0.);
    for (unsigned int i = 0; i < NSS_; i++)
        cSolid_[i] = rho_solid_ * omega_solid_[i] / thermodynamicsSolidMap_.MW(i + NGS_);

    // calculate gas concentrations
    double cTot_gas_ = P_ / (PhysicalConstants::R_J_kmol * T_);
    if (mass_tot_gas_ > 0)
        rho_gas_ = cTot_gas_ * thermodynamicsMap_.MolecularWeight_From_MassFractions(omega_gas_.data());
    else
        rho_gas_ = cTot_gas_ * thermodynamicsMap_.MolecularWeight_From_MassFractions(omega0_gas_.data());

    std::vector<double> cGas_(NGS_, 0.);
    for (unsigned int i = 0; i < NGS_; i++)
        cGas_[i] = cTot_gas_ * x0_gas_[i]; // we only use the inlet gas composition

    // calculate rates
    std::vector<double> R_gas_(NGS_, 0.);
    std::vector<double> R_solid_(NSS_, 0.);
    kineticsSolidMap_.ReactionRates(cGas_.data(), cSolid_.data());
    kineticsSolidMap_.FormationRates(R_gas_.data(), R_solid_.data());

    // calculate residuals
    for (unsigned int i = 0; i < NE_; i++)
    {
        if (i < NGS_)
            dy[i] = R_gas_[i] * thermodynamicsSolidMap_.MW(i) * (mass_tot_solid_ / rho_solid_);
        else if (i < NGS_ + NSS_)
            dy[i] = R_solid_[i - NGS_] * thermodynamicsSolidMap_.MW(i) * (mass_tot_solid_ / rho_solid_);
        else
            dy[i] = heating_rate_;
    }

    return 0;
}

void TGAnalysis::Solve(const double t0, const double tf)
{
    std::cout << std::endl;
    std::cout << "-----------------------------------------------------------------------------" << std::endl;
    std::cout << " Solving the TG analysis...                                                  " << std::endl;
    std::cout << "-----------------------------------------------------------------------------" << std::endl;

    for (unsigned int i = 0; i < NE_; i++)
    {
        if (i < NGS_)
            y0_[i] = 0.;
        else if (i >= NGS_ && i < NGS_ + NSS_)
            y0_[i] = omega0_solid_[i - NGS_] * rho0_solid_ * V0_solid_;
        else
            y0_[i] = T0_solid_;
    }

    // Set the final time, used to print the final status
    final_time_ = tf;

    // Print initial conditions
    {
        std::vector<double> dy0(y0_.size());
        Equations(t0, y0_, dy0);
        Print(t0, y0_);
    }

    // Min and max values
    Eigen::VectorXd yMin(NE_);
    for (unsigned int i = 0; i < NE_; i++)
        yMin(i) = 0.;
    Eigen::VectorXd yMax(NE_);
    for (unsigned int i = 0; i < NE_; i++)
        yMax(i) = 10000.;

    // Initial conditions
    const Eigen::VectorXd y0_eigen = Eigen::Map<Eigen::VectorXd>(y0_.data(), y0_.size());

    // Final solution
    Eigen::VectorXd yf_eigen(y0_eigen.size());

    typedef OdeSMOKE::KernelDense<BioSMOKE::ODESystem_BioSMOKE_TGAnalysis> denseOde;
    typedef OdeSMOKE::MethodGear<denseOde> methodGear;
    OdeSMOKE::MultiValueSolver<methodGear> ode_solver;
    ode_solver.SetReactor(this);

    // Set initial conditions
    ode_solver.SetInitialConditions(t0, y0_eigen);

    // Set linear algebra options
    ode_solver.SetLinearAlgebraSolver(ode_parameters_.linear_algebra());
    ode_solver.SetFullPivoting(ode_parameters_.full_pivoting());

    // Set relative and absolute tolerances
    ode_solver.SetAbsoluteTolerances(ode_parameters_.absolute_tolerance());
    ode_solver.SetRelativeTolerances(ode_parameters_.relative_tolerance());

    // Set minimum and maximum values
    ode_solver.SetMinimumValues(yMin);
    ode_solver.SetMaximumValues(yMax);

    // Set user defined Jacobian
    if (ode_parameters_.analytical_jacobian() == true)
    {
        ode_solver.SetUserDefinedJacobian();
        kineticsMap_.jacobian_sparsity_pattern_map()->SetEpsilon(ode_parameters_.absolute_tolerance() / 5.);
    }

    // Set maximum number of steps
    if (ode_parameters_.maximum_number_of_steps() > 0)
        ode_solver.SetMaximumNumberOfSteps(ode_parameters_.maximum_number_of_steps());

    // Set maximum integration order
    if (ode_parameters_.maximum_order() > 0)
        ode_solver.SetMaximumOrder(ode_parameters_.maximum_order());

    // Set maximum step size allowed
    if (ode_parameters_.maximum_step() > 0)
        ode_solver.SetMaximumStepSize(ode_parameters_.maximum_step());

    // Set minimum step size allowed
    if (ode_parameters_.minimum_step() > 0)
        ode_solver.SetMinimumStepSize(ode_parameters_.minimum_step());

    // Set initial step size
    if (ode_parameters_.initial_step() > 0)
        ode_solver.SetFirstStepSize(ode_parameters_.initial_step());

    // Solve the system
    double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
    OdeSMOKE::OdeStatus status = ode_solver.Solve(tf);
    double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

    // Check the solution
    if (status > 0)
    {
        ode_solver.Solution(yf_eigen);
        yf_ = std::vector<double>(yf_eigen.data(), yf_eigen.data() + yf_eigen.size());
        ode_parameters_.TransferDataFromOdeSolver(ode_solver, tEnd - tStart);
    }

    std::cout << std::endl;
    std::cout << "-----------------------------------------------------------------------------" << std::endl;
    std::cout << " Completed the simulation in " << std::setprecision(6) << tEnd - tStart << " seconds" << std::endl;
    std::cout << "-----------------------------------------------------------------------------" << std::endl;

    CloseAllFiles();
}

int TGAnalysis::Print(const double t, const std::vector<double> &y)
{
    iteration_++;

    if (biosmoke_options_.verbose_video() == true)
    {
        // Video output
        if (iteration_ % biosmoke_options_.n_step_video() == 1 || biosmoke_options_.n_step_video() == 1 ||
            t == final_time_)
        {
            counter_file_video_++;
            if (counter_file_video_ % 100 == 1)
            {
                std::cout << std::endl;
                std::cout << std::setw(10) << std::left << "#Step";
                std::cout << std::setw(16) << std::left << "Time[s]";
                std::cout << std::setw(10) << std::left << "T[K]";
                std::cout << std::setw(10) << std::left << "Ms/Ms0[-]";
                std::cout << std::endl;
            }
            std::cout << std::setw(10) << std::left << iteration_;
            std::cout << std::scientific << std::setw(16) << std::setprecision(6) << std::left << t;
            std::cout << std::setw(10) << std::left << std::fixed << std::setprecision(3) << T_;
            std::cout << std::setw(10) << std::left << std::fixed << std::setprecision(6)
                      << mass_tot_solid_ / mass0_tot_solid_;
            std::cout << std::endl;
        }

        if (biosmoke_options_.verbose_output() == true)
        {
            // ASCII output
            if (biosmoke_options_.verbose_ascii_file() == true)
            {
                if (iteration_ % biosmoke_options_.n_step_file() == 1 || biosmoke_options_.n_step_file() == 1 ||
                    t == final_time_)
                {
                    counter_file_ASCII_++;
                    PrintFinalStatus(fASCII_, t);
                }
            }

            // XML file output
            if (biosmoke_options_.verbose_xml_file() == true)
            {
                if (iteration_ % biosmoke_options_.n_step_file() == 1 || biosmoke_options_.n_step_file() == 1 ||
                    t == final_time_)
                {
                    counter_file_XML_++;
                    fXML_ << t << " ";
                    fXML_ << T_ << " ";
                    fXML_ << P_ << " ";
                    fXML_ << MW_solid_ << " ";
                    fXML_ << rho_solid_ << " ";
                    fXML_ << 0 << " "; // Qr
                    for (unsigned int i = 0; i < NC_; i++)
                        if (i < NGS_)
                            fXML_ << std::setprecision(12) << omega_gas_[i] << " ";
                        else
                            fXML_ << std::setprecision(12) << omega_solid_[i - NGS_] << " ";

                    fXML_ << std::endl;

                    // Write formation rates and reaction rates
                    // if (on_the_fly_post_processing_.is_active() == true)
                    // {
                    //     // Write formation rates (kg/m3/s)
                    //     for (unsigned int j = 1; j <= thermodynamicsMap_.NumberOfSpecies(); j++)
                    //         fXML_formation_rates_ << std::scientific << std::setprecision(9)
                    //                               << thermodynamicsMap_.MW(j - 1) * R_[j] << " ";
                    //     fXML_formation_rates_ << std::endl;

                    //     // Write reaction rates (kmol/m3/s)
                    //     std::vector<double> r = kineticsMap_.GiveMeReactionRates();
                    //     for (unsigned int j = 0; j < r.size(); j++)
                    //         fXML_reaction_rates_ << std::scientific << std::setprecision(9) << r[j] << " ";
                    //     fXML_reaction_rates_ << std::endl;
                    // }
                }
            }
        }
    }

    if (biosmoke_options_.sensitivity_analysis() == true)
        SensitivityAnalysis(t, y);

    return 0;
}

void TGAnalysis::SparseAnalyticalJacobian(const double t, const std::vector<double> &y, Eigen::SparseMatrix<double> &J)
{
    OpenSMOKE::ErrorMessage("TGAnalysis", "SparseAnalyticalJacobian is not yet available for TGAnalysis");
}

void TGAnalysis::DenseAnalyticalJacobian(const double t, const std::vector<double> &y, Eigen::MatrixXd &J)
{
    OpenSMOKE::ErrorMessage("TGAnalysis", "DenseAnalyticalJacobian is not yet available for TGAnalysis");
}

void TGAnalysis::NumericalJacobian(const double t, const std::vector<double> &y, OpenSMOKE::OpenSMOKEMatrixDouble &J)
{
    // Calculated as suggested by Buzzi (private communication)

    const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
    const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
    const double TOLR = 100. * OpenSMOKE::OPENSMOKE_MACH_EPS_FLOAT;
    const double TOLA = 1.e-10;

    std::vector<double> y_plus = y;
    std::vector<double> dy_original(y.size());
    std::vector<double> dy_plus(y.size());

    Equations(t, y, dy_original);

    // Derivatives with respect to y[kd]
    for (unsigned int kd = 0; kd < y.size(); kd++)
    {
        double hf = 1.e0;
        double error_weight = 1. / (TOLA + TOLR * std::fabs(y[kd]));
        double hJ = ETA2 * std::fabs(std::max(y[kd], 1. / error_weight));
        double hJf = hf / error_weight;
        hJ = std::max(hJ, hJf);
        hJ = std::max(hJ, ZERO_DER);

        // This is what is done by Buzzi
        double dy = std::min(hJ, 1.e-3 + 1e-3 * std::fabs(y[kd]));
        double udy = 1. / dy;
        y_plus[kd] += dy;
        Equations(t, y_plus, dy_plus);

        for (int j = 0; j < y.size(); j++)
            J[j + 1][kd + 1] = (dy_plus[j] - dy_original[j]) * udy;

        y_plus[kd] = y[kd];
    }
}

void TGAnalysis::PrepareSensitivityXMLFiles(OpenSMOKE::SensitivityAnalysis_Options &sensitivity_options)
{
    indices_of_sensitivity_species_.resize(sensitivity_options.list_of_species().size());
    for (unsigned int i = 0; i < indices_of_sensitivity_species_.size(); i++)
        indices_of_sensitivity_species_[i] =
            thermodynamicsSolidMap_.IndexOfSpecies(sensitivity_options.list_of_species()[i]);

    const boost::filesystem::path parent_file = biosmoke_options_.output_path() / "Sensitivities.xml";
    fSensitivityParentXML_.open(parent_file.c_str(), std::ios::out);
    OpenSMOKE::SetXMLFile(fSensitivityParentXML_);
    fSensitivityParentXML_ << "<variables>" << std::endl;

    fSensitivityParentXML_ << indices_of_sensitivity_species_.size() << std::endl;
    for (unsigned int j = 0; j < indices_of_sensitivity_species_.size(); j++)
        fSensitivityParentXML_ << thermodynamicsSolidMap_.NamesOfSpecies()[indices_of_sensitivity_species_[j] - 1]
                               << " " << j << " " << indices_of_sensitivity_species_[j] + 1 << std::endl;

    fSensitivityParentXML_ << "</variables>" << std::endl;
    fSensitivityParentXML_ << "<n-parameters> " << std::endl;
    fSensitivityParentXML_ << sensitivityMap_->number_of_parameters() << std::endl;
    fSensitivityParentXML_ << "</n-parameters> " << std::endl;

    fSensitivityChildXML_ = new std::ofstream[indices_of_sensitivity_species_.size()];
    for (unsigned int j = 0; j < indices_of_sensitivity_species_.size(); j++)
    {
        const std::string name = "Sensitivities." +
                                 thermodynamicsSolidMap_.NamesOfSpecies()[indices_of_sensitivity_species_[j] - 1] +
                                 ".xml";
        const boost::filesystem::path child_file = biosmoke_options_.output_path() / name;
        fSensitivityChildXML_[j].open(child_file.c_str(), std::ios::out);
    }

    // {
    //     const boost::filesystem::path child_file = biosmoke_options_.output_path() / "Sensitivities.temperature.xml";
    //     fSensitivityChildXML_[indices_of_sensitivity_species_.size()].open(child_file.c_str(), std::ios::out);
    // }

    for (unsigned int j = 0; j < indices_of_sensitivity_species_.size(); j++)
    {
        OpenSMOKE::SetXMLFile(fSensitivityChildXML_[j]);
        fSensitivityChildXML_[j] << std::setprecision(5);
        fSensitivityChildXML_[j] << "<coefficients>" << std::endl;
    }
}

void TGAnalysis::CloseSensitivityXMLFiles()
{
    fSensitivityParentXML_ << "<points> " << std::endl;
    fSensitivityParentXML_ << counter_sensitivity_XML_ << std::endl;
    fSensitivityParentXML_ << "</points> " << std::endl;
    fSensitivityParentXML_ << "<constant-parameters> " << std::endl;
    for (unsigned int j = 1; j <= sensitivityMap_->number_of_parameters(); j++)
        fSensitivityParentXML_ << sensitivityMap_->parameters()[j] << std::endl;
    fSensitivityParentXML_ << "</constant-parameters> " << std::endl;
    fSensitivityParentXML_ << "</opensmoke>" << std::endl;

    for (unsigned int j = 0; j < indices_of_sensitivity_species_.size(); j++)
    {
        fSensitivityChildXML_[j] << "</coefficients>" << std::endl;
        fSensitivityChildXML_[j] << "</opensmoke>" << std::endl;
    }
}

void TGAnalysis::EnableSensitivityAnalysis(OpenSMOKE::SensitivitySolidMap &sensitivityMap,
                                           OpenSMOKE::SensitivityAnalysis_Options &sensitivity_options)
{
    sensitivityMap_ = &sensitivityMap;

    PrepareSensitivityXMLFiles(sensitivity_options);

    OpenSMOKE::ChangeDimensions(NE_, &scaling_Jp_, true);

    if (sensitivityMap_->dense_solver_type() != OpenSMOKE::SOLVER_DENSE_NONE)
        OpenSMOKE::ChangeDimensions(NE_, NE_, &Jnum_, true);
    else
        Jan_.resize(NE_, NE_);
}

void TGAnalysis::SensitivityAnalysis(const double t, const std::vector<double> &y)
{
    if (iteration_ == 1)
    {
        // Writes the coefficients on file (only on request)
        if (iteration_ % biosmoke_options_.n_step_file() == 1 || biosmoke_options_.n_step_file() == 1)
        {
            counter_sensitivity_XML_++;
            for (unsigned int k = 0; k < indices_of_sensitivity_species_.size(); k++)
            {
                for (unsigned int j = 1; j <= sensitivityMap_->number_of_parameters(); j++)
                    fSensitivityChildXML_[k] << 0. << " ";
                fSensitivityChildXML_[k] << std::endl;
            }
        }
    }
    else
    {
        // Scaling factors
        for (unsigned int j = 1; j <= NC_; j++)
            if (j < NGS_)
                scaling_Jp_[j] = thermodynamicsSolidMap_.MW(j - 1) / rho_gas_;
            else
                scaling_Jp_[j] = thermodynamicsSolidMap_.MW(j - 1) / rho_solid_;

        // Calculates the current Jacobian
        if (sensitivityMap_->dense_solver_type() != OpenSMOKE::SOLVER_DENSE_NONE)
            NumericalJacobian(t, y, Jnum_);
        else
            SparseAnalyticalJacobian(t, y, Jan_);

        //  Recover concentrations
        OpenSMOKE::OpenSMOKEVectorDouble c_(NC_);
        const double cTot_gas_ = P_ / (PhysicalConstants::R_J_kmol * T_);
        for (unsigned int j = 0; j < NGS_; j++)
            c_[j + 1] = cTot_gas_ * x0_gas_[j]; // we only use the inlet gas composition

        for (unsigned int j = 0; j < NSS_; j++)
            c_[j + NGS_ + 1] = rho_solid_ * omega_solid_[j] / thermodynamicsSolidMap_.MW(j + NGS_);

        if (sensitivityMap_->dense_solver_type() != OpenSMOKE::SOLVER_DENSE_NONE)
            sensitivityMap_->CalculateSolid(t, T_, P_, c_, Jnum_, scaling_Jp_);
        else
            OpenSMOKE::FatalErrorMessage("Not yet implemented for TGAnalysis");
        // sensitivityMap_->CalculateSolid(t, T_, P_, c_, Jan_, scaling_Jp_);

        thermodynamicsSolidMap_.SolidMoleFractions_From_SolidMassFractions(x_solid_.data(), MW_solid_,
                                                                           omega_solid_.data());

        // Write the coefficentes on file (only on request)
        if (iteration_ % biosmoke_options_.n_step_file() == 1 || biosmoke_options_.n_step_file() == 1)
        {
            counter_sensitivity_XML_++;
            for (unsigned int k = 0; k < indices_of_sensitivity_species_.size(); k++)
            {
                const unsigned int i = indices_of_sensitivity_species_[k];
                if (i < NGS_)
                {
                    for (unsigned int j = 1; j <= sensitivityMap_->number_of_parameters(); j++)
                    {
                        double sum = 0.;
                        for (unsigned int kk = 1; kk <= NGS_; kk++)
                            sum += sensitivityMap_->sensitivity_coefficients()(kk - 1, j - 1) /
                                   thermodynamicsSolidMap_.MW(kk - 1);
                        sum *= x_gas_[i - 1] * MW_gas_;

                        double coefficient = sensitivityMap_->sensitivity_coefficients()(i - 1, j - 1) * MW_gas_ /
                                                 thermodynamicsSolidMap_.MW(i - 1) -
                                             sum;
                        fSensitivityChildXML_[k] << coefficient << " ";
                    }
                    fSensitivityChildXML_[k] << std::endl;
                }
                else // solid species
                {
                    for (unsigned int j = 0; j < sensitivityMap_->number_of_parameters(); j++)
                    {
                        double sum = 0.;
                        for (unsigned int kk = 0; kk < NSS_; kk++)
                            sum += sensitivityMap_->sensitivity_coefficients()(kk + NGS_, j) /
                                   thermodynamicsSolidMap_.MW(kk + NGS_);
                        sum *= x_solid_[i - NGS_] * MW_solid_;

                        double coefficient = sensitivityMap_->sensitivity_coefficients()(i, j) * MW_solid_ /
                                                 thermodynamicsSolidMap_.MW(i) -
                                             sum;
                        fSensitivityChildXML_[k] << coefficient << " ";
                    }
                    fSensitivityChildXML_[k] << std::endl;
                }
            }
        }
    }
}

} // namespace BioSMOKE