#include "TGAnalysis.h"

#include "BioSMOKE_OdeInterfaces.h"
#include <math/native-ode-solvers/MultiValueSolver>

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
    x0_gas_.resize(NGS_);

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

    T_ = T0_solid_;
    P_ = P0_solid_;
}

int TGAnalysis::Equations(const double t, const std::vector<double> &y, std::vector<double> &dy)
{
    // recover unknowns: mass_gas <> mass_solid <> T
    std::vector<double> mass_gas_current_(NGS_, 0.);
    std::vector<double> mass_solid_current_(NSS_, 0.);
    for (unsigned int i = 0; i < NE_; i++)
    {
        if (i < NGS_)
            mass_gas_current_[i] = y[i];
        else if (i < NGS_ + NSS_)
            mass_solid_current_[i - NGS_] = y[i];
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
    mass_tot_solid_ = std::accumulate(mass_solid_current_.begin(), mass_solid_current_.end(), 0.0);
    mass_tot_gas_ = std::accumulate(mass_gas_current_.begin(), mass_gas_current_.end(), 0.0);

    // caluclate solid mass fractions
    for (unsigned i = 0; i < NSS_; i++)
        omega_solid_[i] = mass_solid_current_[i] / mass_tot_solid_;

    // calculate solid concentrations
    std::vector<double> cSolid_(NSS_, 0.);
    for (unsigned int i = 0; i < NSS_; i++)
        cSolid_[i] = rho_solid_ * omega_solid_[i] / thermodynamicsSolidMap_.MW(i + NGS_);

    // calculate gas concentrations
    double cTot_gas_ = P_ / (PhysicalConstants::R_J_kmol * T_);
    std::vector<double> cGas_(NGS_, 0.);
    for (unsigned int i = 0; i < NGS_; i++)
        cGas_[i] = cTot_gas_ * x0_gas_[i];

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
            dy[i] = R_solid_[i - NGS_] * thermodynamicsSolidMap_.MW(i) * (mass_tot_solid_ / rho_gas_);
        else
            dy[i] = heating_rate_;
    }

    return 0;
}

void TGAnalysis::Solve(const double t0, const double tf)
{
    if (biosmoke_options_.verbose_video() == true)
    {
        std::cout << std::endl;
        std::cout << "-----------------------------------------------------------------------------" << std::endl;
        std::cout << " Solving the TG analysis...                                                  " << std::endl;
        std::cout << "-----------------------------------------------------------------------------" << std::endl;
    }

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

    if (biosmoke_options_.verbose_video() == true)
    {
        std::cout << std::endl;
        std::cout << "-----------------------------------------------------------------------------" << std::endl;
        std::cout << " Completed the simulation in " << std::setprecision(6) << tEnd - tStart << " seconds"
                  << std::endl;
        std::cout << "-----------------------------------------------------------------------------" << std::endl;
    }

    // CloseAllFiles(); //TODO
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
                std::cout << std::setw(10) << std::left << "M/M0[-]";
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
                }
            }

            // XML file output
            if (biosmoke_options_.verbose_xml_file() == true)
            {
                if (iteration_ % biosmoke_options_.n_step_file() == 1 || biosmoke_options_.n_step_file() == 1 ||
                    t == final_time_)
                {
                    counter_file_XML_++;
                }
            }
        }
    }
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

} // namespace BioSMOKE