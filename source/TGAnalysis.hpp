// #include "TGAnalysis.h"

#include <interfaces/Interface_Dense_OpenSMOKEppOde.h>

namespace BioSMOKE
{
// clang-format off

TGAnalysis::TGAnalysis(OpenSMOKE::ThermodynamicsMap_CHEMKIN &thermodynamicsMap,
                       OpenSMOKE::KineticsMap_CHEMKIN &kineticsMap,
                       OpenSMOKE::TransportPropertiesMap_CHEMKIN &transportMap,
                       OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN &thermodynamicsSolidMap,
                       OpenSMOKE::KineticsMap_Solid_CHEMKIN &kineticsSolidMap,
                       OpenSMOKE::ODE_Parameters &ode_parameters,
                       //OpenSMOKE::BioSMOKE_Options &biosmoke_options,
                       const double T0,
                       const double P0,
                       const double rho0_solid,
                       const std::vector<double> &omega0_gas,
                       const std::vector<double> &omega0_solid,
                       const double heating_rate) :
            BaseSolver(thermodynamicsMap, kineticsMap, transportMap, thermodynamicsSolidMap, kineticsSolidMap, ode_parameters)
{ // clang-format on

    iteration_ = 0;
    counter_file_video_ = 0;
    counter_file_ASCII_ = 0;

    V0_solid_ = 1;
    T0_ = T0;
    P0_ = P0;
    omega0_solid_ = omega0_solid;
    omega0_gas_ = omega0_gas;
    heating_rate_ = heating_rate;

    NGS_ = thermodynamicsMap_.NumberOfSpecies();
    NSS_ = thermodynamicsSolidMap_.number_of_solid_species();
    NC_ = NGS_ + NSS_;
    NE_ = NC_ + 1;

    thermodynamicsMap_.SetTemperature(T0_);
    thermodynamicsMap_.SetPressure(P0_);
    kineticsMap_.SetTemperature(T0_);
    kineticsMap_.SetPressure(P0_);

    x0_gas_.resize(NGS_);
    thermodynamicsMap_.MoleFractions_From_MassFractions(x0_gas_.data(), MW0_gas_, omega0_gas_.data());
    rho0_gas_ = P0_ * MW0_gas_ / (PhysicalConstants::R_J_kmol * T0_);
    MW_gas_ = MW0_gas_;
    rho_gas_ = rho0_gas_;
    mass0_tot_gas_ = 0.;

    thermodynamicsSolidMap_.SetTemperature(T0_);
    thermodynamicsSolidMap_.SetPressure(P0_);
    kineticsSolidMap_.SetTemperature(T0_);
    kineticsSolidMap_.SetPressure(P0_);
    x0_solid_.resize(NSS_);
    thermodynamicsSolidMap_.SolidMoleFractions_From_SolidMassFractions(x0_solid_.data(), MW0_solid_,
                                                                       omega0_solid_.data());
    MW_solid_ = MW0_solid_;
    V_solid_ = V0_solid_;
    rho_solid_ = rho0_solid;
    mass0_tot_solid_ = V0_solid_ * rho0_solid_;
    mass_tot_solid_ = mass0_tot_solid_;

    T_ = T0_;
    P_ = P0_;
}

int TGAnalysis::Equations(const double t, const std::vector<double> &y, std::vector<double> &dy)
{
    // recover unknowns: mass_gas <> mass_solid <> T
    std::vector<double> mass_solid_current_(NSS_, 0.);
    std::vector<double> mass_gas_current_(NGS_, 0.);
    for (unsigned int i = 0; i < NE_; i++)
    {
        if (i < NGS_)
            mass_gas_current_[i] = y[i];
        else if (i < NGS_ + NSS_)
            mass_solid_current_[i - NGS_] = y[i];
        else
            T_ = y[i];
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
    kineticsSolidMap_.ReactionEnthalpiesAndEntropies();
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
    typedef OdeSMOKE::KernelDense<OpenSMOKE::ODESystem_OpenSMOKE_ThermogravimetricAnalysis> denseOde;
    typedef OdeSMOKE::MethodGear<denseOde> methodGear;
    OdeSMOKE::MultiValueSolver<methodGear> ode_solver;
    // ode_solver.SetThermogravimetricAnalysis(this);
    //  TODO: bioSMOKE options
    //  if (biosmoke_options_.verbose_video() == true)
    //  {
    //      std::cout << std::endl;
    //      std::cout << "-----------------------------------------------------------------------------" << std::endl;
    //      std::cout << " Solving the TG analysis...                                                  " << std::endl;
    //      std::cout << "-----------------------------------------------------------------------------" << std::endl;
    //  }

    // TODO: Implement OdeInterface for BioSMOKE
    std::cout << "Solving TG analysis..." << std::endl;
}

int TGAnalysis::Print(const double t, const std::vector<double> &y)
{
    iteration_++;

    // if (biosmoke_options_.verbose_video() == true)
    // {
    //     if (iteration_ % biosmoke_options_.n_step_video() == 1 || biosmoke_options_.n_step_video() == 1)
    //     {
    //         counter_file_video_++;
    //         if (counter_file_video_ % 100 == 1)
    //         {
    //             std::cout << std::endl;
    //         }
    //     }
    // }

    // if (biosmoke_options_.sensitivity_analysis() == true)
    //     SensitivityAnalysis(t, y);
    std::cout << "Printing TG analysis..." << std::endl;
    return 0;
}

} // namespace BioSMOKE