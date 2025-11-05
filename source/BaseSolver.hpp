namespace BioSMOKE
{

// clang-format off
BaseSolver::BaseSolver(OpenSMOKE::ThermodynamicsMap_CHEMKIN &thermodynamicsMap,
                       OpenSMOKE::KineticsMap_CHEMKIN &kineticsMap,
                       OpenSMOKE::TransportPropertiesMap_CHEMKIN &transportMap,
                       OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN &thermodynamicsSolidMap,
                       OpenSMOKE::KineticsMap_Solid_CHEMKIN &kineticsSolidMap,
                       OpenSMOKE::ODE_Parameters &ode_parameters,
                       BioSMOKE::BioSMOKE_Options &biosmoke_options) :
                       thermodynamicsMap_(thermodynamicsMap),
                       kineticsMap_(kineticsMap),
                       transportMap_(transportMap),
                       thermodynamicsSolidMap_(thermodynamicsSolidMap), 
                       kineticsSolidMap_(kineticsSolidMap),
                       ode_parameters_(ode_parameters),
                       biosmoke_options_(biosmoke_options)
{
} // clang-format on

BaseSolver::~BaseSolver() {}

void BaseSolver::SetTemperatureProfile(BioSMOKE::BioSMOKE_Profile &biosmoke_profile)
{
    is_temperature_profile_ = true;
    biosmoke_profile_ = &biosmoke_profile;
    std::cout << "The temperature profile has been set in the solver." << std::endl;
};

} // namespace BioSMOKE