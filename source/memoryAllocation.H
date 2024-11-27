// TGA analysis
std::string analysis;   // variable used also for TotalAnalysis
double Heat_Rates;
double final_time;      // variable used also for TotalAnalysis
std::vector<std::string> output_species_;
bool temperature_profile = false;
OpenSMOKE::BiomassTemperature_Profile*					temperature_profile_;


// Elenco variabili che in memoryAllocation.H vengono solo inizializzate e 
// vengono riempite tramite i valori letti da input.dic nel file readInputFile.H
// - porosità  ->  epsi
// - torutuosità -> tau
// - n° di gusci -> intervalli
// - raggio iniziale della sfera -> raggioTot
// - tempo finale di integrazione -> final_time
// - booleana per calcolo eq. di bilacio energia (ON/OFF) -> energyBalance
// - booleana per calcolo variazione di volume sfera (ON/OFF) -> volumeLoss

bool energyBalance;
bool volumeLoss;
double epsi;
double tau;
int intervalli;
double raggioTot;
double Da;

//Elenco delle varaibili che servono nella lettura dei file di input in readInputFile.H
//- fase gas
double Tbulk;
double Pbulk;
OpenSMOKE::OpenSMOKEVectorDouble omegaIn_gas;
double rhoGas;
OpenSMOKE::OpenSMOKEVectorDouble cGas_TGA;
OpenSMOKE::OpenSMOKEVectorDouble xGasIn;

//- fase solida
double Tsolid;
double Psolid;
double Tsurf_solid;	// Temperature of the biomass surface (employed for flux balances at the biomass-gas interface)
double hext;
double rhoSolid;
OpenSMOKE::OpenSMOKEVectorDouble omegaIn_solid;
double lambdaS;


// Elenco di variabili che sono dichiarate come variabili globali in 
// memoryAllocation.H e che vengono definite una volta per tutte alla fine di
// readInputFile.H
int solid_species;
int gas_species;
int solid_reactions;
int gas_reactions;
int Neq;
int Ntot_eq;


OpenSMOKE::OpenSMOKEVectorDouble P;   // update pressure in every shell
OpenSMOKE::OpenSMOKEVectorDouble MW_tot;
OpenSMOKE::OpenSMOKEVectorDouble MW_solid;
OpenSMOKE::OpenSMOKEMatrixDouble nu;


//Elenco delle variabili per caratterizzare la geometria della sfera
OpenSMOKE::OpenSMOKEVectorDouble raggio;
OpenSMOKE::OpenSMOKEVectorDouble r;
OpenSMOKE::OpenSMOKEVectorDouble V;
OpenSMOKE::OpenSMOKEVectorDouble S;
double Vtot;
double Vtot_iniz;
// Update value of epsi 
OpenSMOKE::OpenSMOKEVectorDouble epsi_var;

double massTotsolid_initial;  // massa totale nella particella
OpenSMOKE::OpenSMOKEVectorDouble massTotsolid_initial_parz; // massa totale iniziale di solido in ogni shell

// Variabili di ODESystem

OpenSMOKE::OpenSMOKEMatrixDouble mass;
OpenSMOKE::OpenSMOKEMatrixDouble dmass;
OpenSMOKE::OpenSMOKEVectorDouble T;
OpenSMOKE::OpenSMOKEVectorDouble dT;

// Variabili che servono all'inzio di ODESystem
OpenSMOKE::OpenSMOKEMatrixDouble mole;
OpenSMOKE::OpenSMOKEVectorDouble massSolid;  // vettore con massa delle singole species solide
OpenSMOKE::OpenSMOKEVectorDouble moleSolid;  // vettore con mole delle singole species solide
OpenSMOKE::OpenSMOKEVectorDouble massGas;
OpenSMOKE::OpenSMOKEVectorDouble moleGas;
OpenSMOKE::OpenSMOKEVectorDouble massTot_solid;  // vettore con massa dei singoli intervalliniF
OpenSMOKE::OpenSMOKEVectorDouble massTot_gas;
OpenSMOKE::OpenSMOKEVectorDouble moleTot_solid;
OpenSMOKE::OpenSMOKEVectorDouble moleTot_gas;
OpenSMOKE::OpenSMOKEMatrixDouble massFraction;
OpenSMOKE::OpenSMOKEMatrixDouble moleFraction;
OpenSMOKE::OpenSMOKEMatrixDouble massFraction_Shellmassbased;

// Variabili che servono in properties
OpenSMOKE::OpenSMOKEVectorDouble omegaGasFraction;
OpenSMOKE::OpenSMOKEVectorDouble omegaSolidFraction;
OpenSMOKE::OpenSMOKEVectorDouble xGasFraction;
OpenSMOKE::OpenSMOKEVectorDouble xSolidFraction;
double MW_gas;  // elemento da passare alle openSmoke
double MW_gas_in;
OpenSMOKE::OpenSMOKEVectorDouble MW_gas_mix;
OpenSMOKE::OpenSMOKEVectorDouble MW_solid_mix;
double cP_gas;  // elemento da passare alle openSmoke
OpenSMOKE::OpenSMOKEVectorDouble cPgas_mix ;
double thermal_conductivity;   // elemento da passare alle openSmoke
double viscosity; // elemento da passare alle openSmoke
OpenSMOKE::OpenSMOKEVectorDouble lambda_gas ;
OpenSMOKE::OpenSMOKEVectorDouble dyn_viscosity ;
OpenSMOKE::OpenSMOKEVectorDouble Diff;  // elemento da passare alle openSmoke settato su species gas
OpenSMOKE::OpenSMOKEMatrixDouble Dmix;
OpenSMOKE::OpenSMOKEMatrixDouble Diff_Knudsen;
OpenSMOKE::OpenSMOKEMatrixDouble Diff_eff;
OpenSMOKE::OpenSMOKEVectorDouble Hi;
OpenSMOKE::OpenSMOKEMatrixDouble Hi_mass;
OpenSMOKE::OpenSMOKEVectorDouble Hi_tot;   // gas+solido  // MOMENTANTANEE
OpenSMOKE::OpenSMOKEVectorDouble Si_tot;
OpenSMOKE::OpenSMOKEVectorDouble Hi_tot_corr;   // solido+gas  MOMENTANEE
OpenSMOKE::OpenSMOKEMatrixDouble Hi_mass_tot; //gas+solido     MOMENTANEE
OpenSMOKE::OpenSMOKEVectorDouble Hmix; 
OpenSMOKE::OpenSMOKEVectorDouble C_gas;
OpenSMOKE::OpenSMOKEVectorDouble rho_gas;
OpenSMOKE::OpenSMOKEVectorDouble cp_tot; // elemento da passare alle openSmoke settato su species solid+gas_species
OpenSMOKE::OpenSMOKEVectorDouble cp_gas_species;
OpenSMOKE::OpenSMOKEVectorDouble cp_solid_species; 
OpenSMOKE::OpenSMOKEVectorDouble cPsolid_mix ;
OpenSMOKE::OpenSMOKEVectorDouble lambda_solid_species;
OpenSMOKE::OpenSMOKEVectorDouble lambda_solid ;
OpenSMOKE::OpenSMOKEVectorDouble lambda_effective ;
//Variabili che servono in kinetics
OpenSMOKE::OpenSMOKEMatrixDouble FormationRate;
OpenSMOKE::OpenSMOKEMatrixDouble FormationRate_old;
OpenSMOKE::OpenSMOKEVectorDouble sumRgas;
OpenSMOKE::OpenSMOKEVectorDouble sumRsolid;
OpenSMOKE::OpenSMOKEVectorDouble sumMoleFormationRates;
OpenSMOKE::OpenSMOKEVectorDouble HeatReaction;
double HeatReactionTot;
OpenSMOKE::OpenSMOKEMatrixDouble DHr; //momentanea
OpenSMOKE::OpenSMOKEMatrixDouble Vel_gas;
OpenSMOKE::OpenSMOKEVectorDouble RGG;

// material and energy balance
double ExtGas_cp;
double ExtGas_lambda;

OpenSMOKE::OpenSMOKEMatrixDouble JD_out;
OpenSMOKE::OpenSMOKEMatrixDouble JD_in;
OpenSMOKE::OpenSMOKEVectorDouble massConvection_in;
OpenSMOKE::OpenSMOKEVectorDouble massConvection_out;
OpenSMOKE::OpenSMOKEVectorDouble JC_out;
OpenSMOKE::OpenSMOKEVectorDouble JC_in;
OpenSMOKE::OpenSMOKEVectorDouble JDH_in;
OpenSMOKE::OpenSMOKEVectorDouble JDH_out;
OpenSMOKE::OpenSMOKEVectorDouble flux_in;
OpenSMOKE::OpenSMOKEVectorDouble flux_out;

OpenSMOKE::OpenSMOKEVectorDouble residual;
double norm2;

//Print
OpenSMOKE::OpenSMOKEVectorDouble residual_CELL;
OpenSMOKE::OpenSMOKEVectorDouble residual_CELLA;
OpenSMOKE::OpenSMOKEVectorDouble residual_CHAR;
OpenSMOKE::OpenSMOKEVectorDouble residual_ASH;

OpenSMOKE::OpenSMOKEVectorDouble residual_CO;
OpenSMOKE::OpenSMOKEVectorDouble residual_CO2;
OpenSMOKE::OpenSMOKEVectorDouble residual_H2;
OpenSMOKE::OpenSMOKEVectorDouble residual_H2O;
OpenSMOKE::OpenSMOKEVectorDouble residual_N2;
OpenSMOKE::OpenSMOKEVectorDouble residual_TGL;
OpenSMOKE::OpenSMOKEVectorDouble residual_TANN;
OpenSMOKE::OpenSMOKEVectorDouble massTGL;
OpenSMOKE::OpenSMOKEVectorDouble  massTANN;
OpenSMOKE::OpenSMOKEVectorDouble  massITANN;

std::string gg;  //prova

OpenSMOKE::OpenSMOKEVectorDouble mass_TGA;
double T_TGA;
OpenSMOKE::OpenSMOKEVectorDouble omegaTot_TGA;
OpenSMOKE::OpenSMOKEVectorDouble omegaSolid_TGA;
OpenSMOKE::OpenSMOKEVectorDouble omegaGas_TGA;


//OpenSMOKE::OpenSMOKEVectorDouble enthalpyReaz (24, 0.,109500.,100900.,-268000.,-10900.,-600.,-181500.,77600.,50200.,-43000.,128300.,108600.,-130000.,-109400.,142600.,-92800.,-345900.,-121600.,-56200.,203300.,0.,0.,0.,0.);

OpenSMOKE::OpenSMOKEVectorDouble enthalpyReaz (6, 0.,109500.,100900.,-268000.,0.,0.);


OpenSMOKE::OpenSMOKEVectorDouble GasFormato;


double massFinalGas;
double massGasPerfect;
OpenSMOKE::OpenSMOKEVectorDouble Kc;
OpenSMOKE::OpenSMOKEVectorDouble massFinalGasSpecies;
OpenSMOKE::OpenSMOKEVectorDouble massFinalGasSpecies_total;
OpenSMOKE::OpenSMOKEVectorDouble moleFinalGasSpecies_total;
double massFinalSolid;

OpenSMOKE::OpenSMOKEVectorDouble massSolidSpeciesTotal;
OpenSMOKE::OpenSMOKEVectorDouble massGasSpeciesTotal;
OpenSMOKE::OpenSMOKEVectorDouble massGasSpecies_outIntegrated;
OpenSMOKE::OpenSMOKEVectorDouble massSolidSpeciesTotal_fromR;
OpenSMOKE::OpenSMOKEVectorDouble massGasSpeciesTotal_fromR;
OpenSMOKE::OpenSMOKEVectorDouble massGasSpeciesFlowrate;
OpenSMOKE::OpenSMOKEVectorDouble massGasSpeciesFlowrate_old;
double t_old = 0.;

OpenSMOKE::OpenSMOKEVectorDouble omegaI;

double timeOld = 0.;
double dt= 0.;
double massDiffusion_out = 0.;

double massConv_first = 0.;
double massDiff_first = 0.;

OpenSMOKE::OpenSMOKEVectorDouble Vrecord;

OpenSMOKE::OpenSMOKEVectorDouble massSolid_interval;

OpenSMOKE::OpenSMOKEVectorDouble Vstart;
OpenSMOKE::OpenSMOKEVectorDouble Rstart;
OpenSMOKE::OpenSMOKEVectorDouble Vstore;
OpenSMOKE::OpenSMOKEVectorDouble Ugrid;
OpenSMOKE::OpenSMOKEVectorDouble dmSolid;

double ele;
double massASH;
double massCELL;


OpenSMOKE::OpenSMOKEVectorDouble print(6);

int count = 0.;

int inertIndex;
int tarcount;
std::vector<int> lightGasIndex;
std::vector<int> tarIndex;
std::vector<int> residualIndex;
double gasFlowTempLG;
double gasFlowTempTAR;
double massResidualSolid;


double yieldChar;


unsigned int n_steps_video_ = 100; //10
unsigned int count_ode_video_ = 1;
unsigned int n_steps_file_ = 5; //30
unsigned int count_file_ = 1;   //30