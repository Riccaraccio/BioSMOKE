////////////////////////////////////////////////////////////////////////
//                      FUNCTIONs DEFINITION                          //
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//                         TGA FUNCTIONs                              //
////////////////////////////////////////////////////////////////////////

void OdeTGA(const double* y, const double t, double* dy, void* args)
{
    ChangeDimensions(solid_species + gas_species, &mass_TGA, true);
    ChangeDimensions(solid_species + gas_species, &omegaTot_TGA, true);
    ChangeDimensions(solid_species, &omegaSolid_TGA, true);
    ChangeDimensions(gas_species, &omegaGas_TGA, true);

    if (temperature_profile == true)
        Heat_Rates = temperature_profile_->Get_increase(t);

    int k = 0;
    for (int j = 1; j <= Neq; j++)
    {
        if (j <= (gas_species + solid_species))
            mass_TGA[j] = y[k++];
        else
            T_TGA = y[k++];
    }

    massFinalSolid = 0.;
    for (int j = 1; j <= (solid_species); j++)
        massFinalSolid += mass_TGA[j];

    massFinalGas = 0.;
    for (int j = solid_species + 1; j <= (solid_species + gas_species); j++)
        massFinalGas += mass_TGA[j];

    for (int j = 1; j <= (solid_species + gas_species); j++)
        omegaTot_TGA[j] = mass_TGA[j] / (massFinalSolid + massFinalGas);

    for (int j = 1; j <= (solid_species); j++)
        omegaSolid_TGA[j] = mass_TGA[j] / (massFinalSolid);

    for (int j = solid_species + 1; j <= (solid_species + gas_species); j++)
        omegaGas_TGA[j - solid_species] = mass_TGA[j] / (massFinalGas);

    // Concentrations evaluation

    for (int j = 1; j <= gas_species; j++) 
    {
        cGas_TGA[j] = Psolid / 8.314 / T_TGA * xGasIn[j]; // Cgas calcualted using fixed gas composition = inlet gas composition 
    }

    OpenSMOKE::OpenSMOKEVectorDouble cSolid_TGA(solid_species);
    for (int j = 1; j <= solid_species; j++)
    {
        cSolid_TGA[j] = rhoSolid * omegaSolid_TGA[j] / MW_tot[j]; // LOCASPI NEW FORUMLA, more general but needs correction therm in the dy eqs
        //cSolid_TGA[j] = mass_TGA[j] / MW_tot[j]; // GENTILE OLD FORULA


    }

    const double T_K = T_TGA;
    const double P_Pa = Psolid;

    thermodynamicsSolidMapXML->SetTemperature(T_K);
    thermodynamicsSolidMapXML->SetPressure(P_Pa);
    kineticsSolidMapXML[0]->SetPressure(P_Pa);
    kineticsSolidMapXML[0]->SetTemperature(T_K);
    kineticsSolidMapXML[0]->ReactionEnthalpiesAndEntropies();
    kineticsSolidMapXML[0]->ReactionRates(cGas_TGA.GetHandle(), cSolid_TGA.GetHandle()); //kmol/m3/s

    OpenSMOKE::OpenSMOKEVectorDouble RGas(thermodynamicsSolidMapXML->number_of_gas_species());
    OpenSMOKE::OpenSMOKEVectorDouble RSolid(thermodynamicsSolidMapXML->number_of_solid_species());
    kineticsSolidMapXML[0]->FormationRates(RGas.GetHandle(), RSolid.GetHandle());

    for (int i = 0; i < Neq; i++)
    {
        if (i < solid_species)
            dy[i] = RSolid[i + 1] * MW_tot[i + 1] * massFinalSolid / rhoSolid; // *massFinalSolid/rhoSolid not needed for Gentile's old formula
        else if (i >= solid_species && i < solid_species + gas_species)
            dy[i] = RGas[i + 1 - solid_species] * MW_tot[i + 1] * massFinalSolid / rhoSolid;
        else
            dy[i] = Heat_Rates;
    }
}

void MyStepPrintTGA(const double t, const Eigen::VectorXd& Y)
{
    dt = t - timeOld;
    timeOld = t;

    unsigned int NS_ = thermodynamicsSolidMapXML->NumberOfSpecies();
    unsigned int NSolid_ = thermodynamicsSolidMapXML->number_of_solid_species();
    unsigned int NGas_ = thermodynamicsSolidMapXML->number_of_gas_species();

    if (count_ode_video_ % n_steps_video_ == 1 || t == final_time)
    {
        if (count_ode_video_ % (n_steps_video_ * 1000) == 1)
        {
            std::cout << std::endl;
            std::cout << std::setw(18) << std::left << "#Step";
            std::cout << std::setw(18) << std::left << "Time[s]";
            std::cout << std::setw(18) << std::left << "T[K]";
            std::cout << std::setw(18) << std::left << "TG[-]";
            std::cout << std::endl;
        }

        std::cout << std::fixed << std::setw(18) << std::left << count_ode_video_;
        std::cout << std::scientific << std::setw(18) << std::setprecision(6) << std::left << t;
        std::cout << std::scientific << std::setw(18) << std::setprecision(6) << std::left << T_TGA;
        std::cout << std::scientific << std::setw(18) << std::setprecision(6) << std::left << massFinalSolid / massTotsolid_initial;
        std::cout << std::endl;
    }

    if (count_file_ == n_steps_file_ || t == final_time)
    {

        TG << std::setw(25) << std::left << t << std::setw(25) << std::left << t / 60
            << std::setw(25) << std::left << T_TGA - 273.15 << std::setw(25) << std::left << T_TGA
            << std::setw(25) << std::left << massFinalSolid / massTotsolid_initial << std::endl;


        yieldSpecies << std::setw(25) << std::left << t << std::setw(25) << std::left << t / 60
            << std::setw(25) << std::left << T_TGA - 273.15 << std::setw(25) << std::left << T_TGA;

        for (int k = 0; k < output_species_.size(); k++)
        {
            std::string type = output_species_[k];
            int index = thermodynamicsSolidMapXML->IndexOfSpecies(type);
            if (index <= NGas_)
            {
                yieldSpecies << std::setw(25) << std::left << omegaTot_TGA[index + NSolid_];
            }
            else
            {
                yieldSpecies << std::setw(25) << std::left << omegaTot_TGA[index - NGas_];
            }
        }

        yieldSpecies << std::endl;

        count_file_ = 0;
    }
    count_ode_video_++;
    count_file_++;
}


////////////////////////////////////////////////////////////////////////
//                         1D FUNCTIONs                               //
////////////////////////////////////////////////////////////////////////

void OdeTotal(const double* y, const double t, double* dy, void* args)
{


    // Unknowns reconstruction
    int k = 0;
    for (int i = 1; i <= intervalli; i++)
    {
        for (int j = 1; j <= Neq; j++)
        {
            if (j <= (gas_species + solid_species))
                mass[i][j] = y[k++];
            else
                T[i] = y[k++];
        }
    }

    if (temperature_profile == true)
    {
        //Heat_Rates = temperature_profile_->Get_increase(t);
        Tbulk = temperature_profile_->Get_temperature(t);
    }

    //T[intervalli] += 10000*t;
    //T[intervalli] = (487.9*exp(0.0006094*t) - 404.5*exp(-0.0216*t))+273.15;
    //T[intervalli] = 473.15;
    //T[intervalli] = 866 - (866 - 300)*exp(-t / 50) + 700 * (1 - exp(-t / 50))*exp(-t / 25); //BEST CORRELATION FOR PARK
    //T[intervalli] = 742-(742-300)*exp(-t/50)+500*(1-exp(-t/20))*exp(-t/25);
    //Tsup_solid= 688-(688-300)*exp(-t/70)+700*(1-exp(-t/200))*exp(-t/18);/*300 + 1.667/4*t;*/

    //T[intervalli] = Tsup_solid;
    //T[intervalli] = 875.5*exp(-4.121e-5*t)-567.8*exp(-0.004238*t); //Tsup beech wood diblasi2003


    for (int i = 1; i <= intervalli; i++)
        for (int j = 1; j <= (solid_species + gas_species); j++)
            mole[i][j] = mass[i][j] / MW_tot[j];

    CreateFractions();
    SetVolume(t);
    CalculateProperties();
    CalculateKinetics();

    CreateMaterialBalances();
    CreateEnergyBalances();


    // Equations

    k = 0;
    for (int i = 1; i <= intervalli; i++)
    {
        for (int j = 1; j <= Neq; j++)
        {
            if (j <= (solid_species + gas_species))
                dy[k++] = dmass[i][j];
            else
                dy[k++] = dT[i];
        }
    }

       

    for (int i = 1; i <= (intervalli * Neq); i++)
        residual[i] = dy[i-1];

    norm2 = residual.Norm2();

    massFinalSolid = 0.;
    for (int i = 1; i <= intervalli; i++)
    {
        for (int j = 1; j <= (solid_species); j++)
            massFinalSolid += mass[i][j];
    }

    massFinalGas = 0.;
    for (int i = 1; i <= intervalli; i++)
    {
        for (int j = solid_species + 1; j <= (solid_species + gas_species); j++)
            massFinalGas += mass[i][j];
    }

    massFinalGasSpecies = 0.;
    for (int i = solid_species + 1; i <= (solid_species + gas_species); i++)
    {
        for (int j = 1; j <= intervalli; j++)
        {
            massFinalGasSpecies[i - solid_species] += mass[j][i];
        }
    }


    for (int i = solid_species + 1; i <= (solid_species + gas_species); i++)
    {
        for (int j = 1; j <= intervalli; j++)
        {
            massFinalGasSpecies_total[i - solid_species] += mass[j][i];
        }
    }

    for (int i = solid_species + 1; i <= (solid_species + gas_species); i++)
    {
        for (int j = 1; j <= intervalli; j++)
        {
            moleFinalGasSpecies_total[i - solid_species] += mole[j][i];
        }
    }

    massGasPerfect = 0.;
    for (int i = 1; i <= intervalli; i++)
    {
        massGasPerfect += P[i] * V[i] * epsi_var[i] * MW_gas_mix[i] / PhysicalConstants::R_J_kmol / T[i];
    }

    for (int i = 1; i <= intervalli; i++)
    {
        massSolid_interval[i] = 0.;
        for (int j = 1; j <= (solid_species); j++)
        {
            massSolid_interval[i] += mass[i][j];
        }
    }

    /*massASH = 0.;
    massCELL = 0.;


    for (int i = 1; i<= intervalli;i++)
    {

            massASH += mass[i][thermodynamicsSolidMapXML->IndexOfSpecies("ASH")-gas_species];
            massCELL += mass[i][thermodynamicsSolidMapXML->IndexOfSpecies("CELL")-gas_species];
            //massTANN += mass[i][thermodynamicsSolidMapXML->IndexOfSpecies("TANN")-gas_species];
            //massTGL += mass[i][thermodynamicsSolidMapXML->IndexOfSpecies("TGL")-gas_species];

    }

    residual_CELL = 0.;
    residual_CELLA = 0.;
    residual_CHAR = 0.;
    residual_ASH = 0.;
    massTANN = 0.;
    massTGL = 0.;

    residual_CO = 0.;
    residual_CO2 = 0.;
    residual_H2 = 0.;
    residual_H2O = 0.;
    residual_N2 = 0.;
    for (int i = 1; i<= intervalli; i++)
    {
        residual_CELL[i] = massFraction[i][thermodynamicsSolidMapXML->IndexOfSpecies("CELL")-gas_species];
        residual_CELLA[i] = massFraction[i][thermodynamicsSolidMapXML->IndexOfSpecies("CELLA")-gas_species];
        residual_CHAR[i] = massFraction[i][thermodynamicsSolidMapXML->IndexOfSpecies("CHAR")-gas_species];
        residual_ASH[i] = massFraction[i][thermodynamicsSolidMapXML->IndexOfSpecies("ASH")-gas_species];
        //residual_TANN[i] = massFraction[i][thermodynamicsSolidMapXML->IndexOfSpecies("TANN")-gas_species];
        //residual_TGL[i] = massFraction[i][thermodynamicsSolidMapXML->IndexOfSpecies("TGL")-gas_species];
        //massTANN[i] = mass[i][thermodynamicsSolidMapXML->IndexOfSpecies("TANN")-gas_species];
        //massITANN[i] = mass[i][thermodynamicsSolidMapXML->IndexOfSpecies("ITANN")-gas_species];
        //massTGL[i] = mass[i][thermodynamicsSolidMapXML->IndexOfSpecies("TGL")-gas_species];

        residual_CO[i] = massFraction[i][thermodynamicsMapXML->IndexOfSpecies("CO")+solid_species];
        residual_CO2[i] = massFraction[i][thermodynamicsMapXML->IndexOfSpecies("CO2")+solid_species];
        residual_H2[i] = massFraction[i][thermodynamicsMapXML->IndexOfSpecies("H2")+solid_species];
        residual_H2O[i] = massFraction[i][thermodynamicsMapXML->IndexOfSpecies("H2O")+solid_species];
        residual_N2[i] = massFraction[i][thermodynamicsMapXML->IndexOfSpecies("N2")+solid_species];
    }
    */
}

void MyStepPrintTotal(const double t, const Eigen::VectorXd& x)
{

    unsigned int NS_ = thermodynamicsSolidMapXML->NumberOfSpecies();
    unsigned int NSolid_ = thermodynamicsSolidMapXML->number_of_solid_species();
    unsigned int NGas_ = thermodynamicsSolidMapXML->number_of_gas_species();


    // Integration
    massSolidSpeciesTotal = 0.;
    massGasSpeciesTotal = 0.;
    double tot = 0.;

    for (int i = 1; i <= intervalli; i++)
    {
        for (int j = 1; j <= (solid_species + gas_species); j++)
        {
            if (j <= solid_species)
            {
                massSolidSpeciesTotal[j] += mass[i][j];
                massSolidSpeciesTotal_fromR[j] += (t - t_old) * (FormationRate[i][j] + FormationRate_old[i][j]) / 2. * MW_tot[j] * V[i] * (1 - epsi_var[i]);
            }
            else
            {
                massGasSpeciesTotal[j - solid_species] += mass[i][j];
                massGasSpeciesTotal_fromR[j - solid_species] += (t - t_old) * (FormationRate[i][j] + FormationRate_old[i][j]) / 2. * MW_tot[j] * V[i] * (1 - epsi_var[i]);
            }

            /*if (count_file_ == n_steps_file_ && j>solid_species)
                std::cout << FormationRate[i][j] << " " << FormationRate_old[i][j] << " " << (t - t_old) * (FormationRate[i][j] + FormationRate_old[i][j]) / 2. << " " << massGasSpeciesTotal_fromR[j - solid_species] << std::endl;*/

            tot += mass[i][j];

        }
    }
    //std::cout << std::endl;


    for (int j = 1; j <= gas_species; j++)
    {
        massGasSpeciesFlowrate[j] = massConvection_out[intervalli] * massFraction[intervalli][j + solid_species] + JD_out[intervalli][j];

        massGasSpecies_outIntegrated[j] += (t - t_old) * (massGasSpeciesFlowrate_old[j] + massGasSpeciesFlowrate[j]) / 2.;
        tot += massGasSpecies_outIntegrated[j];
    }

    t_old = t;
    for (int j = 1; j <= (solid_species + gas_species); j++)
    {
        for (int i = 1; i <= intervalli; i++)
            FormationRate_old[i][j] = FormationRate[i][j];
        if (j > solid_species)
            massGasSpeciesFlowrate_old[j - solid_species] = massGasSpeciesFlowrate[j - solid_species];
    }

    if (count_ode_video_ % n_steps_video_ == 1 || t == final_time)
    {
        if (count_ode_video_ % (n_steps_video_ * 1000) == 1)
        {
            std::cout << std::endl;
            std::cout << std::setw(18) << std::left << "#Step";
            std::cout << std::setw(18) << std::left << "Time[s]";
            std::cout << std::setw(18) << std::left << "TG[-]";
            std::cout << std::endl;
        }

        std::cout << std::fixed << std::setw(18) << std::left << count_ode_video_;
        std::cout << std::scientific << std::setw(18) << std::setprecision(6) << std::left << t;
        std::cout << std::scientific << std::setw(18) << std::setprecision(6) << std::left << massFinalSolid / massTotsolid_initial;

        std::cout << std::endl;
    }


    if (count_file_ == n_steps_file_ || t == final_time)
    {
        const int width = 20;

        press << std::setw(25) << std::left << t;
        for (int j = 1; j <= intervalli; j++)
        {
            press << std::setw(25) << std::left << P[j];
        }
        press << "\n";

        fmass << std::setw(25) << std::left << t;
        for (int j = 1; j <= intervalli; j++)
        {
            fmass << std::setw(25) << std::left << massTot_solid[j];
        }
        fmass << "\n";

        por << std::setw(25) << std::left << t;
        for (int j = 1; j <= intervalli; j++)
        {
            por << std::setw(25) << std::left << epsi_var[j];
        }
        por << "\n";

        QreactShell << std::setw(25) << std::left << t;
        for (int j = 1; j <= intervalli; j++)
        {
            QreactShell << std::setw(25) << std::left << HeatReaction[j];
        }
        QreactShell << "\n";

        QreactTot << std::setw(25) << std::left << t
            << std::setw(25) << std::left << HeatReactionTot << "\n";

        temp << std::setw(25) << std::left << t;
        for (int j = 1; j <= intervalli; j++)
        {
            temp << std::setw(25) << std::left << T[j];
        }
        temp << "\n";

        volume << std::setw(25) << std::left << t
            << std::setw(25) << std::left << Vtot / Vtot_iniz
            << std::setw(25) << std::left << pow(Vtot / Vtot_iniz, 0.33333);
        for (int j = 1; j <= intervalli; j++)
        {
            volume << std::setw(25) << std::left << V[j];
        }
        volume << "\n";

        TG << std::setw(25) << std::left << t
            << std::setw(25) << std::left << massFinalSolid / massTotsolid_initial << "\n";


        gasFlowRate << std::setw(25) << std::left << t;
        for (int i = 1; i <= (solid_species + gas_species); i++)
        {
            if (i > solid_species)
                gasFlowRate << std::setw(25) << std::left << massGasSpeciesFlowrate[i - solid_species];

        }
        gasFlowRate << std::endl;


        yieldSpecies << std::setw(25) << std::left << t;
        yieldSpecies << std::setw(25) << std::left << (massSolidSpeciesTotal_fromR.SumElements() + massGasSpeciesTotal_fromR.SumElements()) / massTotsolid_initial;

        for (int j = 1; j <= solid_species; j++)
            yieldSpecies << std::setw(25) << std::left << massSolidSpeciesTotal_fromR[j] / massTotsolid_initial;

        for (int j = 1; j <= gas_species; j++)
            yieldSpecies << std::setw(25) << std::left << massGasSpeciesTotal_fromR[j] / massTotsolid_initial;

        yieldSpecies << std::endl;

        fmassTotal << std::setw(25) << std::left << t;
        fmassTotal << std::setw(25) << std::left << tot;

        for (int j = 1; j <= solid_species; j++)
            fmassTotal << std::setw(25) << std::left << massSolidSpeciesTotal[j];
        for (int j = 1; j <= gas_species; j++)
            fmassTotal << std::setw(25) << std::left << (massGasSpecies_outIntegrated[j] + massGasSpeciesTotal[j]) / massTotsolid_initial;

        fmassTotal << std::endl;


        /*fmassTotal << std::setw(25) << std::left << t;
        for (int j = 1; j <= solid_species; j++)
            fmassTotal<< std::setw(25) << std::left << massSolidSpeciesTotal[j];
        for (int j = 1; j <= gas_species; j++)
            fmassTotal << std::setw(25) << std::left << massGasSpeciesTotal[j];*/


        gasFlowRateSelective << std::setw(25) << std::left << t;
        gasFlowTempLG = 0.0;
        for (int i = 0; i < lightGasIndex.size(); i++)
        {
            gasFlowTempLG += massConvection_out[intervalli] * massFraction[intervalli][lightGasIndex[i]];
        }

        gasFlowRateSelective << std::setw(25) << std::left << gasFlowTempLG;
        gasFlowTempTAR = 0.0;
        for (int i = 0; i < tarIndex.size(); i++)
        {
            gasFlowTempTAR += massConvection_out[intervalli] * massFraction[intervalli][tarIndex[i]];
        }
        gasFlowRateSelective << std::setw(25) << std::left << gasFlowTempTAR;


        gasFlowRateSelective << std::endl;

        for (int k = 0; k < output_species_.size(); k++)
        {
            std::string type = output_species_[k];
            int index = thermodynamicsSolidMapXML->IndexOfSpecies(type);
            species[k] << std::setw(25) << std::left << t;
            if (index <= NGas_)
            {
                for (int j = 1; j <= intervalli; j++)
                {
                    species[k] << std::setw(25) << std::left << massFraction[j][index + NSolid_];
                    //species[k] << std::setw(25) << std::left << massFraction_Shellmassbased[j][index + NSolid_];
                }
            }
            else
            {
                for (int j = 1; j <= intervalli; j++)
                {
                    species[k] << std::setw(25) << std::left << massFraction[j][index - NGas_];
                    //species[k] << std::setw(25) << std::left << massFraction_Shellmassbased[j][index - NGas_];
                }
            }
            species[k] << std::endl;
        }


        count_file_ = 0;

    }

    count_ode_video_++;
    count_file_++;

    //massResidualSolid = 0.0;
    //for (int i = 1; i<= intervalli; i++)
    //{
    //    for (int j = 0; j< residualIndex.size(); j++)
    //    {
    //        massResidualSolid += mass[i][residualIndex[j]];
    //    }
    //}
    //yieldChar = massResidualSolid/massTotsolid_initial;
    //
    //gasFlowRateSelective << std::setw(35) << std::left << massResidualSolid;
    //gasFlowRateSelective << std::setw(35) << std::left << yieldChar;
    //
    //gasFlowRateSelective << std::endl;

}

void ChangeDimensionsFunction()
{

    ChangeDimensions(intervalli, gas_species + solid_species, &mass, true);
    ChangeDimensions(intervalli, gas_species + solid_species, &dmass, true);
    ChangeDimensions(intervalli, &T, true);
    ChangeDimensions(intervalli, &dT, true);

    ChangeDimensions(intervalli, gas_species + solid_species, &mole, true);
    ChangeDimensions(solid_species, &massSolid, true);
    ChangeDimensions(solid_species, &moleSolid, true);
    ChangeDimensions(gas_species, &massGas, true);
    ChangeDimensions(gas_species, &moleGas, true);
    ChangeDimensions(intervalli, &massTot_solid, true);
    ChangeDimensions(intervalli, &massTot_gas, true);
    ChangeDimensions(intervalli, &moleTot_solid, true);
    ChangeDimensions(intervalli, &moleTot_gas, true);
    ChangeDimensions(intervalli, gas_species + solid_species, &massFraction, true);
    ChangeDimensions(intervalli, gas_species + solid_species, &moleFraction, true);
    ChangeDimensions(intervalli, gas_species + solid_species, &massFraction_Shellmassbased, true);

    //Properties
    ChangeDimensions(gas_species, &omegaGasFraction, true);
    ChangeDimensions(solid_species, &omegaSolidFraction, true);
    ChangeDimensions(gas_species, &xGasFraction, true);
    ChangeDimensions(solid_species, &xSolidFraction, true);
    ChangeDimensions(intervalli, &MW_gas_mix, true);
    ChangeDimensions(intervalli, &MW_solid_mix, true);
    ChangeDimensions(intervalli, &cPgas_mix, true);
    ChangeDimensions(intervalli, &lambda_gas, true);
    ChangeDimensions(intervalli, &dyn_viscosity, true);
    ChangeDimensions(gas_species, &Diff, true);
    ChangeDimensions(intervalli, gas_species, &Dmix, true);
    ChangeDimensions(intervalli, gas_species, &Diff_Knudsen, true);
    ChangeDimensions(intervalli, gas_species, &Diff_eff, true);
    ChangeDimensions(gas_species, &Hi, true);
    ChangeDimensions(intervalli, gas_species, &Hi_mass, true);
    ChangeDimensions(gas_species + solid_species, &Hi_tot, true);
    ChangeDimensions(gas_species + solid_species, &Si_tot, true);// momentanee
    ChangeDimensions(gas_species + solid_species, &Hi_tot_corr, true);  //momentanee
    ChangeDimensions(intervalli, &Hmix, true);  //momentanee
    ChangeDimensions(intervalli, gas_species + solid_species, &Hi_mass_tot, true);  //momentanee

    ChangeDimensions(intervalli, &C_gas, true);
    ChangeDimensions(intervalli, &rho_gas, true);
    ChangeDimensions(solid_species + gas_species, &cp_tot, true);
    ChangeDimensions(gas_species, &cp_gas_species, true);
    ChangeDimensions(solid_species, &cp_solid_species, true);
    ChangeDimensions(intervalli, &cPsolid_mix, true);
    ChangeDimensions(solid_species, &lambda_solid_species, true);
    ChangeDimensions(intervalli, &lambda_solid, true);
    ChangeDimensions(intervalli, &lambda_effective, true);

    //Kinetics
    ChangeDimensions(intervalli, gas_species + solid_species, &FormationRate, true);
    ChangeDimensions(intervalli, &HeatReaction, true);
    ChangeDimensions(intervalli, solid_reactions, &DHr, true);
    ChangeDimensions(intervalli, gas_species, &Vel_gas, true);
    ChangeDimensions(solid_reactions, &RGG, true);
    ChangeDimensions(intervalli, &sumRgas, true);
    ChangeDimensions(intervalli, &sumRsolid, true);
    ChangeDimensions(intervalli, &sumMoleFormationRates, true);
    //energy balance
    ChangeDimensions(intervalli, gas_species, &JD_out, true);
    ChangeDimensions(intervalli, gas_species, &JD_in, true);
    ChangeDimensions(intervalli, &massConvection_in, true);
    ChangeDimensions(intervalli, &massConvection_out, true);


    ChangeDimensions(intervalli, &JC_out, true);
    ChangeDimensions(intervalli, &JC_in, true);
    ChangeDimensions(intervalli, &JDH_in, true);
    ChangeDimensions(intervalli, &JDH_out, true);
    ChangeDimensions(intervalli, &flux_in, true);
    ChangeDimensions(intervalli, &flux_out, true);


    ChangeDimensions(intervalli * Neq, &residual, true);

    //PRINT

    ChangeDimensions(intervalli, &residual_CELL, true);
    ChangeDimensions(intervalli, &residual_CELLA, true);
    ChangeDimensions(intervalli, &residual_CHAR, true);
    ChangeDimensions(intervalli, &residual_ASH, true);
    ChangeDimensions(intervalli, &residual_TANN, true);
    ChangeDimensions(intervalli, &residual_TGL, true);
    ChangeDimensions(intervalli, &massTANN, true);
    ChangeDimensions(intervalli, &massITANN, true);
    ChangeDimensions(intervalli, &massTGL, true);

    ChangeDimensions(intervalli, &residual_CO, true);
    ChangeDimensions(intervalli, &residual_CO2, true);
    ChangeDimensions(intervalli, &residual_H2, true);
    ChangeDimensions(intervalli, &residual_H2O, true);
    ChangeDimensions(intervalli, &residual_N2, true);

    ChangeDimensions(intervalli, gas_species, &Vel_gas, true);  //momentanee
    ChangeDimensions(intervalli, &GasFormato, true);  //momentanee
    ChangeDimensions(solid_reactions, &RGG, true);

    ChangeDimensions(gas_species, &massFinalGasSpecies, true);
    ChangeDimensions(gas_species, &massFinalGasSpecies_total, true);
    ChangeDimensions(gas_species, &moleFinalGasSpecies_total, true);
    ChangeDimensions(gas_species, &omegaI, true);
    ChangeDimensions(gas_species, &Kc, true);

    ChangeDimensions(gas_species, &massGasSpeciesTotal, true);
    ChangeDimensions(solid_species, &massSolidSpeciesTotal, true);
    ChangeDimensions(gas_species, &massGasSpeciesFlowrate, true);


    ChangeDimensions(intervalli, &massSolid_interval, true);  //momentanee
    //momentanee
    ChangeDimensions(intervalli, &Ugrid, true);
    ChangeDimensions(intervalli, &Vstore, true);
    ChangeDimensions(intervalli, &dmSolid, true);
}

void CreateFractions()
{
    massTot_solid = 0.;
    moleTot_solid = 0.;
    massTot_gas = 0.;
    moleTot_gas = 0.;

    for (int i = 1; i <= intervalli; i++)
    {
        massSolid = 0.;
        moleSolid = 0.;
        massGas = 0.;
        moleGas = 0.;

        for (int j = 1; j <= (gas_species + solid_species); j++)
        {
            if (j <= solid_species)
            {
                massSolid[j] = mass[i][j];
                moleSolid[j] = mole[i][j];
            }

            else
            {
                massGas[j - solid_species] = mass[i][j];
                moleGas[j - solid_species] = mole[i][j];
            }

        }

        massTot_solid[i] = massSolid.SumElements();
        moleTot_solid[i] = moleSolid.SumElements();
        massTot_gas[i] = massGas.SumElements();
        moleTot_gas[i] = moleGas.SumElements();
    }

    for (int i = 1; i <= intervalli; i++)
    {
        for (int j = 1; j <= (gas_species + solid_species); j++)
        {
            massFraction_Shellmassbased[i][j] = mass[i][j] / (massTot_solid[i] + massTot_gas[i]);
            if (j <= solid_species)
            {
                massFraction[i][j] = mass[i][j] / massTot_solid[i];
                moleFraction[i][j] = mole[i][j] / moleTot_solid[i];
            }
            else
            {
                massFraction[i][j] = mass[i][j] / massTot_gas[i];
                moleFraction[i][j] = mole[i][j] / moleTot_gas[i];
            }
        }
    }
}

void SetVolume(const double t)
{
    if (volumeLoss == true)
    {
        /*for(int i = 1; i<= intervalli; i++)
            Vstore[i] = V[i];*/

        /*
         double raggioTot = 0.;
         std::cout << "dmass  " << ele << endl;
         getchar();

         double rate = ele/rhoSolid/(1.-epsi_var[1])/S[intervalli]; //                   1e-05;   //m/s
         double Riniz = 0.009948;
         raggioTot = Riniz + rate*t;



         raggio[1] = raggioTot/intervalli;

         for (int i = 2; i<= intervalli;i++)
             raggio[i] = raggio[i-1]+raggio[1];



         V[1] = 4./3*M_PI*pow(raggio[1],3);

         for (int i = 2; i<= intervalli;i++)
             V[i] = 4./3*M_PI*pow(raggio[i],3)-4./3*M_PI*pow(raggio[i-1],3);


         r[1] = raggio[1]/2.;

         if (intervalli >1)
         {
             for(int i=2; i<=intervalli; i++)
                 r[i] = (raggio[i]+raggio[i-1])/2.;
         }

         for(int i = 1; i<= intervalli; i++)
                 S[i] = 4.*M_PI*pow(raggio[i],2.);


          Vtot = V.SumElements();
          */

        V = 0.;
        S = 0.;
        raggio = 0.;
        r = 0.;
        Vtot = 0.;

        for (int i = 1; i <= intervalli; i++)
            V[i] = massTot_solid[i] / (rhoSolid * (1 - epsi_var[i]));

        Vtot = V[1];
        raggio[1] = pow((3. * V[1] / 4. / PhysicalConstants::pi), 1. / 3.);

        if (intervalli > 1)
        {
            for (int i = 2; i <= intervalli; i++)
            {
                Vtot += V[i];
                raggio[i] = pow(3. * Vtot / 4. / PhysicalConstants::pi, 1. / 3.);
            }
        }
        
        r[1] = raggio[1] / 2.;

        if (intervalli > 1)
        {
            for (int i = 2; i <= intervalli; i++)
                r[i] = (raggio[i] + raggio[i - 1]) / 2.;
        }

        for (int i = 1; i <= intervalli; i++)
            S[i] = 4. * PhysicalConstants::pi * pow(raggio[i], 2.);

        Vtot = V.SumElements();

        for (int i = 1; i <= intervalli; i++)
            Ugrid[i] = (Rstart[i] - raggio[i]) / (t);


        /*  V = 0.;
         S = 0.;
         raggio = 0.;
         r = 0.;
         Vtot = 0.;
        for(int i = 1; i<= intervalli; i++)
             V[i] = massTot_solid[i]/(rhoSolid*(1-epsi_var[i]));

         Vtot = V[1];
         raggio[1] = V[1];

         if (intervalli >1)
         {
             for(int i=2; i<=intervalli; i++)
             {
                 Vtot += V[i];
                 raggio[i] = Vtot;
             }
         }

         r[1] = raggio[1]/2.;

         if (intervalli >1)
         {
             for(int i=2; i<=intervalli; i++)
                 r[i] = (raggio[i]+raggio[i-1])/2.;
         }

         for(int i = 1; i<= intervalli; i++)
                 S[i] = 4.*(raggio[k]-raggio[k-1])+2;

          Vtot = V.SumElements();

         for(int i = 1; i<= intervalli; i++)
             Ugrid[i]= (Rstart[i]-raggio[i])/(t);
          */
          /*
          for (int i = 1; i<= intervalli;i++)
              std::cout << "V["<<i<<"]  " << V[i] << endl;

          getchar();

          for (int i = 1; i<= intervalli;i++)
              std::cout << "S["<<i<<"]  " << S[i] << endl;

          getchar();

          for (int i = 1; i<= intervalli;i++)
              std::cout << "raggio["<<i<<"]  " << raggio[i] << endl;

          getchar();

          for (int i = 1; i<= intervalli;i++)
             std::cout << "r["<<i<<"]  " << r[i] << endl;

          getchar();

          std::cout << "Vtot  " << V.SumElements() << endl;
          getchar();
          */
    }

}

void CalculateProperties()
{
    if (volumeLoss == false)
    {
        for (int i = 1; i <= intervalli; i++)
            epsi_var[i] = 1 - massTot_solid[i] / V[i] / rhoSolid;
    }

    for (int i = 1; i <= intervalli; i++)
    {

        // Proprietà del gas
        for (int j = 1; j <= gas_species; j++)
            omegaGasFraction[j] = massFraction[i][j + solid_species];

        for (int j = 1; j <= solid_species; j++)
            omegaSolidFraction[j] = massFraction[i][j];

        thermodynamicsMapXML->MoleFractions_From_MassFractions(xGasFraction.GetHandle(), MW_gas, omegaGasFraction.GetHandle());
        MW_gas_mix[i] = MW_gas;

        P[i] = massTot_gas[i] / MW_gas_mix[i] / (V[i] * epsi_var[i]) * PhysicalConstants::R_J_kmol * T[i];


        thermodynamicsMapXML->SetPressure(P[i]);
        thermodynamicsMapXML->SetTemperature(T[i]);
        thermodynamicsSolidMapXML->SetTemperature(T[i]);
        thermodynamicsSolidMapXML->SetPressure(P[i]);

        transportMapXML->SetPressure(P[i]);
        transportMapXML->SetTemperature(T[i]);


        thermodynamicsSolidMapXML->cpMolar_Species(cp_tot.GetHandle()); //Cp_solid_species [J/kmol/K];
        for (int j = 1; j <= gas_species; j++)
            cp_gas_species[j] = cp_tot[j] / MW_tot[j + solid_species];  //Cp_solid_species [J/kg/K];

        //cPgas_mix[i]= Dot(cp_gas_species,omegaGasFraction); // cPgas_mix[j/kg/K]

        cPgas_mix[i] = thermodynamicsMapXML->cpMolar_Mixture_From_MoleFractions(xGasFraction.GetHandle());			//[J/Kmol/K]
        cPgas_mix[i] = cPgas_mix[i] / MW_gas_mix[i];

        thermal_conductivity = transportMapXML->ThermalConductivity(xGasFraction.GetHandle());
        lambda_gas[i] = thermal_conductivity;

        viscosity = transportMapXML->DynamicViscosity(xGasFraction.GetHandle());
        dyn_viscosity[i] = viscosity;

        transportMapXML->MassDiffusionCoefficients(Diff.GetHandle(), xGasFraction.GetHandle());

        for (int j = 1; j <= gas_species; j++)
            Dmix[i][j] = Diff[j];

        for (int j = 1; j <= gas_species; j++)
            Diff_Knudsen[i][j] = 9700e-4 * 0.5e-6 * pow((T[i] / MW_tot[j + solid_species]), 0.5);

        for (int j = 1; j <= gas_species; j++)
            //Diff_eff[i][j] = epsi_var[i]/tau*((Dmix[i][j]*Diff_Knudsen[i][j])/(Dmix[i][j]+Diff_Knudsen[i][j]));
            //Diff_eff[i][j] = pow(epsi_var[i],3/2) * Dmix[i][j] // Darcy diffusion
            Diff_eff[i][j] = Dmix[i][j] * 0.25;

        Hmix[i] = thermodynamicsMapXML->hMolar_Mixture_From_MoleFractions(xGasFraction.GetHandle());
        Hmix[i] /= MW_gas_mix[i];
        thermodynamicsMapXML->hMolar_Species(Hi.GetHandle());  // Hi J/kmol

        for (int j = 1; j <= gas_species; j++)
            Hi_mass[i][j] = Hi[j] / MW_tot[j + solid_species];  // Hi_mass   J/kg


        thermodynamicsSolidMapXML->hMolar_Species(Hi_tot.GetHandle());
        thermodynamicsSolidMapXML->sMolar_Species(Si_tot.GetHandle());


        /*std::cout << "temperature  " << T[i] << endl;
        for (int j = 1; j<=Hi_tot.Size();j++)
        {
            std::cout << thermodynamicsSolidMapXML->NamesOfSpecies()[j-1] << "\t" << Hi_tot[j]/1000/4.186 << "\t" << Si_tot[j]/1000/4.186 << endl;

        }

        getchar();*/


        for (int j = 1; j <= Hi_tot.Size(); j++)
        {
            if (j <= solid_species)
                Hi_tot_corr[j] = Hi_tot[j + gas_species];
            else
                Hi_tot_corr[j] = Hi_tot[j - solid_species];
        }

        for (int j = 1; j <= gas_species + solid_species; j++)
            Hi_mass_tot[i][j] = Hi_tot_corr[j];


        C_gas[i] = P[i] / (PhysicalConstants::R_J_kmol * T[i]);  //moleTot_gas[i]/(V[i]*epsi_var[i]);

        rho_gas[i] = P[i] * MW_gas_mix[i] / (PhysicalConstants::R_J_kmol * T[i]);  //moleTot_gas[i]/(V[i]*epsi_var[i]);

        // Proprietà del solido

        for (int j = 1; j <= solid_species; j++)
            cp_solid_species[j] = cp_tot[j + gas_species] / MW_tot[j];  //Cp_solid_species [J/kg/K];

        cPsolid_mix[i] = Dot(cp_solid_species, omegaSolidFraction);


        for (int j = 1; j <= solid_species; j++)
            lambda_solid_species[j] = lambdaS;//lambdaA[j] + lambdaB[j]*T[i];

        lambda_solid[i] = (Dot(lambda_solid_species, omegaSolidFraction));

        lambda_effective[i] = epsi_var[i] * lambda_gas[i] + (1 - epsi_var[i]) * lambda_solid[i];
    }
}

void CalculateKinetics()
{

    DHr = 0.;
    Vel_gas = 0.;
    HeatReactionTot = 0.;

    OpenSMOKE::OpenSMOKEVectorDouble cGas(thermodynamicsSolidMapXML->number_of_gas_species());
    OpenSMOKE::OpenSMOKEVectorDouble RGasG(thermodynamicsSolidMapXML->number_of_gas_species());

    double cTotGas_ = 1.;
    double MWgas = 1.;
    double T_K = 300.;
    double P_Pa = 1.e5;


    for (int i = 1; i <= intervalli; i++)
    {
        // Gas-phase kinetics

        T_K = T[i];
        P_Pa = P[i];

        for (int j = 1; j <= gas_species; j++)
            omegaGasFraction[j] = massFraction[i][j + solid_species];

        for (int j = 1; j <= solid_species; j++)
            omegaSolidFraction[j] = massFraction[i][j];
        thermodynamicsMapXML->MoleFractions_From_MassFractions(xGasFraction.GetHandle(), MWgas, omegaGasFraction.GetHandle());
        cTotGas_ = P_Pa / PhysicalConstants::R_J_kmol / T_K;
        Product(cTotGas_, xGasFraction, &cGas);


        thermodynamicsMapXML->SetTemperature(T_K);
        thermodynamicsMapXML->SetPressure(P_Pa);
        kineticsMapXML->SetPressure(P_Pa);
        kineticsMapXML->SetTemperature(T_K);
        kineticsMapXML->ReactionEnthalpiesAndEntropies();
        kineticsMapXML->ReactionRates(cGas.GetHandle()); //kmol/m3/s
        kineticsMapXML->FormationRates(RGasG.GetHandle());

        HeatReaction[i] = kineticsMapXML->HeatRelease(RGasG.GetHandle()) * V[i] * (epsi_var[i]); // W 
        for (int j = 1; j <= gas_species + solid_species; j++)
        {
            if (j <= solid_species)
                FormationRate[i][j] = 0.;
            else
                FormationRate[i][j] = RGasG[j - solid_species];
        }

        for (unsigned int k = 0; k < thermodynamicsSolidMapXML->number_of_materials(); k++)
        {

            // Calculating kinetics
            {
                /*thermodynamicsMapXML->MoleFractions_From_MassFractions(xGasFraction.GetHandle(), MWgas, omegaGasFraction.GetHandle());
                cTotGas_ = 1000000000.;
                Product(cTotGas_, xGasFraction, &cGas);*/

                double MWsolid = 0.;
                thermodynamicsSolidMapXML->SolidMoleFractions_From_SolidMassFractions(xSolidFraction.GetHandle(), MWsolid, omegaSolidFraction.GetHandle());
                MW_solid_mix[i] = MWsolid;
                OpenSMOKE::OpenSMOKEVectorDouble cSolid(solid_species);
                double cTotSolid_ = moleTot_solid[i] / V[i] / (1 - epsi_var[i]);
                Product(cTotSolid_, xSolidFraction, &cSolid);

                thermodynamicsSolidMapXML->SetTemperature(T_K);
                thermodynamicsSolidMapXML->SetPressure(P_Pa);
                kineticsSolidMapXML[k]->SetPressure(P_Pa);
                kineticsSolidMapXML[k]->SetTemperature(T_K);
                kineticsSolidMapXML[k]->ReactionEnthalpiesAndEntropies();
                kineticsSolidMapXML[k]->ReactionRates(cGas.GetHandle(), cSolid.GetHandle()); //kmol/m3/s

                OpenSMOKE::OpenSMOKEVectorDouble RGas(thermodynamicsSolidMapXML->number_of_gas_species());
                OpenSMOKE::OpenSMOKEVectorDouble RSolid(thermodynamicsSolidMapXML->number_of_solid_species());
                kineticsSolidMapXML[k]->FormationRates(RGas.GetHandle(), RSolid.GetHandle());

                for (int j = 1; j <= gas_species + solid_species; j++)
                {
                    if (j <= solid_species)
                        FormationRate[i][j] += RSolid[j];
                    else
                        FormationRate[i][j] += RGas[j - solid_species];
                }

                sumMoleFormationRates[i] = 0.0;//V[i]*(1-epsi_var[i])*RGas.SumElements();
                sumRgas[i] = 0.0;
                for (int j = 1; j <= gas_species; j++)
                    sumRgas[i] += RGas[j] * thermodynamicsMapXML->MW(j - 1);

                sumRsolid[i] = 0.0;
                for (int j = 1; j <= gas_species; j++)
                    sumRsolid[i] += RSolid[j] * MW_solid[j];


                // std::cout << sumRgas[i] << "\t" << sumRsolid[i] << "\t" << sumRgas[i]+ sumRsolid[i] << endl;
                // getchar();
                 //HeatReaction[i] = kineticsSolidMapXML[k]->HeatRelease(RGas, RSolid)*V[i]*(1-epsi_var[i]); // W 

                HeatReaction[i] += kineticsSolidMapXML[k]->HeatRelease(RGas.GetHandle(), RSolid.GetHandle()) * V[i] * (1 - epsi_var[i]); // W 

                //HeatReaction[i] = 0.0; // W 

                HeatReactionTot += HeatReaction[i];

                for (int j = 1; j <= solid_reactions; j++)
                    RGG[j] = kineticsSolidMapXML[k]->GiveMeReactionRates()[j-1];

                for (int j = 1; j <= solid_reactions; j++)
                    for (int k = 1; k <= gas_species; k++)
                        Vel_gas[i][k] += nu[j][k + solid_species] * RGG[j];

                //HeatReaction[i] = -Dot(RGG,enthalpyReaz)*1e3*V[i]*(1-epsi_var[i]);

            }

        }

        // Calcolo DHr (UNUSED??? AL)
        for (int j = 1; j <= solid_reactions; j++)
        {
            for (int z = 1; z <= solid_species + gas_species; z++)
            {
                DHr[i][j] += (Hi_mass_tot[i][z] * nu[j][z]);
            }

        }
    }

}

void CreateMaterialBalances()
{
    // Boundary volumes
    // r=0  ---> dm/dr=0    

    massConvection_in = 0.;
    massConvection_out = 0.;
    JD_in = 0.;
    JD_out = 0.;
    GasFormato = 0.;
    massDiffusion_out = 0.;
    dmSolid = 0.;
    ele = 0.;

    /*for (int i = 1; i<=dmSolid.Size();i++)
        std::cout << "dmSolid[" << i << "]  " << dmSolid[i] << endl;
    getchar();*/
    if (intervalli > 1)
    {
        for (int i = 1; i <= (solid_species + gas_species); i++)
        {
            if (i <= solid_species)
            {
                dmass[1][i] = FormationRate[1][i] * MW_tot[i] * (V[1] * (1 - epsi_var[1]));
                dmSolid[1] += FormationRate[1][i] * MW_tot[i] * (V[1] * (1 - epsi_var[1]));
            }
            else
            {
                massConvection_in[1] = 0.;
                massConvection_out[1] = -Da / dyn_viscosity[1] * rho_gas[1] * S[1] * (P[2] - P[1]) / (r[2] - r[1]);

                JD_in[1][i - solid_species] = 0.;
                JD_out[1][i - solid_species] = -Diff_eff[1][i - solid_species] * rho_gas[1] * S[1] * (massFraction[2][i] - massFraction[1][i]) / (r[2] - r[1]);

                dmass[1][i] = JD_in[1][i - solid_species] - JD_out[1][i - solid_species] + (FormationRate[1][i] * MW_tot[i] * V[1] * (1 - epsi_var[1])) - massConvection_out[1] * massFraction[1][i];

            }

            /*for (int j = 1; j<=solid_species;j++)
             std::cout << "dmSolid[1][" << j << "]  " << dmass[1][j] << endl;
            getchar();*/
        }

        for (int i = 2; i <= intervalli - 1; i++)
        {
            for (int j = 1; j <= (solid_species + gas_species); j++)
            {
                if (j <= solid_species)
                {
                    dmass[i][j] = FormationRate[i][j] * MW_tot[j] * V[i] * (1 - epsi_var[i]);
                    dmSolid[i] += dmass[i][j];
                }

                else
                {
                    massConvection_in[i] = massConvection_out[i - 1];
                    massConvection_out[i] = -Da / dyn_viscosity[i] * rho_gas[i] * S[i] * (P[i + 1] - P[i]) / (r[i + 1] - r[i]);

                    JD_in[i][j - solid_species] = JD_out[i - 1][j - solid_species];
                    JD_out[i][j - solid_species] = -Diff_eff[i][j - solid_species] * rho_gas[i] * S[i] * (massFraction[i + 1][j] - massFraction[i][j]) / (r[i + 1] - r[i]);

                    dmass[i][j] = JD_in[i][j - solid_species] - JD_out[i][j - solid_species] + (FormationRate[i][j] * MW_tot[j] * V[i] * (1 - epsi_var[i])) + massConvection_in[i] * massFraction[i - 1][j] - massConvection_out[i] * massFraction[i][j];

                }
            }
        }

    }

    // Calculate interfacial properties
    for (int i = 1; i <= gas_species; i++) {
        omegaI[i] = (Kc[i] * omegaIn_gas[i] + Diff_eff[intervalli][i] * massFraction[intervalli][i] / (raggio[intervalli] - r[intervalli])) /
            (Kc[i] + Diff_eff[intervalli][i] / (raggio[intervalli] - r[intervalli]));
    }
    double MW_gasI = thermodynamicsMapXML->MolecularWeight_From_MassFractions(omegaI.GetHandle());
    double rhoGasI = P[intervalli] * MW_gasI / PhysicalConstants::R_J_kmol * T[intervalli]; // approx

    // r=R  ---> m=mS  
    for (int i = 1; i <= (solid_species + gas_species); i++)
    {
        if (i <= solid_species)
        {
            dmass[intervalli][i] = FormationRate[intervalli][i] * MW_tot[i] * V[intervalli] * (1 - epsi_var[intervalli]);
            dmSolid[intervalli] += dmass[intervalli][i];
        }
        else
        {
            if (intervalli > 1) {
                massConvection_in[intervalli] = massConvection_out[intervalli - 1];
                JD_in[intervalli][i - solid_species] = JD_out[intervalli - 1][i - solid_species];
            }
            
            massConvection_out[intervalli] = -Da / dyn_viscosity[intervalli] * rho_gas[intervalli] * S[intervalli] * (Psolid - P[intervalli]) / (raggio[intervalli] - r[intervalli]);

            JD_out[intervalli][i - solid_species] = -Kc[i - solid_species] * S[intervalli]* (rhoGas * omegaIn_gas[i - solid_species] - rhoGasI * omegaI[i - solid_species]);
            
            dmass[intervalli][i] = JD_in[intervalli][i - solid_species] - JD_out[intervalli][i - solid_species] + (FormationRate[intervalli][i] * MW_tot[i] * V[intervalli] * (1 - epsi_var[intervalli])) + (massConvection_in[intervalli] * massFraction[intervalli - 1][i] - massConvection_out[intervalli] * massFraction[intervalli][i]);

            massDiffusion_out += JD_out[intervalli][i - solid_species];

        }

    } 

    
    for (int i = 1; i <= intervalli; i++)
    {
        for (int j = 1; j <= solid_species; j++)
        {
            ele += dmass[i][j];
        }
    }

    // ele = dmSolid.SumElements();


     /*for (int i = 1; i<=dmSolid.Size();i++)
         std::cout << "dmSolid[" << i << "]  " << dmSolid[i] << endl;
     getchar();*/

     /* std::cout << "dmass2  " << ele << "\t" << dmSolid.SumElements() << endl;
      getchar();*/
}

void CreateEnergyBalances()
{
    if (energyBalance == true)
    {
        flux_in = 0.;
        flux_out = 0.;

        JDH_in = 0.;
        JDH_out = 0.;

        // Boundary volumes
        // r=0  ---> dT/dr=0  

        for (int i = 1; i <= (gas_species); i++)
        {
            // deriva da bilancio materiale
            flux_in[1] = 0.;


            //diffusività in and out
            JDH_in[1] = 0.;
            //JDH_out[1] += Hi_mass[1][i]*(JD_out[1][i]+massConvection_out[1]*massFraction[1][i+solid_species]);

        }

        JC_in[1] = 0.;
        //JC_out[1] = -lambda_solid[1]*(1-epsi_var[1])*S[1]*(T[2]-T[1])/(r[2]-r[1])-lambda_gas[1]*S[1]*epsi_var[1]*(T[2]-T[1])/(r[2]-r[1]);
        JC_out[1] = -lambda_effective[1] * S[1] * (T[2] - T[1]) / (r[2] - r[1]);


        dT[1] = (-(flux_in[1] - flux_out[1]) + HeatReaction[1] + JC_in[1] - JC_out[1] + JDH_in[1] - JDH_out[1]/*+PhysicalConstants::R_J_kmol*T[1]*sumMoleFormationRates[1]*/)
            / ((massTot_solid[1] * cPsolid_mix[1]) + (massTot_gas[1] * cPgas_mix[1])/*-PhysicalConstants::R_J_kmol*C_gas[1]*V[1]*epsi_var[1]*/);





        for (int i = 2; i <= intervalli - 1; i++)
        {

            for (int j = 1; j <= (gas_species); j++)
            {
                flux_in[i] += Hi_mass[i][j] * (JD_in[i][j] + massConvection_in[i] * massFraction[i - 1][j + solid_species]);
                //flux_out[i] += Hi_mass[i][j]*(JD_out[i][j]+massConvection_out[i]*massFraction[i][j+solid_species]);


                JDH_in[i] += Hi_mass[i - 1][j] * (JD_in[i][j] + massConvection_in[i] * massFraction[i - 1][j + solid_species]);
                //JDH_out[i] += Hi_mass[i][j]*(JD_out[i][j]+massConvection_out[i]*massFraction[i][j+solid_species]);


            }

            //JC_in[i] = -lambda_solid[i-1]*S[i-1]*(1-epsi_var[i])*(T[i]-T[i-1])/(r[i]-r[i-1])-lambda_gas[i-1]*epsi_var[i]*S[i-1]*(T[i]-T[i-1])/(r[i]-r[i-1]);
            //JC_out[i] = -lambda_solid[i]*(1-epsi_var[i])*S[i]*(T[i+1]-T[i])/(r[i+1]-r[i])-lambda_gas[i]*S[i]*epsi_var[i]*(T[i+1]-T[i])/(r[i+1]-r[i]);
            JC_in[i] = -lambda_effective[i - 1] * S[i - 1] * (T[i] - T[i - 1]) / (r[i] - r[i - 1]);
            JC_out[i] = -lambda_effective[i] * S[i] * (T[i + 1] - T[i]) / (r[i + 1] - r[i]);

            dT[i] = (-(flux_in[i] - flux_out[i]) + HeatReaction[i] + JC_in[i] - JC_out[i] + JDH_in[i] - JDH_out[i]/*+PhysicalConstants::R_J_kmol*T[i]*sumMoleFormationRates[i]*/)
                / ((massTot_solid[i] * cPsolid_mix[i]) + (massTot_gas[i] * cPgas_mix[i])/*- PhysicalConstants::R_J_kmol*C_gas[i]*V[i]*epsi_var[i]*/);

        }


        // r=R  ---> T=TS

        for (int i = 1; i <= (gas_species); i++)
        {
            //flux_out[intervalli] += Hi_mass[intervalli][i]*(JD_out[intervalli][i]+massConvection_out[intervalli]*massFraction[intervalli][i+solid_species]);
            flux_in[intervalli] += Hi_mass[intervalli][i] * (JD_in[intervalli][i] + massConvection_in[intervalli] * massFraction[intervalli - 1][i + solid_species]);

            JDH_in[intervalli] += Hi_mass[intervalli - 1][i] * (JD_in[intervalli][i] + massConvection_in[intervalli] * massFraction[intervalli - 1][i + solid_species]);
            //JDH_out[intervalli] +=  Hi_mass[intervalli][i]*(JD_out[intervalli][i]+massConvection_out[intervalli]*massFraction[intervalli][i+solid_species]);

        }


        //double Tsurf;
        Tsurf_solid = (hext * Tbulk + lambda_effective[intervalli] * T[intervalli] / (raggio[intervalli] - r[intervalli])) /
            (hext + lambda_effective[intervalli] / (raggio[intervalli] - r[intervalli]));


        //JC_in[intervalli] = -lambda_solid[intervalli-1]*(1-epsi_var[i]_var[intervalli-1])*S[intervalli-1]*(T[intervalli]-T[intervalli-1])/(r[intervalli]-r[intervalli-1])-lambda_gas[intervalli-1]*S[intervalli-1]*epsi_var[i]_var[intervalli-1]*(T[intervalli]-T[intervalli-1])/(r[intervalli]-r[intervalli-1]);
        //JC_out[intervalli] = -lambda_solid[intervalli]*(1-epsi_var[i]_var[intervalli])*S[intervalli]*(Tsup_solid-T[intervalli])/(raggio[intervalli]-r[intervalli])-lambda_gas[intervalli]*S[intervalli]*epsi_var[i]_var[intervalli]*(Tsup_solid-T[intervalli])/(raggio[intervalli]-r[intervalli]);
        JC_in[intervalli] = -lambda_effective[intervalli - 1] * S[intervalli - 1] * (T[intervalli] - T[intervalli - 1]) / (r[intervalli] - r[intervalli - 1]);
        JC_out[intervalli] = -lambda_effective[intervalli] * S[intervalli] * (Tsurf_solid - T[intervalli]) / (raggio[intervalli] - r[intervalli]);

        //dT[intervalli] = (-(flux_in[intervalli]-flux_out[intervalli])+HeatReaction[intervalli]+JC_in[intervalli]-JC_out[intervalli]+JDH_in[intervalli]-JDH_out[intervalli]/*+PhysicalConstants::R_J_kmol*T[intervalli]*sumMoleFormationRates[intervalli]*/)
          //          /((massTot_solid[intervalli]*cPsolid_mix[intervalli]) + (massTot_gas[intervalli]*cPgas_mix[intervalli])/*- PhysicalConstants::R_J_kmol*C_gas[intervalli]*V[intervalli]*epsi_var[intervalli]*/);

        dT[intervalli] = (-(flux_in[intervalli] - flux_out[intervalli]) + HeatReaction[intervalli] + JC_in[intervalli] - JC_out[intervalli] + JDH_in[intervalli] - JDH_out[intervalli]) / ((massTot_solid[intervalli] * cPsolid_mix[intervalli]) + (massTot_gas[intervalli] * cPgas_mix[intervalli]));


        //dT[intervalli] = 0.;

    }
    else
    {
        for (int i = 1; i <= intervalli; i++)
            dT[i] = Heat_Rates;

    }
}