clc, clear all, close all %#ok<*CLALL>

%
% Post processing full DLC 1.2 to get fatigue loads and other metrics for
% the whole life time.
% The inputs are folders with .res files from FAST
% Witlis style. The weibull, seeds and yaw misallignment is taken to
% account and the results are presented also relative to the first case.
% DELs are calcualted for Wohler exponent of 4 and 10. The pitch metrics
% include pitch travel, rate and acceleration which can be evaluates as
% accumulated values, min,max,mean and std. 
% Fatigue calculations can be done for life extension purposes. Damage
% values for the baseline (1st) case can be defined (for now one value for
% all but can be easily modified) and subsequently calculate damage margins
% and lifetime extension based on simplified MIner's linear damage
% hypothesis. Only statistical values are stored and no time-series
%
% Dependencies: rainflow function
% 
% Vasilis Pettas 2.2018 SWE 
%
% update 4.12.2018 VPE: - parametrized safety factors (SF) where 2 was left hardcopied in the code
%                       - Fixed seeds for non-DEL channels 
%                       - In LT cacluations for non-DEL (Energy, pitch travel etc.) the total years are considered variable  
%                       - Pitchrate is now converted to degrees!!!


%% Input files

ResDLC_folders = {
    
%     {'D:\Tasks\Torque2018\Full_DLC_results\Base_DLC12_mYaw';   % all files including all speeds. Make sure all are correct before!!!!
%      'D:\Tasks\Torque2018\Full_DLC_results\Base_DLC12';        % Order mYaw,NoYaw, pYaw !!! if no Yaw put the same file 3 times
%      'D:\Tasks\Torque2018\Full_DLC_results\Base_DLC12_pYaw' }
%
%         {'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE2_Base_DLC12_mYaw'
%         'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE2_Base_DLC12'
%         'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE2_Base_DLC12_pYaw'}
%         %
%         {'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE2_Boost15_Mg_DLC12_mYaw'
%         'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE2_Boost15_Mg_DLC12'
%         'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE2_Boost15_Mg_DLC12_pYaw'}
%         %
%         {'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE2_Boost15_Omega_DLC12_mYaw'
%         'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE2_Boost15_Omega_DLC12'
%         'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE2_Boost15_Omega_DLC12_pYaw'}
%         %
%         {'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE2_Boost15_MgOmega_DLC12_mYaw'
%         'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE2_Boost15_MgOmega_DLC12'
%         'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE2_Boost15_MgOmega_DLC12_pYaw'}

%       {'D:\Tasks\Torque2020\DLC12\Base_DLC12_mYaw'
%         'D:\Tasks\Torque2020\DLC12\Base_DLC12'
%         'D:\Tasks\Torque2020\DLC12\Base_DLC12_pYaw'}
%         %
%         {'D:\Tasks\Torque2020\DLC12\Boost10_Mg_DLC12_mYaw'
%         'D:\Tasks\Torque2020\DLC12\Boost10_Mg_DLC12'
%         'D:\Tasks\Torque2020\DLC12\Boost10_Mg_DLC12_pYaw'}
%         %
%         {'D:\Tasks\Torque2020\DLC12\Boost10_Omega_DLC12_mYaw'
%         'D:\Tasks\Torque2020\DLC12\Boost10_Omega_DLC12'
%         'D:\Tasks\Torque2020\DLC12\Boost10_Omega_DLC12_pYaw'}
%         %
%         {'D:\Tasks\Torque2020\DLC12\Boost10_MgOmega_DLC12_mYaw'
%         'D:\Tasks\Torque2020\DLC12\Boost10_MgOmega_DLC12'
%         'D:\Tasks\Torque2020\DLC12\Boost10_MgOmega_DLC12_pYaw'}

        {'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE3_Base_DLC12'
        'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE3_Base_DLC12_seed2'
        'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE3_Base_DLC12_seed3'}
        %
        {'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE3_Boost15_Mg_DLC12_n'
        'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE3_Boost15_Mg_DLC12_seed2_n'
        'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE3_Boost15_Mg_DLC12_seed3_n'}
        %
%          {'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE3_Boost15_Mg_DLC12'
%         'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE3_Boost15_Mg_DLC12_seed2'
%         'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE3_Boost15_Mg_DLC12_seed3'}
        %
        {'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE3_Boost15_Omega_DLC12_n'
        'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE3_Boost15_Omega_DLC12_seed2_n'
        'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE3_Boost15_Omega_DLC12_seed3_n'}
        %
%         {'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE3_Boost15_Omega_DLC12'
%         'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE3_Boost15_Omega_DLC12_seed2'
%         'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE3_Boost15_Omega_DLC12_seed3'}
        %
        {'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE3_Boost15_MgOmega_DLC12'
        'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE3_Boost15_MgOmega_DLC12_seed2'
        'D:\Tasks\Torque2020\DLC12\DTU10MW_FBSWE3_Boost15_MgOmega_DLC12_seed3'}

};

%% Parameter Definition

Speed_bins = [4:2:24];        %#ok<*NBRAK>   % Mean value of the wind speed bins
NBins      = length(Speed_bins);
P_wind     = [0.1100,0.1412,0.1512,0.1427,0.1214,0.09437,0.06749,0.04463,0.02739,0.01563,0.0083]; % Weibull for Class Ia has to have the same size as the speed bins!
P_yaw      = [0 1 0]; % If yaw cases are included give the weighting (plus yaw, no yaw, minus yaw) otherwise (0,1,0)
T_sim      = 3600;            % [s] Time duration of each simulation. All should have the same length!
years_LT   = 20;
T_life     = years_LT*365*24*60*60; % [s] Total life time of wind turbine
N_seeds    = 1;               % Number of seeds for each case (yaw and no yaw cases should have the same number of seeds)
P_seed     = 1/N_seeds;       % Probability of each seed (considered equally possible...) has to be same for all cases
N_ref      = 1e7 ;            % Reference cycles for fatigue calculations
SF         = 1 ;              % Safety factor for DEL calculation
Dbase      = 1;               % Damage from baseline. It is considered one for lifetime extension calculations
T_baseline = (15/20)*T_life  ;   % Time of operation with baseine setup
T_red      = T_life-T_baseline  ;  % Time of operation with load reduction method
timeinterval = 4801:292801;    % Time to be considerd in simulations (indices)

%% Plotting options
plotperbinABS = 0; %absoute values per bin
plotperbinREL = 1; %relative values per bin
plotLTABS     = 0; %absoute lifetime values
plotLTREL     = 0; %relative  lifetime values
plotEXT       = 0; % lifetime extension plot   

%% Save options

save_DELall = 0; %save all load channel information on probability and rainflow counting for further fatigue calculations
save_name   = 'D:\Tasks\Torque2020\DLC12\DEL_all.mat' ;

%% Channels requested     
ChannReq = {
    % Loads
            'RootMxb1'
            'RootMyb1'            
            'RootMzb1'
%             'RootMxb2'
%             'RootMyb2'                        
%             'RootMzb2'  
            'TwrBsMxt'
            'TwrBsMyt'
            'TwrBsMzt'
            'LSShftTq' 
%            'YawBrMyn'                
            'YawBrMxp'  % Roll TT
            'YawBrMyp'  % Pitch TT
            'YawBrMzn'  % Yaw TT/yaw bearing
%             'LSSTipMys' % Non rotating low speed shaft tip (hub) around yaxis HUB
%             'LSSTipMzs' % Non rotating low speed shaft tip (hub) around yaxisHUB           
    % Stats       
            'BldPitch1'
%             'BldPitch2'
%             'BldPitch3'            
%             'Wind1VelX'
%             'RotSpeed'
%             'TTDspFA'
%             'M_g'
            'GenTq'
%             'Azimuth'
            'GenPwr';
            'GenSpeed'                                                      
            };

%%  Channels for each plot       
DEL10Chan = {
            'RootMzb1';
            'RootMyb1';
            'RootMxb1';
            'RootMzb2';
            'RootMyb2';
            'RootMxb2';
            'RootMzb3';
            'RootMyb3';
            'RootMxb3';
            };
         
DEL4Chan = {
            'TwrBsMxt'
            'TwrBsMyt'
            'TwrBsMzt'            
            'LSShftTq' 
            'YawBrMzn'  % Yaw TT/yaw bearing
            'YawBrMyn'             
            'YawBrMxp'  % Roll TT
            'YawBrMyp'  % Pitch TT
            'LSSTipMys' % Non rotating low speed shaft tip (hub) around yaxis HUB
            'LSSTipMzs' % Non rotating low speed shaft tip (hub) around yaxisHUB
            'YawBrMx2'
            };       

%%
for iFold = 1:length(ResDLC_folders)   
    
    mYawResFold  = ResDLC_folders{iFold,1}{1};
    NoYawResFold = ResDLC_folders{iFold,1}{2};
    pYawResFold  = ResDLC_folders{iFold,1}{3};    
    
    % Get files names for each case   
    mYawFilesN = dir(fullfile(mYawResFold,'*results*.mat'));
    NoYawFilesN = dir(fullfile(NoYawResFold,'*results*.mat'));
    pYawFilesN = dir(fullfile(pYawResFold,'*results*.mat'));    
    
    for i = 1:length(Speed_bins)
        mYawFiles{i}  = [mYawResFold '\' mYawFilesN(i).name];
        NoYawFiles{i} = [NoYawResFold '\' NoYawFilesN(i).name];
        pYawFiles{i}  = [pYawResFold '\' pYawFilesN(i).name];          %#ok<*SAGROW>
    end
    
    Allfiles = [mYawFiles';NoYawFiles';pYawFiles']; % the order is exact!!!
        
    % Get for each case all the files and add a weighting for each case based on the speed the yaw and the seed probabilities!
    for  iFile = 1:length(Allfiles)
        load (Allfiles{iFile})
        
        % Find probability and bin for each speed
        Vsim = Parameter.TurbSim.URef; % Wind speed of Simulation 
        if iFile <= NBins
            P_yawsim  = P_yaw(1);
        elseif iFile > 2*NBins
            P_yawsim  = P_yaw(3);
        else
            P_yawsim  = P_yaw(2);  
        end    
        P_windsim = P_wind((Speed_bins==Vsim));
        P_seedsim = P_seed;        
%         Inp_Sim(iFile,:) = [Vsim P_yawsim P_windsim P_seedsim]; % gather all inputs for each simulation
        
        % Get time series and rainflow       
        %% Load Data for each file 
        try
            DataAll   = logsout.getElement('OutData').Values.Data(timeinterval,:);
            TSobjY    = logsout.getElement('y').Values;
            PitchRateIn  = rad2deg(TSobjY.theta_dot.Data(timeinterval,:));
            TSobju    = logsout.getElement('u').Values;
            GenTq     = TSobju.M_g.Data(timeinterval,:);

        catch
            TSobj    = logsout.getElement('OutData').Values;   
            DataAll  = TSobj.OutData.Data(timeinterval,:);
            TSobjY   = logsout.getElement('y').Values;
            PitchRateIn  = rad2deg(TSobjY.y.theta_dot.Data(timeinterval,:));
        end
        Chanels  = Parameter.FASTInput.OutList  ;
        
        %% Get rainflow and values required from timeseries
        for iChan = 1:length (ChannReq)
            clear indfind indmat
            Chan_Nam = ChannReq{iChan};
            indfind = strfind (Chanels, ChannReq{iChan});

            for i = 1:length (indfind)
                if indfind{i}==1
                    indmat(i) = 1; %#ok<*SAGROW>
                else
                    indmat(i) = 0;
                end
            end
            
            %% Time series
            Ind = find(indmat==1);            
            TS     = DataAll(:,Ind);
            TS_STD = std(TS);
            % Pitch Rate           
            PitchRate_STD_in = std (PitchRateIn);   
            PitchRate_Cum_in = sum(abs(PitchRateIn));              
            % Pitch travel
            if strcmp(Chan_Nam,'BldPitch1')
                Pitchtrav = sum(abs(diff(TS))) ;  
                Pitchtravel(iFile,:) = {[Vsim] [P_yawsim] [P_windsim] [P_seedsim] [Pitchtrav]};
                Pitch_acc = diff(PitchRateIn)*(tout(2)-tout(1));
                Pitch_acc_std_calc = std(Pitch_acc);
                Pitch_acc_std(iFile,:)={[Vsim] [P_yawsim] [P_windsim] [P_seedsim] [Pitch_acc_std_calc]};
            end
            % Energy
            if strcmp(Chan_Nam,'GenPwr')            
                EnergyP = trapz(tout(timeinterval),TS);              
                EnergyProd (iFile,:) = {[Vsim] [P_yawsim] [P_windsim] [P_seedsim] [EnergyP]};
                Power_STD (iFile,:)  = {[Vsim] [P_yawsim] [P_windsim] [P_seedsim] [TS_STD]};
            end
            
            %% Rainflow calculation 
            [ext, exttime] = sig2ext(double(1e3*DataAll(:,Ind)),tout(2)-tout(1));
            rf = rainflow(ext,exttime);
            cyclesAmplitude = rf(1,:);
            cyclesNumber    = rf(3,:);
            Y_DEL4          = (cyclesNumber/N_ref).*cyclesAmplitude.^4;   % creating the common expression (Ni/Nref)*Si^m!!
            Y_DEL10         = (cyclesNumber/N_ref).*cyclesAmplitude.^10;  % creating the common expression (Ni/Nref)*Si^m!!
            
            eval([Chan_Nam '(iFile,:)={[Vsim] [P_yawsim] [P_windsim] [P_seedsim] Y_DEL4 Y_DEL10 [TS_STD]};']);  % gather all inputs for each simulation
            Pitchrate_STD  (iFile,:) = {[Vsim] [P_yawsim] [P_windsim] [P_seedsim] [PitchRate_STD_in]};
            Pitchrate_CUM  (iFile,:) = {[Vsim] [P_yawsim] [P_windsim] [P_seedsim] [PitchRate_Cum_in]};
            
            if save_DELall==1
                DEL_info_out.Loads{iFold}{iChan}(iFile,:) = {[Vsim] [P_yawsim] [P_windsim] [P_seedsim] Y_DEL4 Y_DEL10 [TS_STD]};
                DEL_info_out.Pitchrate_STD{iFold}(iFile,:) = Pitchrate_STD  (iFile,:) ;
                DEL_info_out.Pitchrate_CUM{iFold}(iFile,:) = Pitchrate_CUM  (iFile,:) ;
                if strcmp(Chan_Nam,'BldPitch1')
                    DEL_info_out.Pitch_acc_std{iFold}(iFile,:) = Pitch_acc_std(iFile,:);
                    DEL_info_out.Pitchtravel{iFold}(iFile,:) = Pitchtravel(iFile,:);
                end
                if strcmp(Chan_Nam,'GenPwr')
                    DEL_info_out.EnergyProd{iFold}(iFile,:) = EnergyProd (iFile,:) ;
                    DEL_info_out.Power_STD{iFold}(iFile,:) = Power_STD (iFile,:) ;
                end
                DEL_info_out.Names{iChan} = ChannReq{iChan} ;
                DEL_info_out.Names{iChan} = ChannReq{iChan} ;
                
            end
            clear rf ext exttime Y_DEL4 Y_DEL10 cyclesAmplitude cyclesNumber
        end
    end
    
    if save_DELall==1
        save(save_name,'DEL_info_out')
    end
    
    %% Calculation of quantities per bin
    for iChan = 1:length (ChannReq)          % loop over channels
        LoadChan = eval(ChannReq{iChan});
        speedInp = cell2mat(LoadChan(:,1));
        
        for iBin = 1:NBins            
            Speed     = Speed_bins(iBin);            
            Speed_Ind1 = speedInp==Speed;
            Speed_Ind  = find(Speed_Ind1==1); % Get the indices with this speed 
            
            Y_DEL4_curBin  = LoadChan(Speed_Ind,5);
            Y_DEL10_curBin = LoadChan(Speed_Ind,6);
            Pyaw_curBin    = LoadChan(Speed_Ind,2);
            Pseed_curBin    = LoadChan(Speed_Ind,4);
            
            for iSpeedInd = 1:length(Speed_Ind)
                DEL4_ind_bin(iSpeedInd)  = Pseed_curBin{iSpeedInd}*Pyaw_curBin{iSpeedInd}*sum((T_life/T_sim)*Y_DEL4_curBin{iSpeedInd}) ;
                DEL10_ind_bin(iSpeedInd) = Pseed_curBin{iSpeedInd}*Pyaw_curBin{iSpeedInd}*sum((T_life/T_sim)*Y_DEL10_curBin{iSpeedInd}) ;           
            end            
            DEL4_curBin    = SF*(sum(DEL4_ind_bin))^(1/4);
            DEL10_curBin   = SF*(sum(DEL10_ind_bin))^(1/10);            
            eval([ChannReq{iChan} '_DEL4_perBin' '(iBin,iFold)= DEL4_curBin;']);
            eval([ChannReq{iChan} '_DEL10_perBin' '(iBin,iFold)= DEL10_curBin;']);
            
            
            if strcmp(ChannReq{iChan},'GenPwr')
               EnergyProd_curBin = cell2mat(EnergyProd(Speed_Ind,5));
               EnergyProd_perBin(iBin,iFold) = sum(EnergyProd_curBin.*cell2mat(Pyaw_curBin).*cell2mat(EnergyProd(Speed_Ind,4)));
               Power_STD_curBin = cell2mat(Power_STD(Speed_Ind,5));
               Power_STD_perBin(iBin,iFold) = sum(Power_STD_curBin.*cell2mat(Pyaw_curBin).*cell2mat(Power_STD(Speed_Ind,4)));
            end
            
            if strcmp(ChannReq{iChan},'BldPitch1')               
               Pitchrate_STD_curBin = cell2mat(Pitchrate_STD(Speed_Ind,5));
               Pitchrate_CUM_curBin = cell2mat(Pitchrate_CUM(Speed_Ind,5));
               Pitchtravel_curBin   = cell2mat(Pitchtravel(Speed_Ind,5));
               Pitch_acc_std_curBin = cell2mat(Pitch_acc_std(Speed_Ind,5));
               Pitchrate_STD_perBin(iBin,iFold) = sum(Pitchrate_STD_curBin.*cell2mat(Pyaw_curBin).*cell2mat(Pitchrate_STD(Speed_Ind,4))); 
               Pitchrate_CUM_perBin(iBin,iFold) = sum(Pitchrate_CUM_curBin.*cell2mat(Pyaw_curBin).*cell2mat(Pitchrate_CUM(Speed_Ind,4))); 
               Pitchtravel_perBin(iBin,iFold)   = sum(Pitchtravel_curBin.*cell2mat(Pyaw_curBin).*cell2mat(Pitchtravel(Speed_Ind,4))); 
               Pitch_acc_std_perBin(iBin,iFold) = sum(Pitch_acc_std_curBin.*cell2mat(Pyaw_curBin).*cell2mat(Pitch_acc_std(Speed_Ind,4)));
               
            end
            
            if strcmp(ChannReq{iChan},'GenTq')               
               GenTq_STD_curBin             = cell2mat(GenTq(Speed_Ind,7));
               GenTq_STD_perBin(iBin,iFold) = sum(GenTq_STD_curBin.*cell2mat(Pyaw_curBin).*cell2mat(GenTq(Speed_Ind,4)));                                
            end
            
            if strcmp(ChannReq{iChan},'GenSpeed')               
               GenSpeed_STD_curBin             = cell2mat(GenSpeed(Speed_Ind,7));
               GenSpeed_STD_perBin(iBin,iFold) = sum(GenSpeed_STD_curBin.*cell2mat(Pyaw_curBin).*cell2mat(GenSpeed(Speed_Ind,4)));                                
            end
            
            clear Speed_Ind DEL4_ind_bin DEL10_ind_bin DEL4_curBin DEL10_curBin EnergyProd_curBin iSpeedInd Pyaw_curBin
        end
                
        clear LoadChan Speed_Ind  
    
    %% Lifetime Calculations
        LoadChanLT = eval(ChannReq{iChan});
        
        for iCase = 1:length(Allfiles)       
            P_iCase      = LoadChanLT{iCase,2}*LoadChanLT{iCase,3}*LoadChanLT{iCase,4}; % probability calculation yaw x weibull x seed  
            DEL4_LT_int(iCase)  = P_iCase*sum( (T_life/T_sim)*LoadChanLT{(iCase),5} ); % lifetime intermediate variables only time,probability and intermediate varaible with woehler and rainflow results
            DEL10_LT_int(iCase) = P_iCase*sum( (T_life/T_sim)*LoadChanLT{(iCase),6} ); 
            clear P_iCase 
        end
        
        DEL4_curCase  = SF*(sum(DEL4_LT_int))^(1/4); % adding all weigted DELs for life time with m=4
        DEL10_curCase = SF*(sum(DEL10_LT_int))^(1/10);  % adding all weigted DELs for life time with m=10
        eval([ChannReq{iChan} '_DEL4_LT(iFold)= DEL4_curCase;']);       
        eval([ChannReq{iChan} '_DEL10_LT(iFold)= DEL10_curCase;']);  
        
        if strcmp(ChannReq{iChan},'GenPwr')
            EnergyProd_LT(iFold) = years_LT*sum( cell2mat(EnergyProd(:,2)).*cell2mat(EnergyProd(:,3)).*cell2mat(EnergyProd(:,4)).*cell2mat( EnergyProd(:,5)) );
            Power_STD_LT(iFold) =  sum( cell2mat(Power_STD(:,2)).*cell2mat(Power_STD(:,3)).*cell2mat(Power_STD(:,4)).*cell2mat( Power_STD(:,5)) );
        end
        
        if strcmp(ChannReq{iChan},'BldPitch1')
            Pitchtravel_LT(iFold)   =years_LT*sum( cell2mat(Pitchtravel(:,2)).*cell2mat(Pitchtravel(:,3)).*cell2mat(Pitchtravel(:,4)).*cell2mat( Pitchtravel(:,5)) );
            Pitchrate_STD_LT(iFold) = sum( cell2mat(Pitchrate_STD(:,2)).*cell2mat(Pitchrate_STD(:,3)).*cell2mat(Pitchrate_STD(:,4)).*cell2mat( Pitchrate_STD(:,5)) );
            Pitchrate_CUM_LT(iFold) = years_LT*sum( cell2mat(Pitchrate_CUM(:,2)).*cell2mat(Pitchrate_CUM(:,3)).*cell2mat(Pitchrate_CUM(:,4)).*cell2mat( Pitchrate_CUM(:,5)) );
            Pitch_acc_std_LT(iFold) = sum( cell2mat(Pitch_acc_std(:,2)).*cell2mat(Pitch_acc_std(:,3)).*cell2mat(Pitch_acc_std(:,4)).*cell2mat( Pitch_acc_std(:,5)) );
        end
        
        if strcmp(ChannReq{iChan},'GenTq')
            GenTq_STD_LT(iFold)   =  sum( cell2mat(GenTq(:,2)).*cell2mat(GenTq(:,3)).*cell2mat(GenTq(:,4)).*cell2mat( GenTq(:,7)) );
        end
        
        if strcmp(ChannReq{iChan},'GenSpeed')
            GenSpeed_STD_LT(iFold)   =  sum( cell2mat(GenSpeed(:,2)).*cell2mat(GenSpeed(:,3)).*cell2mat(GenSpeed(:,4)).*cell2mat( GenSpeed(:,7)) );
        end
            
        clear LoadChanLT DEL4_LT DEL10_LT DEL10_curCase DEL4_curCase DEL4_LT_int DEL10_LT_int
    end
end

if save_DELall==1
    save(save_name,'DEL_info_out')
end

%% Life time Extension
if length(ResDLC_folders)>=2  
    for iChanEXT = 1:length(ChannReq)        
        ChanEXT  = ChannReq{iChanEXT};
        
        Dnew_4      = (T_baseline/T_life) + (T_red/T_life).* ( eval([ChanEXT '_DEL4_LT']) ./ eval([ChanEXT '_DEL4_LT(1)']) ).^4 ;    % damage after 20 years
        Dnew_10     = (T_baseline/T_life) + (T_red/T_life).* ( eval([ChanEXT '_DEL10_LT']) ./ eval([ChanEXT '_DEL10_LT(1)']) ).^10 ; % damage after 20 years
        LText_cur4  = years_LT.*(1-Dnew_4)./( eval([ChanEXT '_DEL4_LT']) ./ eval([ChanEXT '_DEL4_LT(1)']) ).^4; %#ok<*NASGU>
        LText_cur10 = years_LT.*(1-Dnew_10)./( eval([ChanEXT '_DEL10_LT']) ./ eval([ChanEXT '_DEL10_LT(1)']) ).^10;  
        eval([ChannReq{iChanEXT} '_DEL4_Ext= LText_cur4;'])
        eval([ChannReq{iChanEXT} '_DEL10_Ext= LText_cur10;'])
        eval([ChannReq{iChanEXT} '_Dnew4_Ext= Dnew_4;'])
        eval([ChannReq{iChanEXT} '_Dnew10_Ext= Dnew_10;'])
        
        if strcmp(ChannReq{iChanEXT},'GenPwr')
            EnergyProd_LT_20yrComb =  (T_baseline/T_life)*EnergyProd_LT(1)+ (T_red/T_life)*EnergyProd_LT;
            Power_STD_LT_20yrComb  = (T_baseline/T_life)*Power_STD_LT(1)+ (T_red/T_life)*Power_STD_LT;
        end
        
        if strcmp(ChannReq{iChanEXT},'BldPitch1')
            Pitchtravel_LT_20yrComb   = (T_baseline/T_life)*Pitchtravel_LT(1)+ (T_red/T_life)*Pitchtravel_LT;
            Pitchrate_STD_LT_20yrComb = (T_baseline/T_life)*Pitchrate_STD_LT(1)+ (T_red/T_life)*Pitchrate_STD_LT;     
            Pitch_acc_std_LT_20yrComb = (T_baseline/T_life)*Pitch_acc_std_LT(1)+ (T_red/T_life)*Pitch_acc_std_LT;            
        end
        
        if strcmp(ChannReq{iChanEXT},'GenTq')
            GenTq_STD_LT_20yrComb    =  (T_baseline/T_life)*GenTq_STD_LT(1)+ (T_red/T_life)*GenTq_STD_LT;
        end
        
        if strcmp(ChannReq{iChanEXT},'GenSpeed')
            GenSpeed_STD_LT_20yrComb = (T_baseline/T_life)*GenTq_STD_LT(1)+ (T_red/T_life)*GenTq_STD_LT;
        end
        
        clear Dnew ChanEXT LText_cur4 LText_cur10
    end
    
    
end

%% Plotting
% Per bin

plotLabelBar= {'BladeEW' 'BladeFW' 'BladeTOR' 'TwrBSS' 'TwrBFA' 'TwrBTOR' 'LSShaftTor' 'TTroll' 'TTpitch' 'TTyaw'};

if plotperbinABS==1
    figure
    count=0;
    for i = 1:length(ChannReq)
        
        curChan = ChannReq{i};
        if find(strcmp(DEL10Chan,curChan)==1)
            count=count+1;
            subplot(4,3,count);
            set(gca, 'ColorOrder', [0.0 0.45 0.74; 1 0 0; 0 0.5 0;0.923 0.694 0.1255], 'NextPlot', 'replacechildren');
            plot (Speed_bins,eval([curChan '_DEL10_perBin']),'x-','Linewidth',2 )
%             xlabel ('Wind Speed Bins [m/s]', 'FontSize', 14)
            ylabel ([ plotLabelBar(count) '[kNm]'], 'Interpreter','latex','FontSize', 14)
            xlim ([4 24])
%             grid on
%             title ([curChan 'DEL Wexp 10 per bin'], 'FontSize', 14)
        end
        
        if find(strcmp(DEL4Chan,curChan)==1)
            count=count+1;
            subplot(4,3,count);
            set(gca, 'ColorOrder', [0.0 0.45 0.74; 1 0 0; 0 0.5 0;0.923 0.694 0.1255], 'NextPlot', 'replacechildren');
            plot (Speed_bins,eval([curChan '_DEL4_perBin']),'x-','Linewidth',2 )
%             xlabel ('Wind Speed Bins [m/s]', 'FontSize', 14)
            ylabel ([ plotLabelBar(count) '[kNm]'],  'Interpreter','latex','FontSize', 14)
            xlim ([4 24])            
%             grid on
%             title ([curChan 'DEL Wexp 4 per bin'], 'FontSize', 14)
        end
        if find(strcmp('BldPitch1',curChan)==1)
            count=count+1;
            subplot(4,3,count);
            %pitchrate std
            
            set(gca, 'ColorOrder', [0.0 0.45 0.74; 1 0 0; 0 0.5 0;0.923 0.694 0.1255], 'NextPlot', 'replacechildren');
            plot(Speed_bins,Pitchrate_STD_perBin,'x-','Linewidth',2)
%             xlabel ('Wind Speed Bins [m/s]', 'FontSize', 14)
            ylabel (['$\dot{\theta}$[deg/s]'], 'Interpreter','latex', 'FontSize', 14)
            xlim ([4 24])              
%             grid on
%             title (['Absolute STD of $\dot{\theta}$ per bin'], 'Interpreter','latex', 'FontSize', 14)
         
            %pitch travel
            count=count+1;
            subplot(4,3,count);
            set(gca, 'ColorOrder', [0.0 0.45 0.74; 1 0 0; 0 0.5 0;0.923 0.694 0.1255], 'NextPlot', 'replacechildren');
            plot(Speed_bins,Pitchtravel_perBin,'x-','Linewidth',2)
%             xlabel ('Wind Speed Bins [m/s]', 'FontSize', 14)
            ylabel (['Pitch travel[deg]'], 'Interpreter','latex', 'FontSize', 14)
            xlim ([4 24])              
%             grid on
%             title ([' Absolute Total Pitch Travel per bin'], 'FontSize', 14)                     
        end
        if find(strcmp('GenPwr',curChan)==1)      
%             % energy
%             figure
%             set(gca, 'ColorOrder', [0.0 0.45 0.74; 1 0 0; 0 0.5 0;0.923 0.694 0.1255], 'NextPlot', 'replacechildren');
%             plot(Speed_bins,EnergyProd_perBin,'x-','Linewidth',2)
%             xlabel ('Wind Speed Bins [m/s]', 'FontSize', 14)
%             ylabel (['Energy Poduced [kWh]'], 'FontSize', 14)
%             grid on
%             title (['Energy Production per bin'], 'FontSize', 14)
        end
         
        if find(strcmp('GenTq',curChan)==1)   
            %generator speed std
            figure
            set(gca, 'ColorOrder', [0.0 0.45 0.74; 1 0 0; 0 0.5 0;0.923 0.694 0.1255], 'NextPlot', 'replacechildren');
            plot(Speed_bins,GenSpeed_STD_perBin,'x-','Linewidth',2)
            xlabel ('Wind Speed Bins [m/s]', 'FontSize', 14)
            ylabel (['Generator speed STD [deg/s]'], 'FontSize', 14)
            grid on
            title (['Absolute STD of generator speed per bin'], 'FontSize', 14)    
            xlim ([4 24])              
        end          
    end  
end

if plotperbinREL==1
    if length(ResDLC_folders)>=2
        for i = 1:length(ChannReq)
            
            curChan = ChannReq{i};
            if find(strcmp(DEL10Chan,curChan)==1)
                Plot_Val=(100*(eval(['bsxfun(@rdivide,' curChan '_DEL10_perBin,' curChan '_DEL10_perBin(:,1))'])-1));
                figure
                plot (Speed_bins, Plot_Val(:,2:end),'x-','Linewidth',2 )
                xlabel ('Wind Speed Bins [m/s]', 'FontSize', 14)
                ylabel (['Relative difference ' curChan ' [%]'], 'FontSize', 14)
                grid on
                title ([curChan 'DEL Wexp 10'], 'FontSize', 14)
            end
            
            if find(strcmp(DEL4Chan,curChan)==1)
                Plot_Val=(100*(eval(['bsxfun(@rdivide,' curChan '_DEL4_perBin,' curChan '_DEL4_perBin(:,1))'])-1));
                figure
                plot (Speed_bins, Plot_Val(:,2:end),'x-','Linewidth',2 )
                xlabel ('Wind Speed Bins [m/s]', 'FontSize', 14)
                ylabel (['Relative difference ' curChan ' [%]'], 'FontSize', 14)
                grid on
                title ([curChan 'DEL Wexp 4'], 'FontSize', 14)
            end
        end
        
        % energy
        figure
        Plot_Val=(100*(bsxfun(@rdivide,EnergyProd_perBin,EnergyProd_perBin(:,1))-1));
        plot(Speed_bins,Plot_Val(:,2:end),'x-','Linewidth',2)
        xlabel ('Wind Speed Bins [m/s]', 'FontSize', 14)
        ylabel (['Relative Energy Poduction [%]'], 'FontSize', 14)
        grid on
        title (['Relative Energy Production per bin'], 'FontSize', 14)
        
        %pitchrate std
        figure
        Plot_Val=(100*(bsxfun(@rdivide,Pitchrate_STD_perBin,Pitchrate_STD_perBin(:,1))-1));
        plot(Speed_bins,Plot_Val(:,2:end),'x-','Linewidth',2)
        xlabel ('Wind Speed Bins [m/s]', 'FontSize', 14)
        ylabel (['Relative Pitch rate STD [%]'], 'FontSize', 14)
        grid on
        title (['Relative STD of $\dot{\theta}$ per bin'], 'Interpreter','latex', 'FontSize', 14)
        
        %generator speed std
        figure
        Plot_Val=(100*(bsxfun(@rdivide,GenSpeed_STD_perBin,GenSpeed_STD_perBin(:,1))-1));
        plot(Speed_bins,Plot_Val(:,2:end),'x-','Linewidth',2)
        xlabel ('Wind Speed Bins [m/s]', 'FontSize', 14)
        ylabel (['Relative Generator Speed STD [%]'], 'FontSize', 14)
        grid on
        title (['Relative STD of $\dot{\theta}$ per bin'], 'Interpreter','latex', 'FontSize', 14)
        
        %pitch travel
        figure
        Plot_Val=(100*(bsxfun(@rdivide,Pitchtravel_perBin,Pitchtravel_perBin(:,1))-1));
        plot(Speed_bins,Plot_Val(:,2:end),'x-','Linewidth',2)
        xlabel ('Wind Speed Bins [m/s]', 'FontSize', 14)
        ylabel (['Pitch total traveled distance [deg]'], 'FontSize', 14)
        grid on
        title ([' Absolute Total Pitch Travel per bin'],'FontSize', 14)
    end
end

% Lifetime
if plotLTABS ==1
    Total_fat = [];
    Tick_Str  = {};
    count     = 0;
    for i = 1:length(ChannReq)
        curChan = ChannReq{i};
        
        if find(strcmp(DEL10Chan,curChan)==1)
            count = count+1;
            Total_fat = [Total_fat, eval([curChan '_DEL10_LT'])];
            Tick_Str{count,1}  =curChan ;
        elseif find(strcmp(DEL4Chan,curChan)==1)
            count = count+1;
            Total_fat = [Total_fat, eval([curChan '_DEL4_LT'])]; %#ok<*AGROW>
            Tick_Str{count,1}  =curChan ;
        end
    end
    if length(ResDLC_folders)>=2
        for i = 1:length(ResDLC_folders)
            Plot_comb(i,:) = Total_fat(i:length(ResDLC_folders):end) ;
            Plot_metrics (i,:) = [EnergyProd_LT(i) Pitchtravel_LT(i) Pitchrate_STD_LT(i) GenTq_STD_LT(i)];
        end
    else
        Plot_comb = Total_fat;
    end
    figure
    bar ([1:length(Tick_Str)],Plot_comb','grouped' )
    %     xlabel ('Wind Speed Bins [m/s]', 'FontSize', 14)
    ylabel (['DEL Lifetime [MNm]'], 'FontSize', 14)
    set(gca, 'XTickLabel', Tick_Str)
    grid on
    title (['Fatigue lifetime '], 'FontSize', 14)
    
    % energy, pitch rate, pitch travel, omega std
    %     figure
    %     bar([1:4],Plot_metrics,'grouped')
    %     set(gca, 'XTickLabel', {'Energy' 'Pitchtravel' 'STD Pitch Rate' 'STD generator Torque'},'XTickLabelRotation',45,'YGrid','on')
    %     ylabel (['Lifetime metrics '], 'FontSize', 14)
    %      xlim ([0.5 4.5])
    %     grid on
    %     title (['Relative Energy Production per bin'], 'FontSize', 14)
end

if plotLTREL ==1
    Total_fat = [];
    Tick_Str  = {};
    count     = 0;
    if length(ResDLC_folders)>=2
        for i = 1:length(ChannReq)
            curChan = ChannReq{i};
            if find(strcmp(DEL10Chan,curChan)==1)
                count = count+1;
                Total_fat = [Total_fat, eval([curChan '_DEL10_LT'])];
                Tick_Str{count,1}  =curChan ;
            elseif find(strcmp(DEL4Chan,curChan)==1)
                count = count+1;
                Total_fat = [Total_fat, eval([curChan '_DEL4_LT'])]; %#ok<*AGROW>
                Tick_Str{count,1}  =curChan ;
            end
        end
        if length(ResDLC_folders)>=2
            for i = 1:length(ResDLC_folders)
                Plot_comb(i,:)=Total_fat(i:length(ResDLC_folders):end) ;
                Plot_metrics1 (i,:) = [Pitchtravel_LT(i) Pitchrate_STD_LT(i) Pitchrate_CUM_LT(i)  Pitch_acc_std_LT(i)];
                Plot_metrics2 (i,:) = [EnergyProd_LT(i)  GenTq_STD_LT(i) GenSpeed_STD_LT(i) Power_STD_LT(i) ];
            end
        else
            Plot_comb = Total_fat;
        end
        plot_val =  (100*(bsxfun(@rdivide,Plot_comb,Plot_comb(1,:))-1));
        figure
        hb=bar ([1:length(Tick_Str)],plot_val(2:end,:)','grouped' );
        %    xlabel ('Wind Speed Bins [m/s]', 'FontSize', 14)
        ylabel (['Relative DEL Lifetime '], 'FontSize', 14)
        xlim ([0.5 length(Tick_Str)+0.5])
        set(gca, 'XTickLabel', plotLabelBar,'XTickLabelRotation',45,'YGrid','on', 'FontSize', 14,'YMinorGrid','on')
        title (['Fatigue lifetime Relative Comparison '], 'FontSize', 14)
        set(hb(1), 'FaceColor',[1 0 0])
        set(hb(2), 'FaceColor',[0 0.5 0])
        set(hb(3), 'FaceColor',[0.923 0.694 0.1255])
        
        
        % energy
        figure
        Plot_val2=(100*(bsxfun(@rdivide,Plot_metrics1,Plot_metrics1(1,:))-1));
        hb=bar([1:4],Plot_val2(2:end,:)','grouped');
        set(gca, 'XTickLabel', {'Pitchtravel' 'STD Pitch Rate' 'Cummulative Pitch rate' 'STD Pitch Acc'},'XTickLabelRotation',45,'YGrid','on')
        ylabel (['Lifetime metrics 100(IPC/Base-1) [%]'], 'FontSize', 14)
        xlim ([0.5 4.5])
        title (['Relative metrics Lifetime compared to baseline'], 'FontSize', 14)
        set(hb(1), 'FaceColor',[1 0 0])
        set(hb(2), 'FaceColor',[0 0.5 0])
        set(hb(3), 'FaceColor',[0.923 0.694 0.1255])       
        
        figure
        Plot_val3=(100*(bsxfun(@rdivide,Plot_metrics2,Plot_metrics2(1,:))-1));
        hb=bar([1:4],Plot_val3(2:end,:)','grouped');
        set(gca, 'XTickLabel', {'Energy' 'STD Gen Torque' 'STD Gen Speed' 'STD Power'}, 'FontSize', 14,'XTickLabelRotation',45,'YGrid','on')
        ylabel (['Lifetime metrics 100(IPC/Base-1) [%]'], 'FontSize', 14)
        xlim ([0.5 4.5])
        title (['Relative metrics Lifetime compared to baseline'], 'FontSize', 14)
        set(hb(1), 'FaceColor',[1 0 0])
        set(hb(2), 'FaceColor',[0 0.5 0])
        set(hb(3), 'FaceColor',[0.923 0.694 0.1255])        
        
        Tick_Str{end+1}='Pitchtravel';
        Tick_Str{end+1}='STD_Pitch_Rate';
        Tick_Str{end+1}='CUM_Pitch_Rate';        
        Tick_Str{end+1}='STD_Pitch_Acc';        
        Tick_Str{end+1}='Energy';
        Tick_Str{end+1}='STD_generator_Torque';
        Tick_Str{end+1}='STD_generator_speed';
        Tick_Str{end+1}='STD_Power';        
        LT_rel_table =  array2table( [plot_val(2:end,:),Plot_val2(2:end,:),Plot_val3(2:end,:)]...
            ,'VariableNames',Tick_Str);
        LT_rel_table %#ok<NOPTS>
    end
end

if plotEXT == 1
    Total_Ext = [];
    Tick_StrExt = {};
    Dnew_Ext  = [];
    count     = 0;
    if length(ResDLC_folders)>=2
        for i = 1:length(ChannReq)
            curChan = ChannReq{i};
            if find(strcmp(DEL10Chan,curChan)==1)
                count = count+1;
                Total_Ext = [Total_Ext, eval([curChan '_DEL10_Ext'])];
                Dnew_Ext = [Dnew_Ext, eval([curChan '_Dnew10_Ext'])];
                Tick_StrExt{count,1}  =curChan ;
            elseif find(strcmp(DEL4Chan,curChan)==1)
                count = count+1;
                Total_Ext = [Total_Ext, eval([curChan '_DEL4_Ext'])]; %#ok<*AGROW>
                Dnew_Ext = [Dnew_Ext, eval([curChan '_Dnew4_Ext'])];                
                Tick_StrExt{count,1}  =curChan ;
            end
        end
        
        for i = 1:length(ResDLC_folders)
            Plot_Ext(i,:)=Total_Ext(i:length(ResDLC_folders):end) ;
            Plot_Dnew_Ext(i,:)=Dnew_Ext(i:length(ResDLC_folders):end) ;            
            Plot_metrics1_ext (i,:) = [Pitchtravel_LT_20yrComb(i) Pitchrate_STD_LT_20yrComb(i) Pitch_acc_std_LT_20yrComb(i)];
            Plot_metrics2_ext (i,:) = [EnergyProd_LT_20yrComb(i)  GenTq_STD_LT_20yrComb(i) GenSpeed_STD_LT_20yrComb(i) Power_STD_LT_20yrComb(i)];
        end
    end
    
    
    figure
    hb=bar ([1:length(Tick_StrExt)],Plot_Ext(2:end,:)','grouped' );
    %     xlabel ('Wind Speed Bins [m/s]', 'FontSize', 14)
    ylabel (['Lifetime Extension  in years'], 'FontSize', 14)
    set(gca, 'XTickLabel', plotLabelBar,'XTickLabelRotation',45,'YGrid','on', 'FontSize', 14,'YMinorGrid','on')
    title (['Fatigue lifetime Extension'], 'FontSize', 14)
    xlim ([0.5 length(Tick_StrExt)+0.5])
    ylim ([-1 26.5])    

    set(hb(1), 'FaceColor',[1 0 0])
    set(hb(2), 'FaceColor',[0 0.5 0])
    set(hb(3), 'FaceColor',[0.923 0.694 0.1255])
    breakyaxis([3.5 19.5])
        legend({'IPC 1P' 'IPC 1P+2P+3P' 'IPC all'})

    
    figure
    hb=bar ([1:length(Tick_StrExt)],Plot_Dnew_Ext(2:end,:)','grouped' );
    %     xlabel ('Wind Speed Bins [m/s]', 'FontSize', 14)
    ylabel (['Damage after nominal 20 yr lifetime ' ], 'FontSize', 14)
    set(gca, 'XTickLabel', plotLabelBar,'XTickLabelRotation',45,'YGrid','on', 'FontSize', 14,'YMinorGrid','on')
    title (['Damage after 20 years with ' num2str(T_baseline/T_life) '% Baseline and ' num2str(1-(T_baseline/T_life)) ' % IPC'  ], 'FontSize', 14)
    xlim ([0.5 length(Tick_StrExt)+0.5])
    ylim ([0.6 1.1])
    set(hb(1), 'FaceColor',[1 0 0])
    set(hb(2), 'FaceColor',[0 0.5 0])
    set(hb(3), 'FaceColor',[0.923 0.694 0.1255])
    
    figure
    Plot_val_ext=(100*(bsxfun(@rdivide,Plot_metrics1_ext,Plot_metrics1_ext(1,:))-1));
    hb=bar([1:3],Plot_val_ext(2:end,:)','grouped');
    set(gca, 'XTickLabel', {'Pitchtravel' 'STD $\overdot{\theta}$' 'STD Pitch Acc' },'XTickLabelRotation',45,'YGrid','on')
    ylabel (['Relative to baseline [%]'], 'FontSize', 14)
    xlim ([0.5 3.5])
    title (['Relative values 20 years with ' num2str(T_baseline/T_life) '% Baseline and ' num2str(1-(T_baseline/T_life)) ' % IPC'], 'FontSize', 14)
    set(hb(1), 'FaceColor',[1 0 0])
    set(hb(2), 'FaceColor',[0 0.5 0])
    set(hb(3), 'FaceColor',[0.923 0.694 0.1255])    

    figure
    Plot_val2_ext=(100*(bsxfun(@rdivide,Plot_metrics2_ext,Plot_metrics2_ext(1,:))-1));
    hb=bar([1:4],Plot_val2_ext(2:end,:)','grouped');
    set(gca, 'XTickLabel', {'Energy' 'STD Gen Torque' 'STD Gen Speed' 'STD Power'}, 'FontSize', 14,'XTickLabelRotation',45,'YGrid','on')
    ylabel (['Relative to baseline [%]'], 'FontSize', 14)
    xlim ([0.5 4.5])
    title (['Relative values 20 years with ' num2str(T_baseline/T_life) '% Baseline and ' num2str(1-(T_baseline/T_life)) ' % IPC'], 'FontSize', 14)
    set(hb(1), 'FaceColor',[1 0 0])
    set(hb(2), 'FaceColor',[0 0.5 0])
    set(hb(3), 'FaceColor',[0.923 0.694 0.1255])        
    
    Ext_rel_table =  array2table( [Plot_Ext(2:end,:)]...
        ,'VariableNames',Tick_StrExt);
    Ext_rel_table %#ok<NOPTS>
    
    Ext_Dnew_table =  array2table( [Plot_Dnew_Ext(2:end,:)]...
        ,'VariableNames',Tick_StrExt);
    Ext_Dnew_table %#ok<NOPTS>
end

        
