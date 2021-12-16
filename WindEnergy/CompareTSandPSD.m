clc,clear all, close all %#ok<CLALL>


% Script for plotting and comparing TS, PSD and statistics for different
% simulations. The input is results file formatted in the Witlis style for
% FAST simulations. The first case us be used as baseline for comparison
% with the other simulations. Channels are named according to Fast
% nomeclature. PLots can be relative or absolute and except DELs values for
% pitch and generator metrics can be also requested. PLots for now are
% selected by commenting the relevant sections. The results that can be
% obtained include time series and frequency domain response as well as
% tthe requested statistics.
% 
% 
% Dependencies: Rainflow stuff
%
% Vasilis Pettas 1.2018 SWE
%
% Update 4.12.2018 VPE: - parametrized safety factors (SF) where 2 was left hardcopied in the code
%                       - Pitchrate is now converted to degrees!!!
% Update 12.1.2019 VPE: Include post processing of FASTexe simulaitons as
%                       they are provided by the new witlis


%%
Uref = 16;


%% 
% Resfiles = {
%  'D:\Tasks\Torque2018\DynResponse\FAST_BladeSpeed_ID_DLC12\FAST_BladeSpeed_ID_DLC12_16_results.mat'
% 'D:\Tasks\Torque2018\IPCresults\16d0_test_IPC_P_13_results.mat'
% 'D:\Tasks\Torque2018\DynResponse\FUCK_me_test\DTU10MW_FBSWE2_FAST_IPC_16d0testfull0d0003_Td0_results.mat'
% 'D:\Tasks\Torque2018\DynResponse\DTU10MW_FBSWE2_FAST_IPC_fullTI_test\DTU10MW_FBSWE2_FAST_IPC_18d0testfull0d0002_Td8e-05_results.mat'
% 

% 'D:\Tasks\Torque2018\DynResponse\FAST_BladeSpeed_ID_DLC12\FAST_BladeSpeed_ID_DLC12_14_results.mat'
% 'D:\Tasks\Torque2018\DynResponse\DTU10MW_FBSWE2_FAST_IPC_Transition_Test2\DTU10MW_FBSWE2_FAST_IPC_16d0_TestScaling_results.mat'
% };
%% Full DLC12 results with the final conroller
Resfiles = {

['D:\code_witlis_FAST\witlis3\Simulation\FAST\Senvion5MW_OpenFAST_TurbTest1_1\Senvion5MW_OpenFAST_TurbTest1_1_Uref_' num2str(Uref) '_results.mat']
['D:\code_witlis_FAST\witlis3\Simulation\FAST\Senvion5MW_OpenFAST_TurbTest1_2\Senvion5MW_OpenFAST_TurbTest1_2_Uref_' num2str(Uref) '_results.mat']
['D:\code_witlis_FAST\witlis3\Simulation\FAST\Senvion5MW_OpenFAST_TurbTest1_3\Senvion5MW_OpenFAST_TurbTest1_3_Uref_' num2str(Uref) '_results.mat']
['D:\code_witlis_FAST\witlis3\Simulation\FAST\Senvion5MW_OpenFAST_TurbTest1_4\Senvion5MW_OpenFAST_TurbTest1_4_Uref_' num2str(Uref) '_results.mat']
['D:\code_witlis_FAST\witlis3\Simulation\FAST\Senvion5MW_OpenFAST_TurbTest1_5\Senvion5MW_OpenFAST_TurbTest1_5_Uref_' num2str(Uref) '_results.mat']
['D:\code_witlis_FAST\witlis3\Simulation\FAST\Senvion5MW_OpenFAST_TurbTest1_6\Senvion5MW_OpenFAST_TurbTest1_6_Uref_' num2str(Uref) '_results.mat']
% 'D:\code_witlis_FAST\witlis3\Simulation\FAST\IEA37_35MW_Steps_all\IEA37_35MW_Steps_all_URef_4_results.mat'
% ['D:\code_witlis_FAST\witlis3\Simulation\FAST\IEA37_35MW_DLC12_1\IEA37_35MW_DLC12_1_URef_' num2str(Uref) '_results.mat']
% ['D:\code_witlis_FAST\witlis3\Simulation\FAST\IEA37_35MW_DLC12_2_pyaw\IEA37_35MW_DLC12_2_pyaw_URef_' num2str(Uref) '_results.mat' ]
% ['D:\code_witlis_FAST\witlis3\Simulation\FAST\IEA37_35MW_DLC12_3_myaw\IEA37_35MW_DLC12_3_myaw_URef_' num2str(Uref) '_results.mat']


};

%% Channels
ChannReq = {
            'RootMzb1';
            'RootMyb1';
%               'RootMyb2';
%              'RootMyb3';                        
            'RootMxb1';
            'TwrBsMxt'
            'TwrBsMyt'
            'TwrBsMzt'         
            'BldPitch1';
% 'BldPitch2';
            'Wind1VelX';
            'RotSpeed';
            'LSShftTq';
            'TTDspFA';
%               'M_g'
            'GenTq'
            'TTDspSS'            
%             'TipClrnc1'
%       'TwrClrnc1'
% 'Azimuth'
% 'BldPitch2'
%  'BldPitch3'
              'GenPwr'
             'GenSpeed'
             'RotThrust'
%              'HydroFxi'
%              'HydroFyi'
%              'HydroFzi'
%              'HydroMxi'
%              'HydroMyi'
%              'HydroMzi'

'YawBrMxp' %Roll TT
'YawBrMyp' %Pitch TT
'YawBrMzp' %Yaw TT/yaw bearing
'OoPDefl1'
'IPDefl1'
%         'TwrClrnc1'
%  'LSSTipMys'% Non rotating low speed shaft tip (hub) around yaxis HUB
% 'LSSTipMzs'% Non rotating low speed shaft tip (hub) around yaxisHUB
% 'YawBrMx2'z
% 'YawBrMyn'                                                        
            };


%% Parameters   

VarName = { ['V = ' num2str(Uref) ' Baseline'];
             ['V = ' num2str(Uref) ' IPC'];            
%             ['V = ' num2str(Uref) ' TSR 7.6'];
%             ['V = ' num2str(Uref) ' TSR 7.8'];
            };
timeinterval = 3001:33001;%4801:292801;%3001:33001;%84001;%100001;%45166;%%4801:292801; % 146401;%  52801;%    %% SET TIME period to be considered has to be same for all simulations
SF           = 1;   %safety factor for fatigue loads
Nref_basis   = 1e7; %reference cycles
FastExe      = 2 ;  %1 with witlis3 Fastexe orifinal postprocessing/ 2 with VP vhange in processing file to output per file results

%%  Channels for each plot       
MeanChan = {
            'GenPwr';
            };

DEL10Chan = {
            'RootMzb1';
            'RootMyb1';
            'RootMxb1';
            };
         
DEL4Chan = {
            'TwrBsMxt'
            'TwrBsMyt'
            'TwrBsMzt'            
            'LSShftTq'; 
            'YawBrMzn' %Yaw TT
            'YawBrMxp' %Roll TT
            'YawBrMyp' %Pitch TT    
'LSSTipMys'
'LSSTipMzs'
'YawBrMxn'
'YawBrMyn'  
            };

TSChan =  {
            'GenPwr';
            'RootMzb1';
            'RootMyb1';
            'RootMyb2';
            'RootMyb3';            
            'RootMxb1';
            'TwrBsMxt'
            'TwrBsMzt'            
            'GenTq'
            'GenSpeed'
            'TipClrnc1'
            'TwrClrnc1'
             'M_g'
            'TwrBsMyt'
            'BldPitch1';
            'BldPitch2';
            'BldPitch3';
            'Wind1VelX';
            'RotSpeed';
            'LSShftTq';
            'TTDspFA';
            'TTDspSS'
            'YawBrMxp'
            'RotThrust'
            'HydroFxi'
            'HydroFyi'
            'HydroFzi'
            'HydroMxi'
            'HydroMyi'
            'HydroMzi'
            'OoPDefl1'
            'IPDefl1'
            };
        
PSDchan =  {
            'GenPwr';
            'RootMzb1';
            'RootMyb1';
            'RootMxb1';
            'TwrBsMxt'
            'TwrBsMyt'
            'TwrBsMzt'            
            'BldPitch1';
            'BldPitch2';
            'BldPitch3';            
            'Wind1VelX';
            'RotSpeed';
            'LSShftTq';
            'TTDspFA';
            'GenSpeed'
            'YawBrMxp'
            'YawBrMyp'
            'GenTq'
            'HydroFxi'
            'HydroFyi'
            'HydroFzi'
             'HydroMxi'
             'HydroMyi'
             'HydroMzi'
            };        
        
%% Get timeseries and create frequency resposne
 
for NamInd=1:length(Resfiles)
    load(Resfiles{NamInd});
    if FastExe ==0
        try
            DataAll  = logsout.getElement('OutData').Values.Data(timeinterval,:);
            TSobjY    = logsout.getElement('y').Values;
            PitchRateIn  = TSobjY.theta_dot.Data(timeinterval,:);
            TSobju    = logsout.getElement('u').Values;
            GenTq     = TSobju.M_g.Data(timeinterval,:);

        catch
            TSobj    = logsout.getElement('OutData').Values;   
            DataAll  = TSobj.OutData.Data(timeinterval,:);
            TSobjY    = logsout.getElement('y').Values;
            PitchRateIn  = TSobjY.y.theta_dot.Data(timeinterval,:);
        end
        Chanels  = Parameter.FASTInput.OutList  ;
    elseif FastExe ==1
        for i = 1:length(TimeResults.Data)
            DataAll(i,:) = TimeResults.Data{:,i};         
        end
        DataAll = DataAll';
        Chanels = PostProcessingConfig.Channels  (:,2);
        tout    = TimeResults.Time{1, 1}  ;
    elseif FastExe ==2
        DataAll  = FASTData(timeinterval,:);
        Chanels =  OutList;
        tout     = FASTData(:,1);
    end
    
    
    %% rotor effective speed     
    try
        Wind(:,NamInd) = logsout.get('d').Values.v_0.data(timeinterval);
    catch
    end
    %% Get timeseries and PSD from FAST results        
    for j =1:length (ChannReq)
        clear indfind indmat
        indfind = strfind (Chanels, ChannReq{j});

        for i = 1:length (indfind)
            if indfind{i}==1
                indmat(i) = 1; %#ok<*SAGROW>
            else
                indmat(i) = 0;
            end
        end
        
        %% Time series
        Ind = find(indmat==1,1);
        Tsname = [ChannReq{j} '_TS'];
        eval([Tsname '(:,NamInd)= DataAll(:,Ind);']);
        Time_TS(:,NamInd) =tout;
        if FastExe == 0
            PitchRate(:,NamInd) = PitchRateIn;
            PitchRete_STD (:,NamInd) = std (PitchRateIn);
            PitchAcc_std (:,NamInd)  = std(diff(PitchRateIn)*(tout(2)-tout(1)) );
        else
            if strcmp(ChannReq{j},'BldPitch1')
                PitchRate(:,NamInd) = diff(eval([Tsname '(:,NamInd)']))*(tout(2)-tout(1));
%                 PitchRate(:,NamInd) = [PitchRate(1,NamInd);PitchRate(:,NamInd)];
                PitchRete_STD (:,NamInd)  = std (PitchRate(:,NamInd));
                PitchAcc_std (:,NamInd)   = std(diff( PitchRate(:,NamInd))*(tout(2)-tout(1)) );
            end
        end
% %         M_g (:,NamInd)      = GenTq;
% %         M_g_STD (:,NamInd)      = std(GenTq);
        
        %% PSD calculation
        vWindow    = hamming(floor((length(tout(timeinterval)))/8/2)*2);
        dt         = diff(tout(1:2))';
        [S1,f1]      = pwelch(detrend(DataAll(:,Ind),'constant'),vWindow,[],[],1/dt,'onesided'); 
        Fname = [ChannReq{j} '_Fr'];
        Sname = [ChannReq{j} '_Sp'];
        eval([Fname '(:,NamInd)= f1;']);  
        eval([Sname '(:,NamInd)= S1;']);
        if FastExe == 0
            [PitchRate_Sp_Out,PitchRate_Fr_Out] = pwelch(detrend(PitchRate(:,NamInd),'constant'),vWindow,[],[],1/dt,'onesided');
            PitchRate_Fr(:,NamInd) = PitchRate_Fr_Out;
            PitchRate_Sp(:,NamInd) = PitchRate_Sp_Out;
        end 
        %% DEL calculation
%         [6,9]   'WoehlerExponent =  4; N_REF = 2e6; SECINLIFETIME = 20*8760*3600; dt=diff(Time(1:2)); TimeOfBlocks=Time(end)-Time(1)'
%         [6,9]   'Statistics.DEL (iResultFile,iChannel) = ComputeDamageEquivalentLoad(Data,dt,N_REF*TimeOfBlocks/SECINLIFETIME,WoehlerExponent)'                    

        [ext, exttime] = sig2ext(double(1e3*DataAll(:,Ind)),tout(2)-tout(1));
        rf = rainflow(ext,exttime);
        cyclesAmplitude = rf(1,:);
        cyclesNumber    = rf(3,:);
        N_REF     = Nref_basis*(tout(timeinterval(end))-tout(timeinterval(1)))/(20*(365*24)*(60*60));
        DEL4ref   = SF*(sum(cyclesAmplitude.^4.*cyclesNumber)/N_REF).^(1/4);
        DEL10ref  = SF*(sum(cyclesAmplitude.^10.*cyclesNumber)/N_REF).^(1/10);        
        
        DEL4nam   = [ChannReq{j} '_DEL4'];
        DEL10nam  = [ChannReq{j} '_DEL10'];  
        
        eval([DEL4nam '(:,NamInd)= DEL4ref;']);
        eval([DEL10nam '(:,NamInd)= DEL10ref;']);%      
        clear rf ext exttime DEL4ref DEL10ref
    end
end

%% Plot
for i = 1:length(ChannReq)
    
%     Plot TS
    if find(strcmp(TSChan,ChannReq{i})==1)
        figure
        plot (tout(timeinterval),eval([ChannReq{i} '_TS']),'Linewidth',2)

        title (ChannReq{i})
        %         legend(VarName)
        xlabel('Time [s]','FontSize',14)
        set(gca,'FontSize',14);
        grid on
    end
    
%     Plot PSD
%       if find(strcmp(PSDchan,ChannReq{i})==1)
%         figure
%         for j =  1:length(Resfiles)
%             semilogy (eval([ChannReq{i} '_Fr(:,j)']),eval([ChannReq{i} '_Sp(:,j)']),'Linewidth',2);
%             hold on
%         end
%         hold off
%         xlim([0 2])
%         title (ChannReq{i})
%         legend(VarName)
%         grid on 
%       end
    


%    Plot DEL10
%     if find(strcmp(DEL10Chan,ChannReq{i})==1)
%         figure
%         bar (eval([ChannReq{i} '_DEL10']))
%         title ([ChannReq{i} ' DEL Wexp10'])
%         set(gca,'XTickLabel',VarName);
%     end
%     

    %Plot relative DEL wohler 10
    if find(strcmp(DEL10Chan,ChannReq{i})==1)    
         figure

%        if CntPos==1
%           figure1 = figure;
%            axes1 = gca;
%             set(axes1,'OuterPosition',[pos{1}])
%        else
%            axes1 = axes('Parent',figure1,'OuterPosition',[0 0 1 0.5]);
%        end
%         CntPos=CntPos+1;       
        bar (eval(['100*' ChannReq{i} '_DEL10/' ChannReq{i} '_DEL10(1)-100']))
        title ([ChannReq{i} ' Relative DEL Wexp10'])
        set(gca,'XTickLabel',VarName);
        axes1 = gca;
    end
    
      
      
    
    %Plot DEL4
%     if find(strcmp(DEL4Chan,ChannReq{i})==1)       
%         figure
%         bar (eval([ChannReq{i} '_DEL4']))
%         title ([ChannReq{i} ' DEL Wexp4'])
%         set(gca,'XTickLabel',VarName);
%     end  
    
%     %Plot relative DEL wohler 4
%     if find(strcmp(DEL4Chan,ChannReq{i})==1)       
%         figure
% %         if CntPos==1
% %             set(axes1,'OuterPosition',[pos{1}])
% %         else
% %            axes1 = axes('Parent',figure1,'OuterPosition',[0 0 1 0.5]);         
% %         end
% %         CntPos=CntPos+1;        
%         bar (eval(['100*' ChannReq{i} '_DEL4/' ChannReq{i} '_DEL4(1)-100']))
%         title ([ChannReq{i} 'Relative DEL Wexp4'])
%         set(gca,'XTickLabel',VarName); 
%     end
    
    %Plot mean    
%     if find(strcmp(MeanChan,ChannReq{i})==1)      
%         figure
%         bar (eval(['mean(' ChannReq{i} '_TS)']))
%         title ([ChannReq{i} ' Mean'])
%         set(gca,'XTickLabel',VarName);
%     end
    
    %Plot relative mean
%     if find(strcmp(MeanChan,ChannReq{i})==1)  
%          figure
% %         if CntPos==1
% %             set(axes1,'OuterPosition',[pos{1}])
% %         else
% %            axes1 = axes('Parent',figure1,'OuterPosition',[0 0 1 0.5]);          
% %         end
%         CntPos=CntPos+1; 
%         bar (eval(['(100*mean(' ChannReq{i} '_TS)/mean(' ChannReq{i} '_TS(:,1)))-100']))
%         title ([ChannReq{i} 'Relative Mean'])
%         set(gca,'XTickLabel',VarName); 
%     end
    
end

% Plot Pitch rate PSD
%         figure
%         for j =  1:length(Resfiles)
%             plot (PitchRate_Fr(:,j),PitchRate_Sp(:,j),'Linewidth',2);
%             hold on
%         end
%         hold off
%         xlim([0 2])
%         legend(VarName)
%         grid on 
%         title ( 'Pitch Rate PSD')


%Plot relative PitchRate std 
         figure
         bar (100*PitchRete_STD/PitchRete_STD(1)-100)
         title ( 'Relative STD of $\dot{\theta}$', 'Interpreter','latex')
%Plot relative PitchRate max
%          figure
%          bar (100*(PitchAcc_std)/(PitchAcc_std(1))-100)
%          title ( 'Relative mean of $\dot{\theta}$', 'Interpreter','latex')               
%Plot relative Omegastd std          
%          figure
%          bar (100*std(RotSpeed_TS)/std(RotSpeed_TS(:,1))-100)
%          title ( 'Relative STD of rotor $\omega$', 'Interpreter','latex')
%Plot relative Omegastd std          
%          figure
%          bar (100*std(GenSpeed_TS)/std(GenSpeed_TS(:,1))-100)
%          title ( 'Relative STD of generator $\omega$', 'Interpreter','latex')         
         
%Pitch total distance travel deg         
         Pitchtravel    = sum(abs(diff(BldPitch1_TS)));
         PitchtravelRel = 100*Pitchtravel/Pitchtravel(1)-100;
         figure
         bar (PitchtravelRel)
         title ( 'Relative total Pitch travel')

% Total penergy ;
         EnergyProd    = trapz(tout(timeinterval),GenPwr_TS);
         EnergyProdRel = 100*EnergyProd/EnergyProd(1)-100;
         figure
         bar (EnergyProdRel)
         title ( 'Relative total Energy produced')

% Kp=[1e-4:0.25e-4:2.5e-4];Kd= [8e-05:2e-05:24e-05];z=[RootMyb1_DEL10(2:10)', RootMyb1_DEL10(11:19)',RootMyb1_DEL10(20:28)',RootMyb1_DEL10(29:37)',RootMyb1_DEL10(38:46)',RootMyb1_DEL10(47:55)']
% [X,Y]=meshgrid(Kp,Kd);
% z=[RootMyb1_DEL10(2:10)', RootMyb1_DEL10(11:19)',RootMyb1_DEL10(20:28)',RootMyb1_DEL10(29:37)',RootMyb1_DEL10(38:46)',RootMyb1_DEL10(47:55)'];
% figure,contour(X,Y,z)