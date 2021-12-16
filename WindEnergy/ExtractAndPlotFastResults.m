clc 
clear all %#ok<*CLALL>
close all

% ResultForIdentify ={
% 'D:\Tasks\Torque2018\DynResponse\FAST_Blades_dyn_resp_steps_V_16_CL1\FAST_Blades_dyn_resp_steps_V_16_CL1_results.mat'
% 'D:\Tasks\Torque2018\DynResponse\FAST_Blades_dyn_resp_steps_V_16_CL1_flatV\FAST_Blades_dyn_resp_steps_V_16_CL1_flatV_results.mat'
% 'D:\Tasks\Torque2018\DynResponse\FAST_Blades_dyn_resp_steps_V_16_OL1_flatV_NoPitch\FAST_Blades_dyn_resp_steps_V_16_OL1_flatV_NoPitch_results.mat'
% 'D:\Tasks\Torque2018\DynResponse\FAST_Blades_dyn_resp_steady_V_16_OL1_flatV_Pitch_steps01\FAST_Blades_dyn_resp_steady_V_16_OL1_flatV_Pitch_steps01_results.mat'
% 'D:\Tasks\Torque2018\DynResponse\FAST_Blades_dyn_resp_steady_V_16_OL1_Pitch_harmonic_f01\FAST_Blades_dyn_resp_steady_V_16_OL1_Pitch_harmonic_f01_results.mat'
% 'D:\Tasks\Torque2018\DynResponse\FAST_Blades_dyn_resp_steps_V_16_OL1_DLC12_NoPitch\FAST_Blades_dyn_resp_steps_V_16_OL1_DLC12_NoPitch_results.mat'
% 'N:\SWE\90_Transfer\Mohammad Salari\Vasilis Torque\FBSWE_FAST_DynIn_Stea_NT_DLC12_16_results.mat'
% 'D:\Tasks\Torque2018\DynResponse\FAST_Blades_dyn_resp_Kaimal_V_16_CL1\FAST_Blades_dyn_resp_Kaimal_V_16_CL1_results.mat'
% 'D:\Tasks\Torque2018\DynResponse\FAST_Blades_dyn_resp_Kaimal_V_CL1\FAST_Blades_dyn_resp_Kaimal_CL1_V17_results.mat'
% 'D:\Tasks\Torque2018\DynResponse\FAST_Blades_dyn_resp_steps_CL1\FAST_Blades_dyn_resp_steps_CL1_16_results.mat'
% 'D:\Tasks\Torque2018\DynResponse\FAST_Blades_dyn_resp_Kaimal_V_CL1\FAST_Blades_dyn_resp_Kaimal_CL1_V16_results.mat'
% 'D:\Tasks\Torque2018\DynResponse\FAST_Blades_dyn_resp_steps_V_16_CL1_No_tower\FAST_Blades_dyn_resp_steps_V_16_CL1_No_Tower_results.mat'
% 'D:\Tasks\Torque2018\DynResponse\FAST_Blades_dyn_resp_steps_V_16_CL1_No_tower_Nobl\FAST_Blades_dyn_resp_steps_V_16_CL1_No_Tower_Nobl_results.mat'
% 'D:\Tasks\Torque2018\DynResponse\FAST_Blades_dyn_resp_steps_CL1_allOff\FAST_Blades_dyn_resp_steps_CL1_allOff17_results.mat'
% 'D:\Tasks\Torque2018\DynResponse\FAST_Blades_dyn_resp_stepsAll_CL1_allOff\FAST_Blades_dyn_resp_stepsAll_CL1_allOff12_results.mat'
% 'D:\Tasks\Torque2018\DynResponse\FAST_Blades_dyn_resp_stepsAll_CL1_allOn\FAST_Blades_dyn_resp_stepsAll_CL1_allOff12_results.mat'
% 'D:\Tasks\BasicSWEController\DTU10MW_FBSWE2_DLC1d2\WIND_DLC12_6_72_results.mat'
% 'D:\Tasks\Torque2018\DynResponse\FAST_BladeSpeed_ID_DLC12\FAST_BladeSpeed_ID_DLC12_16_results.mat'
%'D:\DynamicInflowSims\DLC12simulations\DTU10MW_FB_DLC12_FAST_Equil_Stea_NT\FBSWE_FAST_Equil_Stea_NT_DLC12_12_results.mat'
% };

% Choose channels
ChannReq = {
%      'GenPwr';
%     'TipDxb1';
%     'TipALxb1';    
%     'RootMzb1';
%     'RootMyb1';
%      'RootMxb1';
%     'TwrBsMxt'
%     'TwrBsMyt'
%     'BldPitch1';
    'Wind1VelX';
%     'RotSpeed';
% %     'LSShftTq';
%     'TwrBsMyt'
%      'TTDspFA';
    'Azimuth'
    };

% Parameters
URef_Identify             = 16; % case results from FAST to be identified 
URef_ComparisonInput      = 16; % inputs (M,v,theta) from FAST simulation to be fed in the simulation of the identified system
Uref_SysdtemForSimulation = 16; % system name to be simulated

% Flags
CalculateAzimuth    = 1;
PlotTimeseries      = 1;
PlotPSD             = 0;
GetRotEffWindSspeed = 0;
GetBlade1TipSpeed   = 0;
GetBladeEffSpeed    = 1;
EffBlade_Id         = 'Linear_Mean';%'Linear_Quadratic'; %'NoWeight_Qubic'; % 'CpWeight_Mean';
Getblad3_V_eff      = 0;

%Identification options
IdentifySystem       = 0;
VStr                 = 'VBlEff';%'Vtip';%'Vhub'; 'VRotEff'; 
IdentModelOrder      = 5;
SaveIdentifiedSystem = 1;
Model_Name           = ['BladeSysOrd' num2str(IdentModelOrder) 'V' num2str(URef_Identify) '_VBlEff_LowTI'];

%Simulation Options
SimulateSystem       = 0;
InputsForSystemSimulation = ['D:\Tasks\Torque2018\DynResponse\FAST_BladeSpeed_ID_DLC12\FAST_BladeSpeed_ID_DLC12_' num2str(URef_ComparisonInput) '_results.mat'];
% InputsForSystemSimulation = ['D:\Tasks\Torque2018\DynResponse\FAST_BladeSpeed_IDbase_LowTI\FAST_BladeSpeed_IDbase_LowTI_' num2str(URef_ComparisonInput) '_results.mat'];

% InputsForSystemSimulation = ['..\results\FAST_Blades_dyn_resp_Kaimal_CL1_V' num2str(URef_ComparisonInput) '_results.mat'];
% SystemTobeSimulated       = ['..\Models\BladeSysOrd' num2str(IdentModelOrder) 'V' num2str(Uref_SysdtemForSimulation) '_Vtip'] ; 
SystemTobeSimulated         = [Model_Name EffBlade_Id ];
WindforSimulation         = 'BldEffV';%'RotEff'%'BldTipV'; ; %'Wind1VelX' %  

 ResultForIdentify = {['D:\Tasks\Torque2018\DynResponse\FAST_BladeSpeed_IDbase_LowTI\FAST_BladeSpeed_IDbase_LowTI3_',num2str(URef_Identify),'_results.mat']
%      ['D:\Tasks\Torque2018\DynResponse\FAST_BladeSpeed_ID_DLC12\FAST_BladeSpeed_ID_DLC12_',num2str(URef_Identify),'_results.mat']
%                      {['..\results\FAST_Blades_dyn_resp_Kaimal_CL1_V',num2str(URef_Identify),'_results.mat']};
                       }; 

%% Get timeseries and PSD from FAST results

load(ResultForIdentify{1})
DataAll          = logsout.getElement('OutData').Values.Data;
Chanels          = Parameter.FASTInput.OutList  ;

for j =1:length (ChannReq)
    indfind = strfind (Chanels, ChannReq{j});

    for i = 1:length (indfind)
        if indfind{i}==1
            indmat(i) = 1; %#ok<*SAGROW>
        else
            indmat(i) = 0;
        end
    end

    Ind = find(indmat==1);
    eval([ChannReq{j} '= DataAll(:,Ind);']); %Create variable with time series
    
    % PSD calculation
    vWindow    = hamming(floor((length(tout)-1000)/8/2)*2);
    dt         = diff(tout(1:2))';
    [S,f]      = pwelch(detrend(DataAll(:,Ind),'constant'),vWindow,[],[],1/dt,'onesided'); 
   
    %plot time series
    if PlotTimeseries == 1;
        figure 
        plot(tout,DataAll(:,Ind));
        legend (ChannReq{j})
        grid on
    end
   
    % plot PSD
    if PlotPSD == 1;
        figure 
        semilogy(f,S);
        xlim([0 2])
        grid on
        legend (ChannReq{j})
    end    
end

%% Azimuth calculation

if CalculateAzimuth == 1;
    indfind = strfind (Chanels, 'RotSpeed');
    for i = 1:length (indfind)
        if indfind{i}==1
            indmat(i) = 1;
        else
            indmat(i) = 0;
        end
    end
    Ind     = find(indmat==1);
    Omega   = DataAll(:,Ind);
    dtheta = (tout(2)-tout(1)).* Omega*6; %from rpm to deg/s
    % theta_dif = diff(theta_int);
    Azimuth_V     = 0+ cumsum(dtheta); % azimuth in degrees, 0 upper position, initial
    Azimuth_V     = rem (Azimuth_V,360);
    Azimuth_V     = [0 ; Azimuth_V(1:end-1)] ;

%     figure 
%     scatter(Azimuth,RootMyb1);
%     legend ('Azimuth')
%     figure 
%     plot(tout,Azimuth_V);
%     legend ('Azimuth_V')
    % 
%     figure
%     plotyy(tout,Azimuth,tout,RootMyb1);
%     legend ('Azimuth')
%     grid on
%      figure;plot(tout,Azimuth_V,tout,(rem(9.6*6*tout,360)))
end

%% Get rotor effective wind speed

% the rotor effective wind speed is found from the pre processed wind files directly
if GetRotEffWindSspeed == 1
    
    V_load  = load(['..\Wind\WIND_DLC12_',num2str(URef_Identify),'_results']);
    v0_lowSampleRate = V_load.Disturbance.v_0.signals.values;
    t0 = V_load.Disturbance.v_0.time;

    t_positiveTime = t0(t0>=0); 
    v0_positiveTime = v0_lowSampleRate(t0>=0);

%     t_new = interp1(t_positiveTime,t_positiveTime,tout);
    V_RotEff = interp1(t_positiveTime,v0_positiveTime,tout);
end

%% Get Blade 1 Tip speed
if GetBlade1TipSpeed == 1
    VTip_load  = load(['..\Wind\WIND_DLC12_',num2str(URef_Identify),'_TipBl1']);
    Vtip       = eval(['VTip_load.TipV_' num2str(URef_Identify)]);
    if isnan(Vtip(end))
        Vtip(end)= Vtip(end-1);
    end
%     figure; 
%     plot(tout,Wind1VelX,tout,Vtip,tout,V_RotEff);legend({'Hub Height Wind' 'Blade 1 Tip Wind' 'Rotor Effective Wind Speed'})    
end

%% Get BladeEffective Wind speed
if GetBladeEffSpeed   == 1;
    
    %Get Tower Velocity
    indfind = strfind (Chanels, 'TTDspFA');
    for i = 1:length (indfind)
        if indfind{i}==1
            indmat(i) = 1;
        else
            indmat(i) = 0;
        end
    end
    Ind     = find(indmat==1);
    TTpos   = DataAll(:,Ind);
    TTVy = diff(TTpos)./(tout(2)-tout(1)) ;
    TTVy = [0;TTVy];
    
%     VBlEff_load  = load(['D:\WindFiles\DLC12wind\WIND_DLC12_',num2str(URef_Identify),'_BlEff_',EffBlade_Id]); 
    VBlEff_load  = load(['D:\WindFiles\LowTI_Iref_04\DTU10MW_LowTI3_',num2str(URef_Identify),'_1_BlEff_',EffBlade_Id]);  
    %Locate azimuth and time and interpolate
    VBl_Eff1      = interp2(VBlEff_load.v_BlEff.time,VBlEff_load.v_BlEff.Azimuth,VBlEff_load.v_BlEff.signals.values',tout, Azimuth_V );
    VBl_Eff      = VBl_Eff1-TTVy;  %add tower top velocity
    
    %Calculate 3blades
    if Getblad3_V_eff==1
       VBlEff2    = interp2(VBlEff_load.v_BlEff.time,VBlEff_load.v_BlEff.Azimuth,VBlEff_load.v_BlEff.signals.values',tout, rem((Azimuth_V+120),360) )-TTVy;
       VBlEff3    = interp2(VBlEff_load.v_BlEff.time,VBlEff_load.v_BlEff.Azimuth,VBlEff_load.v_BlEff.signals.values',tout, rem((Azimuth_V+240),360) )-TTVy;
       VBlEff_all.values = [VBl_Eff,VBlEff2,VBlEff3];       
       VBlEff_all.time   = tout;
       save (['D:\DynamicInflowSims\DLC12wind\WIND_DLC12_' num2str(URef_Identify) '_VBlEff_all_Linear_Mean' ],'VBlEff_all')
    end
%     figure; 
%     plot(tout,Wind1VelX,tout,VBl_Eff);legend({'Hub Height Wind' 'Blade 1 Tip Wind' 'Rotor Effective Wind Speed' 'Blade Effective Wind speed'})    
%     ,tout,Vtip,tout,V_RotEff
end

%% Bandpass filter 
w1p               = 0.15*(2*pi);
z1_1p             = 0;
z2_1p             = 1;
NotchFilter_1P    = tf([1 2*z1_1p*w1p w1p^2],[1 2*z2_1p*w1p w1p^2]);
MyBlfil           = lsim(NotchFilter_1P,RootMyb1,tout);
% 
% figure 
% plot(tout,MyBlfil,tout,RootMyb1);
% legend (ChannReq{j})

% figure
% h = bodeplot(NotchFilter_1P); 
% setoptions(h,'FreqUnits','Hz','PhaseVisible','off'); 

%% Low/High pass filter
fc_low = 0.05;
a1     = fc_low*dt/(1+fc_low*dt);
mx_filt_low(1) = 0; 
for i=2:length(tout)
    mx_filt_low(i)=(1-a1)*mx_filt_low(i-1)+a1*(RootMyb1(i)); 
end

% figure 
% plot(tout,RootMyb1,tout,mx_filt_low,'LineWidth',2),grid on
% legend ({'measured' ['LP filtered fc =' num2str(fc_low)]})
fc_high=0.45;
mx_filt_high(1)=0; 
a=1/(1+fc_high*dt);
for i=2:length(tout)
        mx_filt_high(i)=a*mx_filt_high(i-1)+a*(RootMyb1(i)-RootMyb1(i-1)); 
end

%% Identification 
 
if IdentifySystem == 1 
    
    if strcmp(VStr, 'Vtip')
        v0      = Vtip;
    elseif strcmp(VStr,  'Vhub')
    	v0      = WindVelX;    
    elseif strcmp(VStr, 'VRotEff')
        v0      = V_RotEff;    
    elseif strcmp(VStr, 'VBlEff')    
        v0      = VBl_Eff;
    end
    
    if ~strcmp(VStr, 'VBlEff')  
        Model_Name =  [Model_Name VStr ]; 
    else
        Model_Name =  [Model_Name EffBlade_Id ]; 
    end
    theta   = BldPitch1;
    M       = RootMyb1; %mx_filt_high';% mx_filt_low' ;   %MyBlfil; %M_movAv
    M_movAv = moving(M,1500); %moving average from matlab central change nuber of samples (removing steady state fluctuations)
    
    % figure
    % plot(tout,M,tout,M_mean)
    % legend({'measured Mbfl', 'Averaged Mbfl'});     

    % Input data for identification
    Blade            = iddata(M,[theta v0],tout(2)-tout(1));     %identification data to be fed in the system identification algorithm 
    Blade.InputName  = {'theta','v0'};
    Blade.OutputName = {'moment'};
    
    % ns4id
    col1 = [6 6 6 9 9 9 15 15 15]';
    col2 = [35 37 39 35 37 39 35 37 39]';%35*ones(length(col1),1);%[5:5:125]';
    col3 = [35 37 39 35 37 39 35 37 39]';% 35*ones(length(col1),1);%[35:10:125]';
    matr = [col1 col2 col3] ;
    
    Blade_d             = detrend(Blade,0);            % detrend the data                                 
    Blade_t             = Blade(1000:length(tout));%:floor(length(tout)/2));  % training data
    Blade_v             = Blade(floor(length(tout)/2)+1:length(tout));        % validation data
    Opt                 = n4sidOptions('N4Weight','CVA', 'N4Horizon',matr);%[15 39 39]);%matr);   %[10:10:100 35:10:135' 35:10:135' ]);%15 35 35;0 0 0;
    Opt.Focus    = 'simulation';
    Opt.N4Weight = 'CVA';
    Opt.Display  = 'on';   
    eval(['[' Model_Name ',x0] = n4sid(Blade_t,IdentModelOrder,''Ts'',0 ,Opt)']);
    
    if SaveIdentifiedSystem ==1
        %  save('N:\SWE\90_Transfer\Mohammad Salari\Vasilis Torque\Models\Blade_n4s2_V16_Kaimal4th','Blade_n4s2','x0')
        save(['..\Models\Test\' Model_Name],Model_Name,'x0')
    end
    
    figure
    compare(Blade_v,eval(Model_Name),inf)
%     T1 = BladeSysOrd3V16_VBlEff_testNoWeight_Mean(1, 1) ; 
%     isstable(T1)
%     zero(T1)
end

%% Simulating system

if SimulateSystem==1
    
%     TestIdentifiedsystem(InputsForSystemSimulation,['..\Models\Test\' Model_Name],WindforSimulation,Uref_SysdtemForSimulation,['..\Wind\WIND_DLC12_',num2str(URef_Identify),'_BlEff_',EffBlade_Id]);%['..\Wind\DTU10MW_LowTI_',num2str(URef_Identify),'_BlEff_',EffBlade_Id])
  TestIdentifiedsystem(['D:\Tasks\Torque2018\DynResponse\FAST_BladeSpeed_IDbase_LowTI\FAST_BladeSpeed_IDbase_LowTI_',num2str(URef_Identify),'_results.mat'],['..\Models\Test\' Model_Name],WindforSimulation,Uref_SysdtemForSimulation,['..\Wind\DTU10MW_LowTI_',num2str(URef_Identify),'_BlEff_',EffBlade_Id]);%['..\Wind\DTU10MW_LowTI_',num2str(URef_Identify),'_BlEff_',EffBlade_Id])   
%    TestIdentifiedsystem(InputsForSystemSimulation,['..\Models\' Model_Name],WindforSimulation,Uref_SysdtemForSimulation,['..\Wind\WIND_DLC12_',num2str(URef_Identify),'_BlEff_',EffBlade_Id])
%    TestIdentifiedsystem(InputsForSystemSimulation,SystemTobeSimulated,WindforSimulation,Uref_SysdtemForSimulation,['..\Wind\TurbwindHHKaimal_16_BlEff_',num2str(URef_Identify),EffBlade_Id])
end











