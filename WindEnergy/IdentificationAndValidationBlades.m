clc 
clear all %#ok<*CLALL>
close all

% Script to identifying linear model from blades and validating it. The input
% files should contain pitch, moment, tout and Vbl_Eff. There can be
% provided different dataset for identification and validation
%
% Dependencies: ValidateIdentifiedSystem function
%
% Vasilis Pettas 12.2017 SWE

%% Parameters

URef_Identify   = 16; % case results from FAST to be identified 
URef_Validate   = 16; % inputs (M,v,theta) from FAST simulation to be fed in the simulation/validation of the identified system

IdentifySystem  = 0; % flag to identify linear system
IdentModelOrder      = 5;
SaveIdentifiedSystem = 1;

SimulateSystem  = 1; % flag to simulate linear system for validation

%% Identification options


InputForID           = {
%                 ['D:\Tasks\Torque2018\Identification\InputsForID\Chirp_V',num2str(URef_Identify),'.mat']
%                 ['D:\Tasks\Torque2018\Identification\InputsForID\Chirp_V',num2str(URef_Identify),'_Delta.mat']
%                 ['D:\Tasks\Torque2018\Identification\InputsForID\Chirp_V',num2str(URef_Identify),'_HPmom.mat']
%                 
%                 ['D:\Tasks\Torque2018\Identification\InputsForID\IecTI_V',num2str(URef_Identify),'.mat']
%                 ['D:\Tasks\Torque2018\Identification\InputsForID\IecTI_V',num2str(URef_Identify),'_Delta.mat']
%                 ['D:\Tasks\Torque2018\Identification\InputsForID\IecTI_V',num2str(URef_Identify),'_HPmom.mat']

%                 ['D:\Tasks\Torque2018\Identification\InputsForID\LowTI_V',num2str(URef_Identify),'.mat']               
                ['D:\Tasks\Torque2018\Identification\InputsForID\LowTI_V',num2str(URef_Identify),'_Delta.mat']
%                 ['D:\Tasks\Torque2018\Identification\InputsForID\LowTI_V',num2str(URef_Identify),'_HPmom.mat']

%                 ['D:\Tasks\Torque2018\Identification\InputsForID\LowTI_V',num2str(URef_Identify),'_Delta_FD.mat']  % mohammad TF
%                 ['D:\Tasks\Torque2018\Identification\models\LowTI_V',num2str(URef_Identify),'_Delta_TF.mat']
                       }; 
%Create linear model name                   
ind1 =  strfind(InputForID,'\');
ind=max(ind1{1});
Model_Name           = ['BladeSysOrd' num2str(IdentModelOrder) '_' InputForID{1}(ind+1:end-4)]; % resulting name of identified linear model

%% Validation Options
InputsForSystemValidation = {
%                             ['D:\Tasks\Torque2018\Identification\ValidationInputsID\IecTI_V',num2str(URef_Validate),'.mat']
%                             ['D:\Tasks\Torque2018\Identification\ValidationInputsID\IecTI_V',num2str(URef_Validate),'_Delta.mat']
%                             ['D:\Tasks\Torque2018\Identification\ValidationInputsID\IecTI_V',num2str(URef_Validate),'_HPmom.mat']

%                             ['D:\Tasks\Torque2018\Identification\ValidationInputsID\LowTI3_V',num2str(URef_Validate),'.mat']
                            ['D:\Tasks\Torque2018\Identification\ValidationInputsID\LowTI3_V',num2str(URef_Validate),'_Delta.mat']
%                             ['D:\Tasks\Torque2018\Identification\ValidationInputsID\LowTI3_V',num2str(URef_Validate),'_HPmom.mat']
                            };
% SystemTobeSimulated         = ['D:\Tasks\Torque2018\Identification\models\' Model_Name]; %#ok<*NBRAK>

%% Identification 
 
if IdentifySystem == 1 
    load (InputForID{1})
    v0      = VBl_Eff(1000:length(tout));
    theta   = pitch(1000:length(tout));
    M       = moment(1000:length(tout)); 

    % Input data for identification
    Blade            = iddata(M,[theta v0],tout(2)-tout(1));     %identification data to be fed in the system identification algorithm 
    Blade.InputName  = {'theta','v0'};
    Blade.OutputName = {'moment'};
    
    % ns4id
    col1 = [6 6 6 9 9 9 15 15 15]';
    col2 = [35 37 39 35 37 39 35 37 39]';%35*ones(length(col1),1);%[5:5:125]';
    col3 = [35 37 39 35 37 39 35 37 39]';% 35*ones(length(col1),1);%[35:10:125]';
    matr = [col1 col2 col3] ;
    
%     Blade_d             = detrend(Blade,0);            % detrend the data                                 
%     Blade_t             = Blade(1000:length(tout));%:floor(length(tout)/2));  % training data
%     Blade_v             = Blade(floor(length(tout)/2)+1:length(tout));        % validation data
    Opt                 = n4sidOptions('N4Weight','CVA', 'N4Horizon',matr);%[15 39 39]);%matr);   %[10:10:100 35:10:135' 35:10:135' ]);%15 35 35;0 0 0;
    Opt.Focus    = 'simulation';
    Opt.N4Weight = 'CVA';
    Opt.Display  = 'on';   
    eval(['[' Model_Name ',x0] = n4sid(Blade,IdentModelOrder,''Ts'',0 ,Opt)']);
    
    if SaveIdentifiedSystem ==1
        save(['D:\Tasks\Torque2018\Identification\models\' Model_Name],Model_Name,'x0')
    end
    
    figure
    compare(Blade,eval(Model_Name),inf)
    isstable(eval(Model_Name))
%     zero(eval(Model_Name))


    figure,iopzmap(eval(Model_Name))
end

%% Simulating system

LinearSystem = ['D:\Tasks\Torque2018\Identification\models\' Model_Name];

if SimulateSystem==1
  ValidateIdentifiedSystem(InputsForSystemValidation,LinearSystem)%    TestIdentifiedsystem(InputsForSystemSimulation,['..\Models\' Model_Name],WindforSimulation,Uref_SysdtemForSimulation,['..\Wind\WIND_DLC12_',num2str(URef_Identify),'_BlEff_',EffBlade_Id])
end

