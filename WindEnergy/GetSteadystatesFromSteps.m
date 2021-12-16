% Script to get steady state results from FAST simulations with Steps
clc
clear all %#ok<CLALL>
close all

%% Input
RedFolder = 'D:\Dropbox\SWE\PhD_Topic\ControllerAndSimulation\SteadyStates\Results';
ResName = {%'SS_100_04';
                %     'SS_100_04_IPC';
                %     'SS_MinPitchDER04_IPC';
                %     'SS_MinPitchm1DER04';
                %     'SS_TSR78_CP42_DER_04';
%                 'SS_BR_100_04';
%                 'SS_BR_TSR78_CP42_DER_04';
%                 'SS_BR_TSR85_CP4304_OmegaMIn58';
                'SS_all_TSR7804';
                    };
ChanReq = {
    'RotSpeed'
    'GenSpeed'
    'BldPitch1'
    'GenPwr'
    'GenTq'
    'RotThrust'};
StepSize = 120; %Stepsize in seconds
BufferT = 2; %Seconds before the step to probe the time series for SS
rho = 1.225; % density of air kg/m3
R = 89.2;  % rotor radius m
Ploss = 220; %electrical losses kW / 0 if you dont know
eta_el = 0.94; % efficiency of generator / 1 if you dont know


%% Calculations

% Create all file names
for i = 1:length(ResName)
    ResFiles{i} = [RedFolder '\' ResName{i} '_results.mat'];
end

%Extract time series for each channel and
for i=1:length(ResFiles)
    
    % Assign the names of all requested variables (except the coefficents that will be added later)
    for iq =1 :length(ChanReq)
        SteadyState.Names{iq} = ChanReq{iq};
    end
    
    iNamFile = ResFiles{i};
    load(iNamFile) % logsout, parameter and tout are loaded
    
    TSobj    = logsout.getElement('OutData').Values;
    DataAll  = TSobj.OutData.Data;
    for ii=1:length(Parameter.FASTInput.OutList)
        NamList{ii,1}=Parameter.FASTInput.OutList{ii,1};
    end
    
    % get speed vector
    iInd = GetChannelIndex(Parameter.FASTInput.OutList, 'Wind1VelX');
    SpeedVec = DataAll(:,iInd); % this is the actual time series
    
    TimeVec =DataAll(:,1); % time vector, always first in the results file
    StepTimes = StepSize:StepSize:TimeVec(end); %find all the swtching points of steps
    dt = TimeVec(2)-TimeVec(1);
    IndToRemove = ceil(BufferT/dt); % how many time steps to be removed
    [~,TimeStepsInd_init] = intersect(TimeVec,StepTimes);
    TimeStepsInd = TimeStepsInd_init-IndToRemove; % remove time steps
    SpeedSteps = (SpeedVec(TimeStepsInd));
    
    %get the current channel SS values
    for iChan=1:length(ChanReq)
        iInd = GetChannelIndex(Parameter.FASTInput.OutList, ChanReq{iChan});
        curChan = DataAll(:,iInd); % this is the actual time series
        SteadyState.Values(:,iChan) = curChan(TimeStepsInd);
    end
    % Add wind speed
    SteadyState.Names{end+1} = 'WSP';
    SteadyState.Values(:,end+1) = SpeedSteps;
    
    % Add the calculations of the coefficients
    % Calculate Cp
    iInd = GetChannelIndex(Parameter.FASTInput.OutList, 'GenPwr');
    Pwr = DataAll(TimeStepsInd,iInd)*1000; % power in W
    Cp = 2*((Pwr+Ploss)/eta_el)./(rho*pi*R^2*SpeedSteps.^3);
    
    % Calculate CT
    iInd = GetChannelIndex(Parameter.FASTInput.OutList, 'RotThrust');
    Thrust = DataAll(TimeStepsInd,iInd)*1000;
    Ct =Thrust./(0.5*rho*pi*R^2*(SpeedSteps.^2));
    
    % Calculate lambda
    iInd = GetChannelIndex(Parameter.FASTInput.OutList, 'RotSpeed');
    RotSpeed = DataAll(TimeStepsInd,iInd)*(2*pi/60); % rot speed in rad/s
    lambda = RotSpeed*R./SpeedSteps;
    
    % Calculate Cm
    Cm = Cp./(lambda);
    
    % Add them to output
    SteadyState.Names{end+1} = 'Cp';
    SteadyState.Values(:,end+1) = Cp;
    SteadyState.Names{end+1} = 'Ct';
    SteadyState.Values(:,end+1) = Ct;
    SteadyState.Names{end+1} = 'Cm';
    SteadyState.Values(:,end+1) = Cm;
    SteadyState.Names{end+1} = 'lambda';
    SteadyState.Values(:,end+1) = lambda;
    SteadyState.Names{end+1} = 'Thrust';
    SteadyState.Values(:,end+1) = Thrust;
    
    % Save steady states output
    save([RedFolder '\' ResName{i} '_SS'],'SteadyState')
    clear SteadyState
end


function iInd = GetChannelIndex(OutList, Chname)
indfind = strfind (OutList, Chname);
for j = 1:length (indfind)
    if indfind{j}==1
        indmat(j) = 1;
    else
        indmat(j) = 0;
    end
end
iInd = find(indmat==1);
end

