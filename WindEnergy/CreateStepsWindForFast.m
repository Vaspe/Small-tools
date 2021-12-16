clc; clear all; close all

% Create step wind files


filedirectory = 'N:\SWE\71_Simulationsserver\Rok\34_pettas\PhD\Wind\Steps\';
name      = 'WindStepsBelowRated_SS_DTU10MW';
timeintervals = 120;                      % time intervals for change step
const         = 4;
steps         = 4:0.25:11;%[const:1:9, 9.25:0.25:12, 13:1:24];%[const:1:15]';%[const:0.1:const+0.3,const+0.2:-0.1:const-0.2, const-0.1:0.1:const+0.1]' ;%12:0.1:12.3;%[12:0.1:24,23.9:0.1:12]';   % wind values including initial one
dt            = 0.25 ;                  % time step for wind field

%% 
filename    = [filedirectory name num2str(const) '.wnd'] ;
time        = [0:dt:length(steps)*timeintervals]';     
zerocolumns = zeros(length(time),1);

% locate speed changing time stamps
multi     = 1:length(steps)-1;
StepTimes = multi*timeintervals;

initTInd = 1;
Wspeed   = zeros(length(time),1);
for i = 1:length(StepTimes)
    TimeIndEnd = find (time==StepTimes(i))-1;
    Wspeed (initTInd:TimeIndEnd)  = steps(i); 
    initTInd =TimeIndEnd + 1;   
end

Wspeed (initTInd:end)  = steps(end);

plot(time,Wspeed)
finalvar = [time,Wspeed,zerocolumns,zerocolumns,zerocolumns,zerocolumns,zerocolumns,zerocolumns];
dlmwrite(filename,finalvar,'precision','%.4f','delimiter','\t');

