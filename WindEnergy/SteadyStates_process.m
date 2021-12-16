clearvars
clc
close all

%% INPUTS
folder ='D:\Dropbox\SWE\PhD_Topic\ControllerAndSimulation\SteadyStates\Results\';
cases = {
    'SS_100_04_SS';
    %     'SS_100_04_IPC_SS';
%     'SS_MinPitchDER04_IPC_SS';
%     'SS_MinPitchm1DER04_SS';
%     'SS_TSR78_CP42_DER_04_SS'
    %     'SS_BR_100_04_SS';
    %     'SS_BR_TSR78_CP42_DER_04_SS';
    %     'SS_BR_TSR85_CP4304_OmegaMIn58_SS';
    'SS_all_TSR7804_SS';
    };
for i=1:length(cases)
    filenam{i,1} = [folder cases{i} '.mat'];
end


%% Calculations
for iCase = 1:length(filenam)
    aa{iCase}= load(filenam{iCase});
    WSPInd(iCase) = GetChannelIndex(aa{iCase}.SteadyState.Names, 'WSP');
    GenTqInd(iCase) = GetChannelIndex(aa{iCase}.SteadyState.Names, 'GenTq');
    GenSpInd(iCase) = GetChannelIndex(aa{iCase}.SteadyState.Names, 'GenSpeed'); %#ok<*SAGROW>
end

%% Plotting


% plot all
for i=1:length(aa{1}.SteadyState.Names)
    figure('Name',[aa{1}.SteadyState.Names{i},'_SS'],'NumberTitle','off')
    for  iCase = 1:length(filenam)
        plot(aa{1,iCase}.SteadyState.Values(:,WSPInd(iCase)),aa{1,iCase}.SteadyState.Values(:,i),'linewidth',2),grid on,hold on
    end
    ylabel(aa{1}.SteadyState.Names{i})
    xlabel(aa{1}.SteadyState.Names{WSPInd(1)})
    legend(cases,'Interpreter', 'none')
    hold off
end

% gentq vs gensp
figure('Name','GenChar_SS','NumberTitle','off')
for  iCase = 1:length(filenam)
    plot(aa{1,iCase}.SteadyState.Values(:,GenSpInd(iCase)),aa{1,iCase}.SteadyState.Values(:,GenTqInd(iCase)),'linewidth',2),grid on,hold on
end
legend(cases,'Interpreter', 'none')
ylabel('Gen Torque')
xlabel('Gen Speed')


figure('Name','CT-TSR','NumberTitle','off')
for  iCase = 1:length(filenam)
    plot(aa{1,iCase}.SteadyState.Values(:,11),aa{1,iCase}.SteadyState.Values(:,9),'linewidth',2),grid on,hold on
end
ylabel('CT')
xlabel('TSR')
legend(cases,'Interpreter', 'none')
hold off

figure('Name','CP-TSR','NumberTitle','off')
for  iCase = 1:length(filenam)
    plot(aa{1,iCase}.SteadyState.Values(:,11),aa{1,iCase}.SteadyState.Values(:,8),'linewidth',2),grid on,hold on
end
ylabel('CP')
xlabel('TSR')
legend(cases,'Interpreter', 'none')
hold off

figure('Name','CP-TSR-thera','NumberTitle','off')
for  iCase = 1:length(filenam)
    plot3(aa{1,iCase}.SteadyState.Values(:,11),aa{1,iCase}.SteadyState.Values(:,3),aa{1,iCase}.SteadyState.Values(:,8),'linewidth',2),grid on,hold on
end
ylabel('theta')
xlabel('TSR')
zlabel('CP')
legend(cases,'Interpreter', 'none')
hold off

figure('Name','CT-TSR-theta','NumberTitle','off')
for  iCase = 1:length(filenam)
    plot3(aa{1,iCase}.SteadyState.Values(:,11),aa{1,iCase}.SteadyState.Values(:,3),aa{1,iCase}.SteadyState.Values(:,9),'linewidth',2),grid on,hold on
end
ylabel('theta')
xlabel('TSR')
zlabel('CT')
legend(cases,'Interpreter', 'none')
hold off


%% Functions
function iInd = GetChannelIndex(OutList, Chname)
indfind = strfind (OutList, Chname);
for j = 1:length (indfind)
    if indfind{j}==1
        indmat(j) = 1; %#ok<*AGROW>
    else
        indmat(j) = 0;
    end
end
iInd = find(indmat==1);
end
