clc,clear all, close all  %#ok<*CLALL>

% script for designing filter for the IPC controller. filters applied in
% moment and pitch signal are checked with time series from FAST
% simulation in order to visualize and verify the filter implementation.
% The parameters of filters are guven in the parameter section and
% different plotting flags are used for figure definition. Input time
% series can be defined in order to see the effect of the filters in TS and
% PSD.
%
% Vasilis Pettas 2.2018 SWE


%% get inputs 
Uref   = 16; 

% load time series for testing the filteers
% load  (['D:\Tasks\Torque2018\Identification\ValidationInputsID\IecTI_V',num2str(Uref) '.mat'])
load  (['D:\Tasks\Torque2018\Identification\ValidationInputsID\IecTI_V',num2str(Uref) '_IPC.mat'])

%% Parameters for filter design 

LP_fc   = 2.5;    %Hz 1st order LP pitch commanded input to pitch actuator
HP_fc   = 0.05; %Hz HP 1st order blade root moment input to ctrl
fNO_s   = 0.32; % Hz target frequency of simple notch filter pitch commanded input to pitch actuator
fNO2_s  = 0.48; % Hz target frequency of simple notch filter pitch commanded input to pitch actuator
% notch2 3parameter david
fNO4_s      = 0.32;  % Hz targeted frequency
fNO4_depth  = 0.4;   % depth of notch
NO4_width   = 0.5;   % width of notch
fNO42_s      = 0.48;  % targeted frequency second notch

% notch 3parameter david high frequencies instead of LP
fNO5_s      = 0.9;  % Hz targeted frequency 0.8
fNO5_depth  = 0.3;   % depth of notch  0.4
NO5_width   = 7;   % width of notch 5

% second notch 3parameter david for even higher frequencies instead of LP
fNO6_s      = 1.5;  % Hz targeted frequency 0.8
fNO6_depth  = 0.4;   % depth of notch  0.4
NO6_width   = 5;   % width of notch 5

QNO        = 25;   % Q factor for notch filter, determines the steepness of the filter (see wikipedia)
LP_fc_comb = 2; % Hz 1st order LP combined with notch for pitch commanded input to pitch actuator
D_LP2      = 1.2; % damping of second order LP filter for pitch 

%% flags for plotting
plotLP   = 0; % low pass filter for pitch after controller and before actuator. Used to cut off high frequency actuation
plotHP   = 0; % high pass filter for moments before entering the controller. Used to remove low freqs (steady state loading), preventing interaction with baseline
plotNO1  = 0; % notch filter for removing frequencies of pitch actuation. Targets controller desired bandwidth 
plotNO4  = 0; % notch filter for removing frequencies of pitch actuation. 3Parametwer from david. Targets controller desired bandwidth 
plotBode = 1; % flag to plot the bode response of the filters
plotNO3  = 0; % multiple frequency notch filter for removing frequencies of pitch actuation. Targets controller desired bandwidth 
plotNO42 = 0; % multiple frequency notch filter for removing frequencies of pitch actuation using david's 3 parameter filter
plotNO5  = 0;  % high frequency focused 3DOF notch filter for replacing LP 
plotNO6  = 1;  % second in series high frequency focused 3DOF notch filter for replacing LP 
plotCOMB = 0; % Combined filter with double noth and pitch for pitch actuation bandwidth restriction
plotCOMB2 = 0; % Combined with Davids notch
plotLP2  = 0; % low pass 2nd order filter for pitch after controller and before actuator. Used to cut off high frequency bandwidth actuation

%% IFFt inputs

vWindow    = hamming(floor((length(tout)-1000)/8/2)*2);
dt         = diff(tout(1:2))';
[Sraw_mom,fraw_mom]     = pwelch(detrend(moment,'constant'),vWindow,[],[],1/dt,'onesided'); 
[Sraw_pitch,fraw_pitch] = pwelch(detrend(pitch,'constant'),vWindow,[],[],1/dt,'onesided'); 

%% HP filter first order

% first Order HP filter
s       = tf('s');
tauHP   = 1/(HP_fc*2*pi);
F_HP    = (tauHP*s ) / (1 + tauHP*s );
HP1_mom = lsim(F_HP,moment,tout);
[S_HP1_mom,f_HP1_mom]  = pwelch(detrend(HP1_mom,'constant'),vWindow,[],[],1/dt,'onesided');

%% LP filter first order
% first Order HP filter

tauLP = 1/(LP_fc*2*pi);    % s
F_LP  = 1 / (1 + s*tauLP); % TF

LP1_mom   = lsim(F_LP,moment,tout);
LP1_pitch = lsim(F_LP,pitch,tout);
[S_LP1_mom,f_LP1_mom]      = pwelch(detrend(LP1_mom,'constant'),vWindow,[],[],1/dt,'onesided'); 
[S_LP1_pitch,f_LP1_pitch]  = pwelch(detrend(LP1_pitch,'constant'),vWindow,[],[],1/dt,'onesided'); 

%% LP filter second order

wn_LP2 = (LP_fc*2*pi) ; % rad/s
F_LP2 = wn_LP2^2 / (s^2 + s*2*wn_LP2*D_LP2 + wn_LP2^2);

LP2_mom   = lsim(F_LP2,moment,tout);
LP2_pitch = lsim(F_LP2,pitch,tout);
[S_LP2_mom,f_LP2_mom]      = pwelch(detrend(LP2_mom,'constant'),vWindow,[],[],1/dt,'onesided'); 
[S_LP2_pitch,f_LP2_pitch]  = pwelch(detrend(LP2_pitch,'constant'),vWindow,[],[],1/dt,'onesided'); 

%% Notch filter 1

% simple 2 pole 2 zero only parameter: Wn 

wn   = 2*pi*fNO_s;
F_NO1 = (s^2+wn^2) / (s^2 + s*(1/QNO)*wn + wn^2);

Noch1_mom   = lsim(F_NO1,moment,tout);
Noch1_pitch = lsim(F_NO1,pitch,tout);
[S_NO1_mom,f_NO1_mom]     = pwelch(detrend(Noch1_mom,'constant'),vWindow,[],[],1/dt,'onesided'); 
[S_NO1_pitch,f_NO1_pitch] = pwelch(detrend(Noch1_pitch,'constant'),vWindow,[],[],1/dt,'onesided'); 

%% Notch filter 2(3DOF) David

wn_4_Dav   = 2*pi*fNO4_s;
FNO4_Dav = (s^2 + s*NO4_width*fNO4_depth + wn_4_Dav^2) / (s^2 + s*NO4_width + wn_4_Dav^2);

Noch4_mom   = lsim(FNO4_Dav,moment,tout);
Noch4_pitch = lsim(FNO4_Dav,pitch,tout);
[S_NO4_mom,f_NO4_mom]     = pwelch(detrend(Noch4_mom,'constant'),vWindow,[],[],1/dt,'onesided'); 
[S_NO4_pitch,f_NO4_pitch] = pwelch(detrend(Noch4_pitch,'constant'),vWindow,[],[],1/dt,'onesided');

%% Multi Notch filter 2(3DOF) David 

wn_42_Dav   = 2*pi*fNO42_s;
FNO42_Dav = FNO4_Dav*(s^2 + s*NO4_width*fNO4_depth + wn_42_Dav^2) / (s^2 + s*NO4_width + wn_42_Dav^2);

Noch42_mom   = lsim(FNO42_Dav,moment,tout);
Noch42_pitch = lsim(FNO42_Dav,pitch,tout);
[S_NO42_mom,f_NO42_mom]     = pwelch(detrend(Noch42_mom,'constant'),vWindow,[],[],1/dt,'onesided'); 
[S_NO42_pitch,f_NO42_pitch] = pwelch(detrend(Noch42_pitch,'constant'),vWindow,[],[],1/dt,'onesided');

%% Notch filter (3DOF) for high frequencies (replacing LP)

wn_5_Dav   = 2*pi*fNO5_s;
FNO5_Dav = (s^2 + s*NO5_width*fNO5_depth + wn_5_Dav^2) / (s^2 + s*NO5_width + wn_5_Dav^2);

Noch5_mom   = lsim(FNO5_Dav,moment,tout);
Noch5_pitch = lsim(FNO5_Dav,pitch,tout);
[S_NO5_mom,f_NO5_mom]     = pwelch(detrend(Noch5_mom,'constant'),vWindow,[],[],1/dt,'onesided'); 
[S_NO5_pitch,f_NO5_pitch] = pwelch(detrend(Noch5_pitch,'constant'),vWindow,[],[],1/dt,'onesided');

%% Second in series Notch filter (3DOF) for high frequencies (replacing LP)

wn_6_Dav   = 2*pi*fNO6_s;
FNO6_Dav = FNO5_Dav*(s^2 + s*NO6_width*fNO6_depth + wn_6_Dav^2) / (s^2 + s*NO6_width + wn_6_Dav^2);

Noch6_mom   = lsim(FNO6_Dav,moment,tout);
Noch6_pitch = lsim(FNO6_Dav,pitch,tout);
[S_NO6_mom,f_NO6_mom]     = pwelch(detrend(Noch6_mom,'constant'),vWindow,[],[],1/dt,'onesided'); 
[S_NO6_pitch,f_NO6_pitch] = pwelch(detrend(Noch6_pitch,'constant'),vWindow,[],[],1/dt,'onesided');

%% multi Notch 1

wn2   = 2*pi*fNO2_s;

F_NO2 = (s^2+wn2^2) / (s^2 + s*(1/QNO)*wn2 + wn2^2);
F_NO3 = F_NO1*F_NO2;
Noch3_mom   = lsim(F_NO3,moment,tout);
Noch3_pitch = lsim(F_NO3,pitch,tout);

[S_NO3_mom,f_NO3_mom]     = pwelch(detrend(Noch3_mom,'constant'),vWindow,[],[],1/dt,'onesided'); 
[S_NO3_pitch,f_NO3_pitch] = pwelch(detrend(Noch3_pitch,'constant'),vWindow,[],[],1/dt,'onesided'); 

%% Combined LP with multiple notch

tauLP_comb = 1/(LP_fc_comb*2*pi);
F_LP_comb  = 1 / (1 + s*tauLP_comb);

F_COMB = F_NO3*F_LP_comb;
COMB_mom   = lsim(F_COMB,moment,tout);
COMB_pitch = lsim(F_COMB,pitch,tout);

[S_COMB_mom,f_COMB_mom]     = pwelch(detrend(COMB_mom,'constant'),vWindow,[],[],1/dt,'onesided'); 
[S_COMB_pitch,f_COMB_pitch] = pwelch(detrend(COMB_pitch,'constant'),vWindow,[],[],1/dt,'onesided');

%% Combined David's filter with LP

F_COMB2 = FNO42_Dav*F_LP_comb;
COMB2_mom   = lsim(F_COMB2,moment,tout);
COMB2_pitch = lsim(F_COMB2,pitch,tout);

[S_COMB2_mom,f_COMB2_mom]     = pwelch(detrend(COMB2_mom,'constant'),vWindow,[],[],1/dt,'onesided'); 
[S_COMB2_pitch,f_COMB2_pitch] = pwelch(detrend(COMB2_pitch,'constant'),vWindow,[],[],1/dt,'onesided');

%% plotting 

opts = bodeoptions;
opts.Title.FontSize = 12;
opts.FreqUnits = 'Hz';
opts.MagUnits = 'abs';
    
if plotLP == 1

    figure,
    bode(F_LP,opts);
    xlim([0.0001 10]),grid on, 
    legend({['LP 1st order filter fc =' num2str(LP_fc) ' Hz']})    

%     figure
%     plot (tout,moment,tout,LP1_mom,'LineWidth',2),grid on
%     set(gca,'FontSize',14);
%     legend({'Raw' ['LP filter fc =' num2str(LP_fc) ' Hz']})
%     ylabel ('Moment')
%     xlabel ('Time')

%     figure
%     semilogy(fraw_mom,Sraw_mom,f_LP1_mom,S_LP1_mom,'LineWidth',2), xlim([0 4]),grid on,  
%     set(gca,'FontSize',14);
%     legend({'Raw' ['LP filter fc =' num2str(LP_fc) ' Hz']})
%     ylabel ('Moment')
%     xlabel ('Frequency')

    figure
    plot (tout,pitch,tout,LP1_pitch,'LineWidth',2),grid on
    set(gca,'FontSize',14);
    legend({'Raw' ['LP filter fc =' num2str(LP_fc) ' Hz']})
    ylabel ('Pitch')
    xlabel ('Time')

    figure
    semilogy(fraw_pitch,Sraw_pitch,f_LP1_pitch,S_LP1_pitch,'LineWidth',2), xlim([0 4]),grid on,  
    set(gca,'FontSize',14);
    legend({'Raw' ['LP filter fc =' num2str(LP_fc) ' Hz']})
    ylabel ('Pitch')
    xlabel ('Frequency')

end

if plotLP2 == 1

    figure,
    bode(F_LP2,opts);
    xlim([0.0001 10]),grid on, 
    legend({['LP 2nd order filter fc =' num2str(LP_fc) ' Hz and damping D =' num2str(D_LP2)]})    

%     figure
%     plot (tout,moment,tout,LP2_mom,'LineWidth',2),grid on
%     set(gca,'FontSize',14);
%     legend({'Raw' ['LP 2nd order filter fc =' num2str(LP_fc) ' Hz and damping D =' num2str(D_LP2)]})
%     ylabel ('Moment')
%     xlabel ('Time')

%     figure
%     semilogy(fraw_mom,Sraw_mom,f_LP1_mom,S_LP1_mom,f_LP2_mom,S_LP2_mom,'LineWidth',2), xlim([0 4]),grid on,  
%     set(gca,'FontSize',14);
%     legend({'Raw' ['LP filter fc =' num2str(LP_fc) ' Hz'] ['LP 2nd order filter fc =' num2str(LP_fc) ' Hz and damping D =' num2str(D_LP2)]})
%     ylabel ('Moment')
%     xlabel ('Frequency')

    figure
    plot (tout,pitch,tout,LP1_pitch,tout,LP2_pitch,'LineWidth',2),grid on
    set(gca,'FontSize',14);
    legend({'Raw' ['LP filter fc =' num2str(LP_fc) ' Hz'] ['LP 2nd order filter fc =' num2str(LP_fc) ' Hz and damping D =' num2str(D_LP2)]})
    ylabel ('Pitch')
    xlabel ('Time')

    figure
    semilogy(fraw_pitch,Sraw_pitch,f_LP1_pitch,S_LP1_pitch,f_LP2_pitch,S_LP2_pitch,'LineWidth',2), xlim([0 4]),grid on,  
    set(gca,'FontSize',14);
    legend({'Raw' ['LP filter fc =' num2str(LP_fc) ' Hz']  ['LP 2nd order filter fc =' num2str(LP_fc) ' Hz and damping D =' num2str(D_LP2)]})
    ylabel ('Pitch')
    xlabel ('Frequency')

end

if plotHP == 1

    figure,
    bode(F_HP,opts)
    xlim([0.0001 10]),grid on, 
    legend({['HP 1st order filter fc =' num2str(HP_fc) ' Hz']})    
    
%     figure
%     plot (tout,moment,tout,HP1_mom,'LineWidth',2),grid on
%     set(gca,'FontSize',14);
%     legend({'Raw' ['HP 1st order filter fc =' num2str(HP_fc) ' Hz']})
%     ylabel ('Moment')
%     xlabel ('Time')    

    figure
    semilogy(fraw_mom,Sraw_mom,f_HP1_mom,S_HP1_mom,'LineWidth',2), xlim([0 2]),grid on,  
    set(gca,'FontSize',14);
    legend({'Raw' ['HP 1st order filter fc =' num2str(HP_fc) ' Hz']})
    ylabel ('Moment')
    xlabel ('Frequency')
    
end

if plotNO1 == 1

    figure,
    bode(F_NO1,opts)
    xlim([0.001 10]),grid on, 
    legend({['Notch filter fn = ' num2str(fNO_s) ' Hz']})    
    
%     figure
%     plot (tout,moment,tout,Noch1_mom,'LineWidth',2),grid on
%     set(gca,'FontSize',14);
%     legend({'Raw' ['Notch filter fn =0,32' num2str(HP_fc) ' Hz']})
%     ylabel ('Moment')
%     xlabel ('Time')    

%     figure
%     semilogy(fraw_mom,Sraw_mom,f_NO1_mom,S_NO1_mom,'LineWidth',2), xlim([0 2]),grid on,  
%     set(gca,'FontSize',14);
%     legend({'Raw' ['Notch filter fn = ' num2str(fNO_s) ' Hz']})
%     ylabel ('Moment')
%     xlabel ('Frequency')

    figure
    plot (tout,pitch,tout,Noch1_pitch,'LineWidth',2),grid on
    set(gca,'FontSize',14);
    legend({'Raw' ['Notch filter fn = ' num2str(fNO_s) ' Hz']})
    ylabel ('Pitch')
    xlabel ('Time')    

    figure
    semilogy(fraw_pitch,Sraw_pitch,f_NO1_pitch,S_NO1_pitch,'LineWidth',2), xlim([0 2]),grid on,  
    set(gca,'FontSize',14);
    legend({'Raw' ['Notch filter fn = ' num2str(fNO_s) ' Hz']})
    ylabel ('Pitch')
    xlabel ('Frequency')

end

if plotNO3 == 1

    figure,
    bode(F_NO3,opts)
    xlim([0.001 10]),grid on, 
    legend({['Notch filter fn = ' num2str(fNO_s) ' Hz'] ['Notch filter fn = ' num2str(fNO_s) 'and ' num2str(fNO2_s) 'Hz']})    
    
    figure
    plot (tout,pitch,tout,Noch3_pitch,'LineWidth',2),grid on
    set(gca,'FontSize',14);
    legend({'Raw' ['Notch filter fn = ' num2str(fNO_s) 'and ' num2str(fNO2_s) 'Hz']})
    ylabel ('Pitch')
    xlabel ('Time')
    
    figure
    semilogy(fraw_mom,Sraw_mom,f_NO3_mom,S_NO3_mom,'LineWidth',2), xlim([0 2]),grid on,  
    set(gca,'FontSize',14);
    legend({'Raw' ['Notch filter fn = ' num2str(fNO_s) 'and ' num2str(fNO2_s) 'Hz']})
    ylabel ('Moment')
    xlabel ('Frequency')
    
    figure
    semilogy(fraw_pitch,Sraw_pitch,f_NO3_pitch,S_NO3_pitch,'LineWidth',2), xlim([0 2]),grid on,  
    set(gca,'FontSize',14);
    legend({'Raw' ['Notch filter fn = ' num2str(fNO_s) 'and ' num2str(fNO2_s) 'Hz']})
    ylabel ('Pitch')
    xlabel ('Frequency')    

end

if plotNO4 == 1

    figure,
    bode(FNO4_Dav,opts)
    xlim([0.001 10]),grid on, 
    legend({['Notch 3Param David filter fn = ' num2str(fNO4_s) ' Hz']})
    
%     figure
%     plot (tout,moment,tout,Noch1_mom,'LineWidth',2),grid on
%     set(gca,'FontSize',14);
%     legend({'Raw' ['Notch filter fn =0,32' num2str(HP_fc) ' Hz']})
%     ylabel ('Moment')
%     xlabel ('Time')    

%     figure
%     semilogy(fraw_mom,Sraw_mom,f_NO4_mom,S_NO4_mom,'LineWidth',2), xlim([0 2]),grid on,  
%     set(gca,'FontSize',14);
%     legend({'Raw' ['Notch filter fn = ' num2str(fNO4_s) ' Hz']})
%     ylabel ('Moment')
%     xlabel ('Frequency')

    figure
    plot (tout,pitch,tout,Noch1_pitch,'LineWidth',2),grid on
    set(gca,'FontSize',14);
    legend({'Raw' ['Notch 3Param David filter fn = ' num2str(fNO4_s) ' Hz']})
    ylabel ('Pitch')
    xlabel ('Time')    

    figure
    semilogy(fraw_pitch,Sraw_pitch,f_NO4_pitch,S_NO4_pitch,f_NO1_pitch,S_NO1_pitch,'LineWidth',2), xlim([0 2]),grid on,  
    set(gca,'FontSize',14);
    legend({'Raw' ['Notch 3Param David filter fn = ' num2str(fNO4_s) ' Hz']})
    ylabel ('Pitch')
    xlabel ('Frequency')

end

if plotNO5 == 1

    figure,
    bode(FNO5_Dav,opts)
    xlim([0.001 10]),grid on, 
    legend({['Notch 3Param David filter fn = ' num2str(fNO5_s) ' Hz']})
    
%     figure
%     plot (tout,moment,tout,Noch1_mom,'LineWidth',2),grid on
%     set(gca,'FontSize',14);
%     legend({'Raw' ['Notch filter fn =0,32' num2str(HP_fc) ' Hz']})
%     ylabel ('Moment')
%     xlabel ('Time')    

%     figure
%     semilogy(fraw_mom,Sraw_mom,f_NO4_mom,S_NO4_mom,'LineWidth',2), xlim([0 2]),grid on,  
%     set(gca,'FontSize',14);
%     legend({'Raw' ['Notch filter fn = ' num2str(fNO4_s) ' Hz']})
%     ylabel ('Moment')
%     xlabel ('Frequency')

    figure
    plot (tout,pitch,tout,Noch5_pitch,'LineWidth',2),grid on
    set(gca,'FontSize',14);
    legend({'Raw' ['Notch 3Param David filter fn = ' num2str(fNO5_s) ' Hz']})
    ylabel ('Pitch')
    xlabel ('Time')    

    figure
    semilogy(fraw_pitch,Sraw_pitch,f_NO5_pitch,S_NO5_pitch,'LineWidth',2), xlim([0 2]),grid on,  
    set(gca,'FontSize',14);
    legend({'Raw' ['Notch 3Param David filter fn = ' num2str(fNO5_s) ' Hz']})
    ylabel ('Pitch')
    xlabel ('Frequency')

end

if plotNO6 == 1

    figure,
    bode(FNO6_Dav,opts)
    xlim([0.001 10]),grid on, 
    legend({['Notch 3Param David filter fn = ' num2str(fNO5_s) ' Hz']})
    
%     figure
%     plot (tout,moment,tout,Noch1_mom,'LineWidth',2),grid on
%     set(gca,'FontSize',14);
%     legend({'Raw' ['Notch filter fn =0,32' num2str(HP_fc) ' Hz']})
%     ylabel ('Moment')
%     xlabel ('Time')    

%     figure
%     semilogy(fraw_mom,Sraw_mom,f_NO4_mom,S_NO4_mom,'LineWidth',2), xlim([0 2]),grid on,  
%     set(gca,'FontSize',14);
%     legend({'Raw' ['Notch filter fn = ' num2str(fNO4_s) ' Hz']})
%     ylabel ('Moment')
%     xlabel ('Frequency')

    figure
    plot (tout,pitch,tout,Noch6_pitch,'LineWidth',2),grid on
    set(gca,'FontSize',14);
    legend({'Raw' ['Notch 3Param David filter fn = ' num2str(fNO5_s) ' Hz']})
    ylabel ('Pitch')
    xlabel ('Time')    

    figure
    semilogy(fraw_pitch,Sraw_pitch,f_NO5_pitch,S_NO5_pitch,f_NO6_pitch,S_NO6_pitch,'LineWidth',2), xlim([0 2]),grid on,  
    set(gca,'FontSize',14);
    legend({'Raw' ['Notch 3Param David filter fn = ' num2str(fNO5_s) ' Hz']})
    ylabel ('Pitch')
    xlabel ('Frequency')

end

if plotNO42 == 1

    figure,
    bode(FNO42_Dav,opts)
    xlim([0.001 10]),grid on, 
    legend({['Notch 3Param David filter fn = ' num2str(fNO42_s) ' Hz']})
    
%     figure
%     plot (tout,moment,tout,Noch1_mom,'LineWidth',2),grid on
%     set(gca,'FontSize',14);
%     legend({'Raw' ['Notch filter fn =' num2str(fNO42_s) ' Hz']})
%     ylabel ('Moment')
%     xlabel ('Time')    

%     figure
%     semilogy(fraw_mom,Sraw_mom,f_NO4_mom,S_NO4_mom,'LineWidth',2), xlim([0 2]),grid on,  
%     set(gca,'FontSize',14);
%     legend({'Raw' ['Notch filter fn = ' num2str(fNO42_s) ' Hz']})
%     ylabel ('Moment')
%     xlabel ('Frequency')

    figure
    plot (tout,pitch,tout,Noch42_pitch,'LineWidth',2),grid on
    set(gca,'FontSize',14);
    legend({'Raw' ['Notch 3Param David filter fn = ' num2str(fNO42_s) ' Hz']})
    ylabel ('Pitch')
    xlabel ('Time')    

    figure
    semilogy(fraw_pitch,Sraw_pitch,f_NO42_pitch,S_NO42_pitch,f_NO3_pitch,S_NO3_pitch,'LineWidth',2), xlim([0 2]),grid on,  
    set(gca,'FontSize',14);
    legend({'Raw' ['Notch 3Param David filter fn = ' num2str(fNO4_s) ' Hz']})
    ylabel ('Pitch')
    xlabel ('Frequency')

end

if plotCOMB == 1
    
    figure,
    bode(F_COMB,opts)
    xlim([0.001 10]),grid on, 
    legend({['Notch filter fn = ' num2str(fNO_s) ' Hz']})    
    
    figure
    plot (tout,pitch,tout,COMB_pitch,'LineWidth',2),grid on
    set(gca,'FontSize',14);
    legend({'Raw' ['Combined Notch filter fn = ' num2str(fNO_s) ' and ' num2str(fNO2_s) 'Hz and in series LP with fc =' num2str(LP_fc_comb) 'Hz']})
    ylabel ('Pitch')
    xlabel ('Time')
    
    figure
    semilogy(fraw_mom,Sraw_mom,f_COMB_mom,S_COMB_mom,'LineWidth',2), xlim([0 2]),grid on,  
    set(gca,'FontSize',14);
    legend({'Raw' ['Combined Notch filter fn = ' num2str(fNO_s) ' and ' num2str(fNO2_s) 'Hz and in series LP with fc =' num2str(LP_fc_comb) 'Hz']})
    ylabel ('Moment')
    xlabel ('Frequency')
    
    figure
    semilogy(fraw_pitch,Sraw_pitch,f_COMB_pitch,S_COMB_pitch,f_NO3_pitch,S_NO3_pitch,'LineWidth',2), xlim([0 2]),grid on,  
    set(gca,'FontSize',14);
    legend({'Raw' ['Combined Notch filter fn = ' num2str(fNO_s) ' and ' num2str(fNO2_s) 'Hz and in series LP with fc =' num2str(LP_fc_comb) 'Hz'] 'Only double notch'})
    ylabel ('Pitch')
    xlabel ('Frequency')
end

if plotCOMB2 == 1
    
    figure,
    bode(F_COMB2,opts)
    xlim([0.001 10]),grid on, 
    legend({['Combined Notch David fn = ' num2str(fNO_s) ' and ' num2str(fNO2_s) 'Hz and in series LP with fc =' num2str(LP_fc_comb) 'Hz']})    
    
    figure
    plot (tout,pitch,tout,COMB2_pitch,tout,COMB_pitch,'LineWidth',2),grid on
    set(gca,'FontSize',14);
    legend({'Raw' ['Combined Notch David fn = ' num2str(fNO_s) ' and ' num2str(fNO2_s) 'Hz and in series LP with fc =' num2str(LP_fc_comb) 'Hz']})
    ylabel ('Pitch')
    xlabel ('Time')
    
%     figure
%     semilogy(fraw_mom,Sraw_mom,f_COMB2_mom,S_COMB2_mom,'LineWidth',2), xlim([0 2]),grid on,  
%     set(gca,'FontSize',14);
%     legend({'Raw' ['Combined Notch David fn = ' num2str(fNO_s) ' and ' num2str(fNO2_s) 'Hz and in series LP with fc =' num2str(LP_fc_comb) 'Hz']})
%     ylabel ('Moment')
%     xlabel ('Frequency')
    
    figure
    semilogy(fraw_pitch,Sraw_pitch,f_COMB2_pitch,S_COMB2_pitch,f_COMB_pitch,S_COMB_pitch,'LineWidth',2), xlim([0 2]),grid on,  
    set(gca,'FontSize',14);
    legend({'Raw' ['Combined  Davids Notch filter fn = ' num2str(fNO_s) ' and ' num2str(fNO2_s) 'Hz and in series LP with fc =' num2str(LP_fc_comb) 'Hz'] 'Double notch with LP'})
    ylabel ('Pitch')
    xlabel ('Frequency')
end

if plotBode==1
%     
%     figure,
%     bode(F_LP,opts);
%     xlim([0.0001 10]),grid on, 
%     legend({['LP 1st order filter fc =' num2str(LP_fc) ' Hz']})
%     
%     figure,
%     bode(F_HP,opts)
%     xlim([0.0001 10]),grid on, 
%     legend({['HP 1st order filter fc =' num2str(HP_fc) ' Hz']})
    
    figure,
    bode(F_NO1,opts)
    xlim([0.001 10]),grid on, 
    legend({['Notch filter fn = ' num2str(fNO_s) ' Hz']})
    
%     figure,
%     bode(F_NO3,opts)
%     xlim([0.001 10]),grid on, 
%     legend({['Notch filter fn = ' num2str(fNO_s) ' Hz']})
%     
%     figure,
%     bode(F_COMB,opts)
%     xlim([0.001 10]),grid on, 
%     legend({['Notch filter fn = ' num2str(fNO_s) ' Hz']})
%     
    figure,
    bode(FNO4_Dav,opts)
    xlim([0.001 10]),grid on, 
    legend({['Notch filter fn = ' num2str(fNO_s) ' Hz']})
    
end
