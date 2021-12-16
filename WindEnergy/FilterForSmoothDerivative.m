% With this code we want to test different possibilities for filtering
% noisy data in order to get smooth derivative.
% VPE 28.11.2018

clc, clear all,close all %#ok<*CLALL>

load ('C:\Users\pettas\Desktop\TestDerFilt.mat')

%% Parameters
dt = (time(2)-time(1));
LP_fc_pos  = 0.1;   % LP filter corner frequency for position  
LP_fc_der  = 0.15 ; % LP filter corner frequency for derivative   

mov_av_pos = 300;   % number of points for moving average for position
mov_av_der = 40;   % number of points for moving average for derivative

sgf_ord_pos = 5;   % Savitzky-Golay Filter order of polunomial for position 
sgf_len_pos = 401; % Savitzky-Golay Filter window length for position (odd!)
sgf_ord_vel = 5;   % Savitzky-Golay Filter order of polunomial for velocity
sgf_len_vel = 10001; % Savitzky-Golay Filter window length for velocity (odd!)

%% Raw data 

% derivative time from raw
raw_der = diff(sig)/dt;                    %simple derivative dy/dx
raw_der = [raw_der(1,1);raw_der];          %correcting length
raw_cent_der = (sig(3:end)-sig(1:end-2))/(2*dt); %derivative with central differences
raw_cent_der = [raw_cent_der(1);raw_cent_der;raw_cent_der(end)];

% frequency domain
vWindow              = hamming(floor((length(time))/4/2)*2);
[S_raw,fraw]         = pwelch(detrend(sig,'constant'),vWindow,[],[],1/dt,'onesided'); 
[S_raw_der,fraw_der] = pwelch(detrend(raw_der,'constant'),vWindow,[],[],1/dt,'onesided'); 

%% Filtering LP
s = tf('s');

% LP for position
tauLP = 1/(LP_fc_pos*2*pi);    % s
F_LP  = 1 / (1 + s*tauLP); % TF
LP1_pos   = lsim(F_LP,sig,time);
[S_LP1_pos,f_LP1_pos]      = pwelch(detrend(LP1_pos,'constant'),vWindow,[],[],1/dt,'onesided'); 
LP1_der= diff(LP1_pos)/dt;
LP1_der =[LP1_der(1,1);LP1_der];
LP1_der_cent= (LP1_pos(3:end)-LP1_pos(1:end-2))/(2*dt);
LP1_der_cent = [LP1_der_cent(1);LP1_der_cent;LP1_der_cent(end)];
[S_LP1_vel,f_LP1_vel]      = pwelch(detrend(LP1_der,'constant'),vWindow,[],[],1/dt,'onesided'); 

% LP filter for derivative
tauLP_vel = 1/(LP_fc_der*2*pi);    % s
F_LP_vel=1 / (1 + s*tauLP_vel); % TF
LP1_vel      = lsim(F_LP_vel,LP1_der,time);
LP1_vel_cent = lsim(F_LP_vel,LP1_der_cent,time);

%% Savitzky-Golay Filter
SG_pos = sgolayfilt(sig,sgf_ord_pos,sgf_len_pos); %smoothened position
 
SG_der = diff(SG_pos)/dt;  % differentiate for deriative
SG_der =[SG_der(1);SG_der];
SG_vel = sgolayfilt(SG_der,sgf_ord_pos,sgf_len_pos);    %apply smoothening on the velocity
 
%% Moving average filter

%Filtering  position
b_mov = (1/mov_av_pos)*ones(1,mov_av_pos);
a_mov = 1;

Mov_pos     = filter(b_mov,a_mov,sig); % applying a simple movinf average filter on the position
Mov_pos_Delay = dt*(length(b_mov)-1)/2;     %  the delay introduced to the signal

Mov_der  = diff(Mov_pos)/dt;
Mov_der  = [Mov_der(1,1);Mov_der];

% Filter on the derivative
b_mov_der = (1/mov_av_der)*ones(1,mov_av_der);
a_mov_der = 1;
Mov_vel  = filter(b_mov_der,a_mov_der,Mov_der); % applying a simple movinf average filter on the derivative

%% Plotting

%%% plot raw data 
%Timeseries
% figure,
% subplot(2,1,1)
% plot(time,sig)
% title('Raw Data')
% ylabel ('Position', 'FontSize', 14)
% xlabel('time [s]','FontSize', 14)
% grid on
% subplot(2,1,2)
% plot(time,raw_der,time,raw_cent_der)
% ylabel ('Velocity', 'FontSize', 14)
% xlabel('time [s]','FontSize', 14)
% legend({'dy/dx' 'central differences'})
% grid on

%PSD
figure
subplot(2,1,1)
semilogy (fraw,S_raw,'Linewidth',2);
title('PSD Position raw')
xlabel('Frequency [Hz]','FontSize', 14)
ylabel ('PSD postion', 'FontSize', 14)
grid on
xlim([0 2])
subplot(2,1,2)
semilogy (fraw_der,S_raw_der,'Linewidth',2);
title('PSD Velocity ')
xlabel('Frequency [Hz]','FontSize', 14)
ylabel ('Velocity', 'FontSize', 14)
grid on
xlim([0 2])

%%% Plot filtered data
% Time series LP position
figure,
subplot(2,1,1)
plot(time,sig,time-1,LP1_pos,time-Mov_pos_Delay,Mov_pos,time,SG_pos,'Linewidth',2)
title('LP filtered Data')
ylabel ('Position', 'FontSize', 14)
xlabel('time [s]','FontSize', 14)
grid on
legend({'raw','LP1','Moving average+Delay','SG filter'},'FontSize', 14)
subplot(2,1,2)
plot(time,raw_der,time,LP1_der,time,LP1_vel,time-Mov_pos_Delay,Mov_der,time-Mov_pos_Delay,Mov_vel,time,SG_der,time,SG_vel,'Linewidth',2)
ylabel ('Velocity', 'FontSize', 14)
xlabel('time [s]','FontSize', 14)
grid on
legend({'raw','LP1 on pos','LP on pos+vel','Moving average on Pos','Moving average on Pos+Vel','SG filter on Pos','SG filter on Pos+Vel'},'FontSize', 14)

%PSD after LP
% figure
% subplot(2,1,1)
% semilogy (f_LP1_pos,S_LP1_pos,'Linewidth',2);
% title('PSD Position LP1')
% xlabel('Frequency [Hz]','FontSize', 14)
% ylabel ('PSD postion', 'FontSize', 14)
% grid on
% xlim([0 2])
% subplot(2,1,2)
% semilogy (f_LP1_vel,S_LP1_vel,'Linewidth',2);
% title('PSD Velocity LP1 ')
% xlabel('Frequency [Hz]','FontSize', 14)
% ylabel ('Velocity', 'FontSize', 14)
% grid on
% xlim([0 2])
