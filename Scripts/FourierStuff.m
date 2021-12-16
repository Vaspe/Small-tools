clear all
close all

% FFT and IFFT testing
%VPE 3.2020
% https://rosettacode.org/wiki/Fast_Fourier_transform#Perl
%% Inputs
Sig_mean = 10;
dt   = 0.1;
Tmax = 100;



%% Create signal
t   = [0:dt:Tmax]';
nt_init  = length(t);
sig_n = randn(nt_init,1);
% sig = Sig_mean+0.8*Sig_mean*sig_n+1.7*sin(2*pi*50*t);
sig = Sig_mean+1.7*sin(2*pi*2*t);


% for fft the signal's size has to be in the power of 2^n
% check and correct with zerro padding
if mod(log(nt_init)/log(2),1)~=0
    cur_pow     = log(nt_init)/log(2);
    Targ_pow    = ceil(cur_pow);
    targ_length = 2^Targ_pow;
    sig(end+1:targ_length) = 0;
    t(end+1:targ_length) = t(end)+[1:targ_length-nt_init]*dt;
end
nt  = length(t);
Fs = 1/max(t);
f_vec = Fs*(0:(nt/2));
figure,plot(t,sig),title ('Original signal'),xlabel('t [s]'),ylabel('y [sth]')



%% Do FFT with the Cooley–Tukey algorithm according to https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
tic
Odd_sig = sig(2:2:end);
Even_sig = sig(1:2:end);
k_cunt = 0;
for ik= 0:nt-1
    m_cunt = 0;
    k_cunt = k_cunt+1;
    for m=0:nt/2-1
        m_cunt = m_cunt+1;
        Exp_k(m_cunt,1) = exp( -2*pi*1i*m*ik/(nt/2));
        % E_k(i) = Odd_sig*exp( -2*pi*1i*m*i/(nt/2));
        % O_k = 1;
    end
    Exp_com = exp(-2*pi*1i*ik/nt);
    E_k(k_cunt,1) = sum(Odd_sig.*Exp_k);
    O_k(k_cunt,1) = sum(Even_sig.*Exp_k);
    X_k(k_cunt,1) = E_k(k_cunt,1) + Exp_com*O_k(k_cunt,1);
    X_kpN2(k_cunt,1) = E_k(k_cunt,1) - Exp_com*O_k(k_cunt,1);
end

t_fft_VP =toc;
disp(['Cooley Turkey implementation time = ' num2str(t_fft_VP) 's'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Geneeral fft
tic
cunt_gen = 0;
for i=0:nt-1
    cunt_gen =cunt_gen+1;
    X_k_gen(cunt_gen) = sum( sig.*exp(-2.*pi.*1i.*[0:nt-1]'.*i./(nt)) );
end
t_fft_general = toc;
figure, plot(f_vec,abs(X_k_gen(1:nt/2+1)/nt))
disp(['General iplementation time = ' num2str(t_fft_VP) 's'])

% matlab fft
tic
X_mat = fft(sig);
t_fft_mat = toc;
figure, plot(f_vec,abs(X_mat(1:nt/2+1)/nt))
disp(['Matlab fft implementation time = ' num2str(t_fft_mat) 's'])
abs(X_mat);

%% Do IFFT


% Do IFFT with random phase



