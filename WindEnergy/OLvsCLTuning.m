clear all, clc, close all %#ok<*DUALC,*CLALL>

%
% Work in progress script for plotting open loop and closed loop response
% od a system with a PID controller. Requires a plant model identified as
% TF or SS as input.
%
% TO DOs: automate open loop shaping by stability margins assesment
%
% Vasilis Pettas 2.2018 SWE


%%
Td= 0;
Kp= [1e-5 5e-5 1e-4 3e-4 5e-4];
for i = 1:length (Kp)
   Kpstr{i}=num2str(Kp(i)) ; %#ok<*SAGROW>
end
Ti= 0.3*Kp;
s = tf('s');
opts=  bodeoptions;
opts.FreqUnits = 'Hz';
% Plant= 1/(s*(s+1)*(s+5));% Ogata ex p572
% load ('D:\Tasks\Torque2018\Identification\models\BladeSysOrd5_LowTI_V16_Delta.mat');
load('D:\Tasks\Torque2018\Identification\models\BladeSysOrd4_LowTI_V16_Delta_TF.mat')
Plant= BladeSysOrd4_LowTI_V16_Delta_TF(1,1);
    
figure
for i =1:length(Kp)  
%     Control= Kp(i)*(1+(1/Ti(i))*(1/s)+Td*s); %standard PID
    Control= Kp(i); %standard PID
    tauHP   = 1/(0.05*2*pi);
     F_HP    = (tauHP*s ) / (1 + tauHP*s );
     OL=F_HP*Control*Plant;
   CL=Control*Plant/(1+Control*Plant); %#ok<NASGU>
    % allmargin(Plant)
    % figure,pzmap(Plant)
    % allmargin(OL)
    % figure,pzmap(OL)
    % allmargin(CL)
    % figure,pzmap(CL)
%     bode(Plant,'b',OL,'r',CL,'g',opts)
    bode(OL,opts)    
hold on
 clear OL CL Control
end
hold off
legend(Kpstr)

% figure
% for i =1:length(Kp)  
%         Control= Kp(i)*(1+(1/Ti(i))*(1/s)+Td*s); %standard PID
% %     Control= Kp(i); %standard PID
%     OL=Control*Plant;
%         CL=Control*Plant/(1+Control*Plant);
%     impulse(CL)    
% hold on
%  clear OL CL Control
% end
% hold off
% legend(Kpstr)
% 
% figure
% for i =1:length(Kp)  
%     Control= Kp(i)*(1+(1/Ti(i))*(1/s)+Td*s); %standard PID    
% %     Control= Kp(i); %standard PID
% %     OL=Control*Plant;
%     CL=Control*Plant/(1+Control*Plant);    
%     step(CL)    
% hold on
%  clear OL CL Control
% end
% legend(Kpstr)
% hold off
