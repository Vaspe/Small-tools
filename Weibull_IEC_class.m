clear all,  clc, close all %#ok<*DUALC,*CLALL>

% Weibul (or Rayleigh) distributions Pdf for the different IEC classes.
% Bins can be defined by the input vector. Not much more than that....
%
% Vasilis Pettas 2.2018 SWE

%%
u=3:2:25;
k=2;
Vref=[50 42.5 37.5];

for j=1:length(Vref)
    Vave=0.2*Vref(j);
    A = exp(2);
    
    for i=1:length(u)
        cdf(j,i) = 1-exp(-(u(i)/A)^k); %#ok<*SAGROW>
        P_r(j,i) = 1-exp( -pi*(u(i)/(2*Vave))^k );
    end
    %
    % figure,bar(u,cdf)
    % figure,bar(u,P_r)
    bins=4:2:24;
    pdf(j,:)=diff(P_r(j,:)); 

end

figure,bar(bins,pdf(1,:)),grid on, legend({'PDF class I'})
% figure,bar(bins,pdf(2,:)),grid on, legend({'PDF class II'})
% figure,bar(bins,pdf(3,:)),grid on, legend({'PDF class III'})
figure,plot(bins,pdf(1,:),bins,pdf(3,:),bins,pdf(2,:)),grid on,legend({'Class I' 'Class II' 'Class III'})

PDFaboverated=P_r(1,end)-P_r(1,5);
sum(pdf(1,1:4));
sum(pdf(1,5:end));
PDFbelowrated=P_r(1,5)-P_r(1,1);

ClassIa_pdf =pdf(1,:);
