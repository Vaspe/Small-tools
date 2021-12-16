clearvars
clc
close all 

% Weibul (or Rayleigh) distributions Pdf for the different IEC classes.
% Bins can be defined by the input vector. Not much more than that....
%
% Vasilis Pettas 10.2020 SWE

%% Speed distribution
bins = 4:1:25; %velocities of bin centers in m/s
b = 2; % For Rayleigh distributions according to IEC and https://www.amazon.com/Renewable-Efficient-Electric-Power-Systems/dp/1118140621 and https://www.researchgate.net/publication/267839737_Rayleigh_Distribution-Based_Model_for_Prediction_of_Wind_Energy_Potential_of_Cameroon
Vref = [50 42.5 37.5]; % Class I II III
Iref = [0.18 0.16 0.14 0.12]; % Class A+ A B C

u = [bins(1)- (bins(2)-bins(1))/2, bins(1:end)+(bins(2)-bins(1))/2]; % velocity bin edges in m/s
for j = 1:length(Vref)
    Vave = 0.2*Vref(j);       
    CDF_Ray_u(j,:) = 1-exp(-(pi/4).*(u./Vave).^b); %#ok<*SAGROW>
    PDF_Ray_singleU(j,:) = (pi/2)*([u(1):1:u(end)]./(Vave^b)).*exp( -(pi/4).*([u(1):1:u(end)]/Vave).^b );  %#ok<*NBRAK> %according to IEC and https://www.amazon.com/Renewable-Efficient-Electric-Power-Systems/dp/1118140621 and https://www.researchgate.net/publication/267839737_Rayleigh_Distribution-Based_Model_for_Prediction_of_Wind_Energy_Potential_of_Cameroon
    
    PDF_bin(j,:) = diff(CDF_Ray_u(j,:));
end
% figure;plot(u,PDF_Ray_singleU(1,:),bins,PDF_bin(1,:))

figure,bar(bins,PDF_bin(1,:)),grid on,hold on,plot (bins,PDF_bin(1,:),'Linewidth',2), title('PDF class I'),xlabel('Wind Speeds m/s'),ylabel('Probability');
figure,bar(bins,PDF_bin(2,:)),grid on,hold on,plot (bins,PDF_bin(2,:),'Linewidth',2), title('PDF class II'),xlabel('Wind Speeds m/s'),ylabel('Probability');
figure,bar(bins,PDF_bin(3,:)),grid on,hold on,plot (bins,PDF_bin(3,:),'Linewidth',2), title('PDF class III'),xlabel('Wind Speeds m/s'),ylabel('Probability');

figure,plot(bins,PDF_bin(1,:),bins,PDF_bin(2,:),bins,PDF_bin(3,:),'Linewidth',2),grid on,legend({'Class I' 'Class II' 'Class III'}),xlabel('Wind Speeds m/s'),ylabel('PDF'),title('PDF all classes')
figure,plot(u,CDF_Ray_u(1,:),u,CDF_Ray_u(2,:),u,CDF_Ray_u(3,:),'Linewidth',2),grid on,legend({'Class I' 'Class II' 'Class III'}),xlabel('Wind Speeds m/s'),ylabel('CDF'),title('CDF all classes')

% PDFaboverated=PDF_Ray(1,end)-PDF_Ray(1,5);
% sum(pdf(1,1:4));
% sum(pdf(1,5:end));
% PDFbelowrated=PDF_Ray(1,5)-PDF_Ray(1,1);

%% Turbulece calcualted for the center of the bin

for j=1:length(Iref)
   sigma(j,:) = Iref(j)*(0.75.*bins+5.6); % IEC ed4
   TI(j,:) = sigma(j,:)./bins;    
end

% figure,plot(bins,sigma(1,:),bins,sigma(2,:),bins,sigma(3,:),bins,sigma(4,:),'Linewidth',2),grid on,legend({'Class A+' 'Class A' 'Class B' 'Class C'}), title('Standard deviation per turbulence class'),xlabel('Wind Speeds m/s'),ylabel('Standard deviation (m/s)')
figure,plot(bins,TI(1,:),bins,TI(2,:),bins,TI(3,:),bins,TI(4,:),'Linewidth',2),grid on,legend({'Class A+' 'Class A' 'Class B' 'Class C'}), title('Turbulence intensity per turbulence class'),xlabel('Wind Speeds m/s'),ylabel('Turbulence Intencity [-]')

%% In case of export...
IEC_class.ClassI.pdf(:,1) = bins' ;
IEC_class.ClassI.pdf(:,2) = PDF_bin(1,:)' ;

IEC_class.ClassII.pdf(:,1) = bins' ;
IEC_class.ClassII.pdf(:,2) = PDF_bin(2,:)' ;

IEC_class.ClassIII.pdf(:,1) = bins' ;
IEC_class.ClassIII.pdf(:,2) = PDF_bin(3,:)' ;



IEC_class.ClassI.cdf(:,1) = u' ;
IEC_class.ClassI.cdf(:,2) = CDF_Ray_u(1,:)' ;

IEC_class.ClassII.cdf(:,1) = u' ;
IEC_class.ClassII.cdf(:,2) = CDF_Ray_u(2,:)' ;

IEC_class.ClassIII.cdf(:,1) = u' ;
IEC_class.ClassIII.cdf(:,2) = CDF_Ray_u(3,:)' ;


IEC_class.ClassAplus(:,1) = bins' ;
IEC_class.ClassAplus(:,2) = TI(1,:)' ;

IEC_class.ClassA(:,1) = bins' ;
IEC_class.ClassA(:,2) = TI(2,:)' ;

IEC_class.ClassB(:,1) = bins' ;
IEC_class.ClassB(:,2) = TI(3,:)' ;

IEC_class.ClassC(:,1) = bins' ;
IEC_class.ClassC(:,2) = TI(4,:)' ;




