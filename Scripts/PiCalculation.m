clear all %#ok<*CLALL>
close all
clc

circpointsX=0:0.01:1;
circpointsY=(1-circpointsX.^2).^0.5;
for i=1:40
   
samples = 100*i;
Tot_points(i) = samples^2;

% gridded approach
[gridX,gridY]=meshgrid(linspace(0,1,samples),linspace(0,1,samples))  ;
rP=(gridX.^2+gridY.^2).^0.5;
k=find(rP<=1);
Circ=numel(k);

piCalc(i)= 4*Circ/Tot_points(i) ; %#ok<*SAGROW>

% figure
% plot(gridX,gridY,'o')
% hold on
% Xlim=[0 1];
% Ylim=[0 1];
% plot (circpointsX,circpointsY);
% title(['Points = ' num2str(Tot_points(i)) ' Pi = ' num2str(piCalc(i))])
% hold off

% random approach
% [gridXr,gridYr]=meshgrid(rand(1,samples),rand(1,samples))  ; % gird spacing with uniformly distributed random numbers

% completely random points from a uniform distribution
gridXr = rand(samples,samples);
gridYr = rand(samples,samples);  
rPr=(gridXr.^2+gridYr.^2).^0.5;
kr=find(rPr<=1);
Circr=numel(kr);
piCalcr(i)= 4*Circr/Tot_points(i) ;

% with normally distirubted nuumbers   DOESN't WORK!!!!
% X=randn(samples,samples);
% Y=randn(samples,samples);
% % [gridXr,gridYr]=meshgrid( abs(X./max(abs(X))),abs(Y./max(abs(Y))) )  ;
% gridXrn = abs(X./max(max(abs(X))));
% gridYrn = abs(Y./max(max(abs(Y))));  
% rPrn=(gridXrn.^2+gridYrn.^2).^0.5;
% krn=find(rPrn<=1);
% Circrn=numel(krn);
% piCalcrn(i)= 4*Circrn/Tot_points(i) ;


% figure
% plot(gridXr,gridYr,'o')
% hold on
% Xlim=[0 1];
% Ylim=[0 1];
% plot (circpointsX,circpointsY);
% title(['Points = ' num2str(Tot_points(i)) ' Pi = ' num2str(piCalc(i))])
% hold off

end


figure
plot(1:i,piCalc,'-X',1:i,piCalcr,'-o',1:i,piCalcrn,'-P')
xlabel('grid points')
set(gca,'xtick',1:i);
set(gca,'xticklabel',Tot_points);
% xticklabel(num2str(Tot_points))
ylabel('Pi value')
grid on
