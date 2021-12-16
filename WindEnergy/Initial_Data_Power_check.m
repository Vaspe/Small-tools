clearvars
clc
close all

%% Define power curve and responses
% This section will be exchanged with simulaiton data from FAST
v = [0:28]'; %#ok<*NBRAK> %m/s velocity bin venters
% quick and dirty power curve
cp = 0.48;   % optimal Cp of turbine for calculation of below rated power
r = 178.3/2; % rotor radius [m]
rho = 1.225;  % air denstity kg/m3
V1 = 4:11 ;  % below rated wind speed range
P1=0.5*cp*rho*pi*(r^2)*(V1.^3); % below rated power
P_Base=[0,0,0, 0, P1, 10*10^6 10*10^6 10*10^6 10*10^6 10*10^6 10*10^6 10*10^6 10*10^6 10*10^6 10*10^6 10*10^6 10*10^6 10*10^6 10*10^6  0 0 0]' ; % full power curve including the rated power for above rated power
P_80 = P_Base *0.80; % Crude simplification to obtain changed power curves for different controllers
P_60 = P_Base *0.60; % Crude simplification to obtain changed power curves for different controllers
P_120 = [0,0,0,0, P1, 11*10^6 12*10^6 12*10^6 12*10^6 12*10^6 12*10^6 12*10^6 12*10^6 12*10^6 12*10^6 12*10^6 12*10^6 12*10^6 12*10^6  0 0 0]'; % Crude simplification to obtain changed power curves for different controllers
%  figure,plot(v,P_Base,v,P_80,v,P_120),grid on
% quick and dirty 
% DEL_Base= ;
% DEL_80 = ;
% DEL_120 = ;
P = [P_80 P_Base P_120]; % Add all power curves in one variable

%% Read  in power ans wind data
% this section needs to be updated to take in multiple years
% data here : D:\Dropbox\SWE\PhD Topic\Data Sets\Prices and wind FarmConners 2020 and same data 2030 hourly\Update_FarmConners_corres_results\Results
% FarmConners data
dataprice = readtable('D:\Dropbox\SWE\Tasks\FarmConners\Datasets\extratreat\farm_conners_DayAhead_Prices_2020.csv');
dataWSP = readtable('D:\Dropbox\SWE\Tasks\FarmConners\Datasets\extratreat\Results_CorRES\WSPD.csv');
dataUST = readtable('D:\Dropbox\SWE\Tasks\FarmConners\Datasets\extratreat\Results_CorRES\UST.csv');
% next two lines are adjustments for the pricer and wsp to match on time 
Price =  dataprice{:,2};
Price = [Price(1:24); Price; Price(end-23:end)];%eur/Mwh
WSP = dataWSP{:,2}; % wind speed in m/s
WSP(1) = WSP(2); % throw initial nan
UST = dataUST{:,2};
UST(1) = UST(2);
TI = 2.5*UST./WSP;

% time for plotting
t1 = datetime(2012,1,1,0,0,0);
t2 = datetime(2012,12,31,23,0,0);
t = [t1:hours(1):t2]'; % date time array

% Prices analysis
Price_m = mean(Price); %eur
Price_std= std(Price);
[N,edges] = histcounts(Price,0:2:60) ; % bin definition for prices
pdf_price = N/sum(N);
for i = 1:length(pdf_price)
    cdf_price(i) = sum(pdf_price(1:i));
end

% Wind analysis
[N1,edges1] = histcounts(WSP,0:2:30) ; % bin definition for wind speed
pdf_WSP = N1/sum(N1);
for i = 1:length(pdf_WSP)
    cdf_WSP(i) = sum(pdf_WSP(1:i));
end

%% Calculate produced energy and revenue
% Define index of baseline for comparison
Base_ind = 2; % find a better  way to do it...
% Precalculate power power time series from wind for each controller
for i =1:size(P,2)
    P_ofV(:,i) = interp1(v,P(:,i),WSP); %#ok<*SAGROW>
end
P_ofV=P_ofV/(10^6); %MW times series of power for all controllers

% Calculate Energy
for i =1:size(P,2)
    E_tot(:,i) = sum(P_ofV(:,i),'omitnan'); % MWh
end
% Calculate revenue
for i =1:size(P,2)
    Rev(:,i) = P_ofV(:,i).*Price;
    Rev_tot(i) = sum(Rev(:,i),'omitnan') ;
end

%% Basic optimization 1 threshold
thres_price1 = 14; % price threshold for single threshold price based optimization above this we boost
thres_DEL1 = 13; % price threshold for single threshold DEL based optimization below this we curtail

% Maximize revenue
[E_prof1,Rev_prof1,E_prof1_tot,Rev_prof1_tot,count_maxP] = price_thres1_cont1_opt(thres_price1,Price,P_ofV,2,3,'bigger',1.8);

% Minimize damage
[E_DEL,Rev_DEL1,E_DEL1_tot,Rev_DEL1_tot,count_minDEL] = price_thres1_cont1_opt(thres_DEL1,Price,P_ofV,2,1,'smaller',1.8);

%% Basic optimization 2 thresholds one low for curtasilment and one high for boosting
thres_low =  9; % low price threshold eur to curtail
thres_high = 18; % high price threshold eur to boost
[E_thres2,Rev_thres2,E_thres2_tot,Rev_thres2_tot,count_low1,count_up1] = price_thres2_cont2_opt(thres_low,thres_high,Price,P_ofV,1,3,2,1.8);

%% Basic optimization multiple thresholds

%% Add optimized values to the matrices
P_ofV = [P_ofV,E_prof1',E_DEL',E_thres2']; % MW or MWh
Delta_PofV = P_ofV./P_ofV(:,Base_ind); % relative power to baseline
Rev = [Rev,Rev_prof1',Rev_DEL1',Rev_thres2']; %Eur
E_tot = [E_tot,E_prof1_tot,E_DEL1_tot,E_thres2_tot]/10^3; % GWh
Rev_tot = [Rev_tot,Rev_prof1_tot,Rev_DEL1_tot,Rev_thres2_tot]/10^3; % 1000Eur re
Delta_E = E_tot./E_tot(:,Base_ind); % relative energy to base line
Delta_rev = Rev./Rev(:,Base_ind); % relative instantaneous revenue to baseline
Delta_rev_tot = Rev_tot./Rev_tot(Base_ind); % relative total revenue to baseline
bardata = [Delta_E',Delta_rev_tot']*100; % for the bar plot only

for ii =1:size(Rev,2)
    for i = 1:size(Rev,1)
        cum_rev(i,ii) = sum(Rev(1:i,ii))/10^3;  % thousand Eur Adding the cumulative revenue up to each time instance
    end
end

%% sensitivity on single thresholds
thres_range_single = 2:1:50; % variable space of price thresholds to cutrail or boost [eur]
for i =1:length(thres_range_single)
    [~,~,E_prof1_tot_sens(i),Rev_prof1_tot_sens(i),count_maxP_sens(i)] = price_thres1_cont1_opt(thres_range_single(i),Price,P_ofV,2,3,'bigger',1.8);
    [~,~,E_DEL1_tot_sens(i),Rev_DEL1_tot_sens(i),count_minDEL_sens(i)] = price_thres1_cont1_opt(thres_range_single(i),Price,P_ofV,2,1,'smaller',1.8);
end


%% sensitivity on double thrsholds
thres_low_range = [1:0.5:25]; % variable space of low price thresholds to cutrail [eur]
thres_up_range = [1:0.5:25]; % variable space of high price thresholds to boost [eur]
init_2thres_space = combvec(thres_low_range,thres_up_range);
init_2thres_space = init_2thres_space';
perm_vec_2thres=[];
for i=1:length(init_2thres_space)
    if init_2thres_space(i,1) < init_2thres_space(i,2) % values to boost must be higher than values to curtail
        perm_vec_2thres = [perm_vec_2thres;init_2thres_space(i,:)]; % final variable space for 2 thresholds
    end
end

for i =1:length(perm_vec_2thres)
    [~,~,E_2thres_tot(i),Rev_2thres_tot(i),count_low(i),count_up(i)] = price_thres2_cont2_opt(perm_vec_2thres(i,1),perm_vec_2thres(i,2),Price,P_ofV,1,3,Base_ind,1.8);
end

% find  values in a specific region of total revenue
rev_low_thres = Rev_tot(Base_ind)*1.00;
rev_up_thres = Rev_tot(Base_ind)*1.04;
Ind_opt_int = (Rev_2thres_tot/10^3>rev_low_thres) &  (Rev_2thres_tot/10^3<rev_up_thres);
Ind_rev_2thres_opt = find(Ind_opt_int);

% find values within a range of usage of each controller
use_per_low1 = 50;
use_per_high1 = 50;
Ind_use_low_int1 = find(100*count_low/length(t) < use_per_low1);
Ind_use_high_int1 = find(100*count_up/length(t) < use_per_high1);
ind_usage1 = intersect(Ind_use_low_int1,Ind_use_high_int1) ;

% sensitivity analysis on usage versus total revenue
use_per_low = 2.5:2.5:100; % max percentage of time the curtailment controller is on
use_per_high = 2.5:2.5:100;  % max percentage of time the boosting controller is on
for i=1:length(use_per_low)
    for ii=1:length(use_per_high)
        Ind_use_low_int = find(100*count_low/length(t) < use_per_low(i));
        Ind_use_high_int = find(100*count_up/length(t) < use_per_high(ii));
        ind_usage{i,ii} = intersect(Ind_use_low_int,Ind_use_high_int) ;
        if isempty(ind_usage{i,ii})
            Rev_usage_sense(i,ii) = nan;
        else
            Rev_usage_sense(i,ii) = max( Rev_2thres_tot(ind_usage{i,ii}));
        end
        clear Ind_use_low_int Ind_use_high_int 
    end
end

%% Plotting baseline and basic stuff

figure % price and wind trime series
yyaxis left
plot(table2array(dataWSP(:,1)),table2array(dataWSP(:,2)),'Linewidth',2)
ylabel('wind speed')
xlabel('time')
yyaxis right
plot(table2array(dataprice(:,1)),table2array(dataprice(:,2)),'Linewidth',2)
ylabel('price')
xlabel('time')
grid on
%

figure % bar plotcomparing cumulative revenue and power
bar(bardata)
legend('Enegy produced','Total revenue')
grid on
set(gca, 'XTickLabel', {'P80','P100','P120','EmaxCont','DelMinCont','2Thres'});
ylabel('% Difference')
set(gca, 'YMinorTick','on', 'YMinorGrid','on')


% cases running one controller all the time
figure % energy produced timeseries
subplot(2,1,1);
for i=1:size(P_ofV,2)
    plot(t,P_ofV(:,i),'linewidth',2)
    hold on
end
legend(['P80 Etot=' num2str(round(E_tot(1),2)) ' GWh'],['P100 Etot=' num2str(round(E_tot(2),2)) ' GWh'],...
    ['P120 Etot=' num2str(round(E_tot(3),2)) ' GWh'],['Emax Etot=' num2str(round(E_tot(4),2)) ' GWh'],...
    ['DELmin Etot=' num2str(round(E_tot(5),2)) ' GWh'],['2 Thres Etot=' num2str(round(E_tot(6),2)) ' GWh'])
title('Energy produced')
ylim ([0 20])
ylabel('Energy produced [MWh]')
xlabel('Time')
set(gca, 'YMinorTick','on', 'YMinorGrid','on')

subplot(2,1,2) % revenue time series
for i=1:size(P_ofV,2)
    plot(t,Rev(:,i),'linewidth',2)
    hold on
end
legend(['P80 Revtot=' num2str(round(Rev_tot(1))) ' thous Eur'],['P100 Revtot=' num2str(round(Rev_tot(2))) ' thous Eur'],...
    ['P120 Revtot=' num2str(round(Rev_tot(3))) ' thous Eur'],['EMax Revtot=' num2str(round(Rev_tot(4))) ' thous Eur'],...
    ['DELmIn Revtot=' num2str(round(Rev_tot(5))) ' thous Eur'],['2Thres Revtot=' num2str(round(Rev_tot(6))) ' thous Eur'])
title('Revenue')
ylabel('Revenue [eur]')
xlabel('Time')
set(gca, 'YMinorTick','on', 'YMinorGrid','on')

figure % Delta energy produced timeseries
subplot(2,1,1);
for i=1:size(Delta_PofV,2)
    plot(t,Delta_PofV(:,i),'linewidth',2)
    hold on
end
legend('P80','P100','P120',['EmaxCont ' num2str(round(count_maxP/length(WSP),2))],['DelMinCont ' num2str(round(count_minDEL/length(WSP),2))],['2Thres ' num2str(round(count_low1/length(WSP),2)) ',' ,num2str(round(count_up1/length(WSP),2))])
title('Delta energy produced')
ylabel('Relative to baseline [-]')
xlabel('Time')
set(gca, 'YMinorTick','on', 'YMinorGrid','on')

subplot(2,1,2); % Delta energy produced timeseries
for i=1:size(Delta_rev,2)
    plot(t,Delta_rev(:,i),'linewidth',2)
    hold on
end
legend('P80','P100','P120','EmaxCont','DelMinCont','2Thres')
title('Delta revenue')
ylabel('Relative to baseline [-]')
xlabel('Time')
set(gca, 'YMinorTick','on', 'YMinorGrid','on')

%% plotting optimization relevant
kk = 0;
for i=1:2:59
    kk=kk+1;
    xticklab_hist{kk}=num2str(i); %creates a vector of strings to be used for the x axis. Change the manual range to the relevant variable
end

figure  % price analysis
subplot(2,1,1);
bar(pdf_price)
xticks([1:30])
set(gca, 'XTickLabel', xticklab_hist);
ylabel('Probability of Occurence')
xlabel('MWh Price in eur')
grid on
subplot(2,1,2);
plot(1:2:59,cdf_price)
% xticks([1:30])
% set(gca, 'XTickLabel', xticklab_hist);
ylabel('CDF')
xlabel('MWh Price in eur')
grid on

nn = 0;
for i=1:2:29
    nn=nn+1;
    xticklab_hist2{nn}=num2str(i);
end

figure  % wind analysis
subplot(2,1,1);
bar(pdf_WSP)
xticks([1:15])
set(gca, 'XTickLabel', xticklab_hist2);
ylabel('Probability of Occurence')
xlabel('Wind Speed [m/s]')
grid on
subplot(2,1,2);
plot(1:2:29,cdf_WSP)
% xticks([1:30])
% set(gca, 'XTickLabel', xticklab_hist);
ylabel('CDF')
xlabel('Wind Speed [m/s]')
grid on

figure
for i=1:size(cum_rev,2)
    plot(t,cum_rev(:,i),'linewidth',2)
    hold on
end
legend(['P80 Revtot=' num2str(round(Rev_tot(1))) ' thous Eur'],['P100 Revtot=' num2str(round(Rev_tot(2))) ' thous Eur'],...
    ['P120 Revtot=' num2str(round(Rev_tot(3))) ' thous Eur'],['EMax Revtot=' num2str(round(Rev_tot(4))) ' thous Eur'],...
    ['DELmIn Revtot=' num2str(round(Rev_tot(5))) ' thous Eur'],['2thres Revtot=' num2str(round(Rev_tot(6))) ' thous Eur'],'Location','NW')
title('Revenue Accumulation')
ylabel('Total Revenue [1000eur]')
xlabel('Time')
set(gca, 'YMinorTick','on', 'YMinorGrid','on')

%% plotting sensitivity
figure  % price analysis
subplot(3,1,1);
plot(thres_range_single,(E_prof1_tot_sens/1000)/E_tot(Base_ind),thres_range_single,(E_DEL1_tot_sens/1000)/E_tot(Base_ind),'Linewidth',2)
ylabel('Total Energy produced GWh')
xlabel('Threshold in price eur/MWh')
grid on
legend('EMax','DELmIn')
subplot(3,1,2);
plot(thres_range_single,(Rev_prof1_tot_sens/1000)/Rev_tot(Base_ind),thres_range_single,(Rev_DEL1_tot_sens/1000)/Rev_tot(Base_ind),'Linewidth',2)
% xticks([1:30])
% set(gca, 'XTickLabel', xticklab_hist);
ylabel('Total revenue x1000eur')
xlabel('Threshold in price eur/MWh')
grid on
subplot(3,1,3);
plot(thres_range_single,100*count_maxP_sens/length(t),thres_range_single,100*count_minDEL_sens/length(t),'Linewidth',2)
% xticks([1:30])
% set(gca, 'XTickLabel', xticklab_hist);
ylabel('% of time controlled')
xlabel('Threshold in price eur/MWh')
grid on

figure % two threshold sensitivity
subplot(3,1,1)
scatter3(perm_vec_2thres(:,1),perm_vec_2thres(:,2),round((Rev_2thres_tot/10^3)/Rev_tot(Base_ind),2),'CData',round((Rev_2thres_tot/10^3)/Rev_tot(Base_ind),2))
colorbar
xlabel('Thres low')
ylabel('Thres up')
zlabel('Revenue relative to baseline')
grid on
title ('Total revenue with 2 thresholds (1000eur)')
view(gca,[0.710 90])
subplot(3,1,2)
scatter3(perm_vec_2thres(:,1),perm_vec_2thres(:,2),100*count_low/length(t),'CData',count_low/length(t))
colorbar
xlabel('Thres low')
ylabel('Thres up')
zlabel('Revenue x1000eur')
grid on
title ('Percentage of time down')
view(gca,[0.710 90])
subplot(3,1,3)
scatter3(perm_vec_2thres(:,1),perm_vec_2thres(:,2),100*count_up/length(t),'CData',count_up/length(t))
colorbar
xlabel('Thres low')
ylabel('Thres up')
zlabel('Revenue x1000eur')
grid on
title ('Percentage of time up')
view(gca,[0.710 90])

figure % two threshold sensitivity on specific ctiteria
subplot(2,1,1)
scatter3(perm_vec_2thres(Ind_rev_2thres_opt,1),perm_vec_2thres(Ind_rev_2thres_opt,2),round(Rev_2thres_tot(Ind_rev_2thres_opt)/10^3,2)/Rev_tot(Base_ind),'CData',round(Rev_2thres_tot(Ind_rev_2thres_opt)/10^3,2)/Rev_tot(Base_ind))
colorbar
title(['Revenue in range ' num2str(rev_low_thres) '-' num2str(rev_up_thres) 'thous eur'])
xlabel('Thres low')
ylabel('Thres up')
zlabel('Revenue x1000eur')
view(gca,[0.710 90])
grid on
subplot(2,1,2)
scatter3(perm_vec_2thres(ind_usage1,1),perm_vec_2thres(ind_usage1,2),round(Rev_2thres_tot(ind_usage1)/10^3,2)/Rev_tot(Base_ind),'CData',round(Rev_2thres_tot(ind_usage1)/10^3,2)/Rev_tot(Base_ind))
colorbar
title(['Revenue for usage <' num2str(use_per_low1) '% for low and <' num2str(use_per_high1) '% for high thres'])
xlabel('Thres low')
ylabel('Thres up')
zlabel('Revenue x1000eur')
view(gca,[0.710 90])
grid on
use_per_low1 = 30;
use_per_high1 = 30;



%% Functions
function [production,revenue,prod_tot,rev_tot,count_act] = price_thres1_cont1_opt(thres,inp1_thres,inp2,base_contr_ind,extr_contr_ind,flag_logic,shut_switch)
count_act =0;
for i = 1:length(inp1_thres)
    if inp1_thres(i) < shut_switch
        production(i) = 0;
        revenue(i) = 0;
    else
        if strcmp(flag_logic,'bigger')
            if inp1_thres(i)> thres
                production(i) = inp2(i,extr_contr_ind); %#ok<*AGROW>
                revenue(i) = inp2(i,extr_contr_ind)*inp1_thres(i);
                count_act = count_act +1;
            else
                production(i) = inp2(i,base_contr_ind);
                revenue(i) = inp2(i,base_contr_ind)*inp1_thres(i);
            end
        elseif strcmp(flag_logic,'smaller')
            if inp1_thres(i)<thres
                production(i) = inp2(i,extr_contr_ind); %#ok<*AGROW>
                revenue(i) = inp2(i,extr_contr_ind)*inp1_thres(i);
                count_act = count_act +1;
            else
                production(i) = inp2(i,base_contr_ind);
                revenue(i) = inp2(i,base_contr_ind)*inp1_thres(i);
            end
        end
    end
end
prod_tot = sum(production,'omitnan'); %GWh
rev_tot = sum(revenue,'omitnan'); %1000eur
end

function [production,revenue,prod_tot,rev_tot,count_act_low,count_act_up] = price_thres2_cont2_opt(thres1_low,thres2_up,inp1_thres,inp2,contr_ind_low,contr_ind_up,contr_ind_base,shut_switch)
count_act_low = 0;
count_act_up = 0;
for i = 1:length(inp1_thres)
    if inp1_thres(i) < shut_switch
        production(i) = 0;
        revenue(i) = 0;
    else
        if inp1_thres(i)> thres2_up
            production(i) = inp2(i,contr_ind_up); %#ok<*AGROW>
            revenue(i) = inp2(i,contr_ind_up)*inp1_thres(i);
            count_act_up = count_act_up +1;
        elseif inp1_thres(i)< thres1_low
            production(i) = inp2(i,contr_ind_low); %#ok<*AGROW>
            revenue(i) = inp2(i,contr_ind_low)*inp1_thres(i);
            count_act_low = count_act_low +1;
        else
            production(i) = inp2(i,contr_ind_base);
            revenue(i) = inp2(i,contr_ind_base)*inp1_thres(i);
        end
    end
end
prod_tot = sum(production,'omitnan'); %GWh
rev_tot = sum(revenue,'omitnan'); %1000eur
end

% load franke
% sf = fit([x, y],z,'poly23');
% plot(sf,[x,y],z)
