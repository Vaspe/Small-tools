clearvars
clc
close all
years= 2009:2018;
% dir_path = 'C:\Users\pettas\Downloads\NEWAfiles\'; 
dir_path = 'D:\Dropbox\Papers and Presentations\WES 2021 Alpha Ventus Data and inter farm interactions\Review Process\NEWA_stuff\NEWAfiles\';


% Extract hourly wind speeds for all years.
for i = 1:length(years) % loop over years
    FilesAndFolders = dir([dir_path '*' num2str(years(i)) '*']);
    WS_years_hour =[];
    for ii= 1:numel(FilesAndFolders) %loop over days
        ncfile =[dir_path FilesAndFolders(ii).name];
%         ncdisp(ncfile) % show information of the file
        WS = ncread(ncfile,'WS') ;
        WSvec= squeeze(WS(10,10,1,:));
        WS_fin_day = mean(reshape(WSvec,2,[])); % resampling
        WS_years_hour = [WS_years_hour,WS_fin_day]; %#ok<*AGROW>
        clear WS WSvec WS_fin_day
    end
    WS_years_all(i,:) = WS_years_hour(1:24*364); %#ok<*SAGROW> sum up all hours in a year (365 day assumption)
    [N,edges] = histcounts(WS_years_all(i,:),[0.5:1:30.5]); %#ok<*NBRAK>
    Nfreq(i,:) = N/8760;  % get frequency of occurance ND
    clear WS_years_hour
end

% Statistics from time series
Stats_yearly.mean = mean(WS_years_all');
Stats_yearly.std = std(WS_years_all'); %#ok<*UDIM>
Stats_yearly.median = median(WS_years_all'); %#ok<*UDIM>

%% plotting
figure % plot freququency of occurence per year
for i = 1:size(WS_years_all,1)
    plot(1:30,100*Nfreq(i,:),'linewidth',2,'MarkerSize',6,'Marker','o')
    hold on 
end
aa= [num2str(years') repmat(' mean= ',size(years,2),1) num2str(Stats_yearly.mean',3)   ];
legend (aa)
grid on
xlabel('Wind speed bin center m/s')
ylabel('Frequecny of occurence %')
title(['Data for 1000m height MeanWsp =' num2str(mean(Stats_yearly.mean'),3) ' st.dev =  '  num2str(std(Stats_yearly.mean'),3) ])

figure % plot wind speed time series
for i = 1:size(WS_years_all,1)
    plot(WS_years_all(i,:))
    hold on
end
legend (num2str(years'))
