%Script to post process statistics,DEL and Spectra for DLC Dynamic inflow simulations

%VPE 21.03.2017

clear all, close all, clc %#ok<DUALC,CLALL>
addpath('D:\code_witlis_FAST\witlis\functions\tools')
run('D:\code_witlis_FAST\witlis\trunkSC\AddingWitlisPaths.m')

%% Explanation Statistics files

% FAST statisic files contain Channels in the following order: [7x11 --> speeds x channels]
  %1 v_0
  %2 theta
  %3 M_G
  %4 Omega rotor
  %5 X_T tower displacement
  %6 M_yt fore aft moment
  %7 pitch acceleration
  %8 P_el electric power
  %9 LSS torque
  %10 theta_c requested pitch angle
  %11 Angular acceleration rotor
  
% SLOW statisic files contain Channels in the following order: [7x11 --> speeds x channels]
%1 v_0
%2 theta_c
%3 M_g
%4 Omega
%5 P_el
%6 M_a
%7 x_T
%8 x_T_dot
%9 theta
%10 M_yT
%11 M_LSS

% Each of the statisitcs struct files has the following fields
%MEAN
%STD 
%MIN 
%MAX 
%S 
%f 
%DEL 

%% INPUTS
FASTstatPath={
            'D:\DynamicInflowSims\DLC12simulations\DTU10MW_FB_DLC12_FAST_Equil_Stea_NT\DTU10MW_FB_DLC12_FAST_Equil_Stea_NT_statistics.mat'
            'D:\DynamicInflowSims\DLC12simulations\DTU10MW_FB_DLC12_FAST_DynIn_Stea_NT\DTU10MW_FB_DLC12_FAST_DynIn_Stea_NT_statistics.mat'   
            'D:\DynamicInflowSims\DLC12simulations\DTU10MW_FB_DLC12_FAST_Equil_BLst_NT\DTU10MW_FB_DLC12_FAST_Equil_BLst_NT_statistics.mat'            
            'D:\DynamicInflowSims\DLC12simulations\DTU10MW_FB_DLC12_FAST_DynIn_BLst_NT\DTU10MW_FB_DLC12_FAST_DynIn_BLst_NT_statistics.mat'
            'D:\DynamicInflowSims\DLC12simulations\DTU10MW_FB_DLC12_FAST_Equil_Stea_NT_NoBl\DTU10MW_FB_DLC12_FAST_Equil_Stea_NT_NoBl_statistics.mat' 
            'D:\DynamicInflowSims\DLC12simulations\DTU10MW_FB_DLC12_FAST_DynIn_Stea_NT_NoBl\DTU10MW_FB_DLC12_FAST_DynIn_Stea_NT_NoBl_statistics.mat'      
            'D:\DynamicInflowSims\DLC12simulations\DTU10MW_FB_DLC12_FAST_Equil_BLst_NT_NoBl\DTU10MW_FB_DLC12_FAST_Equil_BLst_NT_NoBl_statistics.mat'            
            'D:\DynamicInflowSims\DLC12simulations\DTU10MW_FB_DLC12_FAST_DynIn_BLst_NT_NoBl\DTU10MW_FB_DLC12_FAST_DynIn_BLst_NT_NoBl_statistics.mat'    
            };
        
SLOWstatPath={
            'D:\DynamicInflowSims\DLC12simulations\DTU10MW_FBswe_DLC12_SLOW\DTU10MW_FBswe_DLC12_SLOW_statistics.mat'
             };

Statfields={'MEAN'
            'STD'
            'MIN'
            'MAX'
            'DEL'};
        
Freqfields={'S'                    
            'f'
            };
                
ChansFAST={ 
            'Vo'                   %Hard copied fileds expected in the correct order FAST
            'theta' 
            'M_g'   
            'Omega'
            'x_T' 
            'M_yT'
            'theta_dot'
            'P_el' 
            'M_LSS'
            'theta_c'
            'Omega_acc'    
            'MyBl'  
            'MxBl'                                                      
            'MzBl'    
            'MxTow' 
            'MzTow' 
            'Trot'
           };
       
ChansSLOW={ 
            'Vo'                    %Chans include in SLOW sims
            'theta_c' 
            'M_g'   
            'Omega'
            'P_el'
            'M_a'
            'x_T' 
            'x_Tdot'
            'theta'
            'M_yT' 
            'M_LSS'         
            };  
      
DegCell={
        'theta'           % cells tp be converted to deg
        'theta_dot'
        'theta_c'
       };

 MegaCell={
           'M_yT'          % chans to be converted in Mega...
           'P_el'
           'M_LSS'
           'M_a'
           };
       
RPMCell=  {
           'Omega'         % chans to be converted to RPM 
           'Omega_acc'
           };     
       
LegendCell ={
             'Equil Stea NT'                %Legend for the plots
             'DynIn Stea NT'  
             'Equil BLst NT'        
             'DynIn BLst NT'
             'Equil Stea NT NoBl'
             'DynIn Stea NT NoBl '
             'Equil BLst NT NoBl'   
             'DynIn BLst NT NoBl'
              'SLOW'
                };
            
plotStatfields = {
                  'theta'                  % Fields of Stat to be plotted
                  'M_g'
                  'Omega'
                  'x_T'
                  'M_yT'
                  'P_el'
                  'M_LSS'
%                   'MyBl'  
%                   'MxBl'                                                      
%                   'MzBl'    
%                   'MxTow' 
%                   'MzTow' 
%                   'Trot'
                  }; 
              
plotFreqfields ={ 
                  'theta'                  % Fields of Freq to be plotted
                  'M_g'
                  'Omega'
                  'M_yT'
                  'P_el'
                  'M_LSS'
%                   'MyBl'  
%                   'MxBl'                                                      
%                   'MzBl'    
%                   'MxTow' 
%                   'MzTow' 
%                   'Trot'
                    };   
              
 plotDELfields = {
                  'M_yT'                  % Fields of DELs to be plotted
                  'M_LSS'
%                   'MyBl'  
%                   'MxBl'                                                      
%                   'MzBl'    
%                   'MxTow' 
%                   'MzTow' 
                  };             
  
col_ord= [    0         0    1.0000       % Colour's order to be displayed
              0    0.5000         0
         1.0000         0         0
              0    0.7500    0.7500
         1.0000         0    1.0000     
         1.0000    0.6900    0.3900     
         0.7500         0    0.7500
         0.7500    0.7500         0
         0.2500    0.2500    0.2500];
     
MarkerCell={'o'                          % Line styles to be used
            '+'
            '*'
            'x'
            'd'
            '^'
            '>'
            '<'
            's'};      
        
PathToSave=  'D:\DynamicInflowSims\Comparison plots\DLC12\' ;       
                             
%% OPTIONS
        
Normalize      = 0 ;      % Choose to normalize data  
Baseline       = 1;       % Choose which simulation is the base
PlotStatistics = 1;       % Plot statistics
PlotFreqs      = 0;       % Plot Frequencies
SLOWinc        = 1;       % Flag to exclude SLOW for DOFs that are not included
SaveStatPlots  = 0;       % Save plots of Statistics
SaveFreqPlots  = 0;       % Save plots of Frequencies
speeds         = 12:2:24; % Speed vector to plot the statistics over
VisiblePlots   = 1 ;      % If 0 plots will not be dispalyed but still created and saved
%-------------------------------------------------------------------------%        
         
%% Load all files in one block 

%the order is defined in the previous section!!!

%Add FAST input
for SimNo = 1:length(FASTstatPath)
    load (FASTstatPath{SimNo})

    %Statitstics
    for StatField  = 1:length(Statfields) 
        for ChanNo = 1:length(ChansFAST)
           if ~strcmp(Statfields{StatField},'DEL') || (strcmp(Statfields{StatField},'DEL') && ChanNo<=size(Statistics.DEL,2))
                if SimNo == 1
                    StatsToT.(Statfields{StatField}).(ChansFAST{ChanNo}) (1,:) = [(Statistics.(Statfields{StatField})(:,ChanNo))']; %#ok<*NBRAK>
                else
                    StatsToT.(Statfields{StatField}).(ChansFAST{ChanNo})  = [StatsToT.(Statfields{StatField}).(ChansFAST{ChanNo}) ;(Statistics.(Statfields{StatField})(:,ChanNo))'];
                end
           end
        end
    end
% end  
% for SimNo=1:length(FASTstatPath)
  if PlotFreqs == 1   
    %Frequencies
    for FreqField  = 1:length(Freqfields) 
        for ChanNo = 1:length(ChansFAST)
%            if ~strcmp(Statfields{StatField},'DEL') || (strcmp(Statfields{StatField},'DEL') && ChanNo<=size(Statistics.DEL,2))
                if SimNo == 1
                    FreqToT.(Freqfields{FreqField}).(ChansFAST{ChanNo}) (1,:) = [(Statistics.(Freqfields{FreqField})(:,ChanNo))']; %#ok<*NBRAK>
                else
                    FreqToT.(Freqfields{FreqField}).(ChansFAST{ChanNo})  = [FreqToT.(Freqfields{FreqField}).(ChansFAST{ChanNo}) ;(Statistics.(Freqfields{FreqField})(:,ChanNo))'];
                end
%            end
        end
    end  
  end
  clearvars Statistics
end

%Add SLOW input
if SLOWinc == 1
for SimNo = 1:length(SLOWstatPath)
   load (SLOWstatPath{SimNo})
   
   %Statistics
   for StatField = 1:length(Statfields) 
        for ChanNo = 1:length(ChansSLOW)
           if ~strcmp(Statfields{StatField},'DEL') || ( strcmp(Statfields{StatField},'DEL') && (sum(strcmp(ChansSLOW{ChanNo},fieldnames(StatsToT.DEL)))~=0 ))
                if SimNo == 1 && (sum(strcmp(ChansSLOW{ChanNo},ChansFAST)) == 0)
                    StatsToT.(Statfields{StatField}).(ChansSLOW{ChanNo}) (1,:) = [(Statistics.(Statfields{StatField})(:,ChanNo))']; %#ok<*NBRAK>
                else
                    StatsToT.(Statfields{StatField}).(ChansSLOW{ChanNo})  = [StatsToT.(Statfields{StatField}).(ChansSLOW{ChanNo}) ;(Statistics.(Statfields{StatField})(:,ChanNo))'];
                end
           end
        end
   end

if PlotFreqs == 1   
   %Frequencies
   for FreqField = 1:length(Freqfields) 
        for ChanNo = 1:length(ChansSLOW)
%            if ~strcmp(Statfields{StatField},'DEL') || (strcmp(Statfields{StatField},'DEL') && ChanNo<=size(Statistics.DEL,2))
                if SimNo == 1 && (sum(strcmp(ChansSLOW{ChanNo},ChansFAST)) == 0)
                    FreqToT.(Freqfields{FreqField}).(ChansSLOW{ChanNo})(1,:) = [(Statistics.(Freqfields{FreqField})(:,ChanNo))']; %#ok<*NBRAK>
                else
                    FreqToT.(Freqfields{FreqField}).(ChansSLOW{ChanNo})  = [FreqToT.(Freqfields{FreqField}).(ChansSLOW{ChanNo}) ;(Statistics.(Freqfields{FreqField})(:,ChanNo))'];
                end
%            end
        end
   end  
end    
        clearvars Statistics
end

end

%% Convert all units appropriately

% for i = 1:length(fieldnames(StatsToT.MEAN))
%     for jj = 1:length(Statfields)
%         interfields = fieldnames(StatsToT.(Statfields{jj}));
%         if i<=length(interfields)  
%             if  ~strcmp (Statfields{jj},'DEL') || ( strcmp (Statfields{jj},'DEL') && ( strcmp(interfields{i},'M_yT')||strcmp(interfields{i},'M_LSS') )  )
%                 if any(strcmp(interfields{i},DegCell))
% %                    StatsToT.(Statfields{jj}).(interfields{i}) = rad2deg(StatsToT.(Statfields{jj}).(interfields{i})); 
%                 elseif any(strcmp(interfields{i},MegaCell))
%                     StatsToT.(Statfields{jj}).(interfields{i}) = (StatsToT.(Statfields{jj}).(interfields{i}))/1e6;  
%                 elseif any(strcmp(interfields{i},RPMCell))
%                     StatsToT.(Statfields{jj}).(interfields{i}) = radPs2rpm(StatsToT.(Statfields{jj}).(interfields{i}));      
%                 end
%             end
%         end
%     end
% end

%% Normalize statistics

if Normalize == 1     
interfields = fieldnames((StatsToT.MEAN))  ;  
for mm = 1:length(interfields) 
   for ii = 1:length(Statfields)-1 %assumes del is last 
    NormRow1{ii}(mm,:) = StatsToT.(Statfields{ii}).(interfields{mm})(Baseline,:);  %#ok<*SAGROW>
   end
    if mm <= length( fieldnames(StatsToT.DEL))
    DELRow1(mm,:)  = StatsToT.DEL.(interfields{mm})(Baseline,:);  %#ok<*SAGROW>
    end
end 
for ll = 1:length(Statfields)
     if sum(strcmp(Statfields{ll},'DEL')) ~=1
        for mm = 1:length(interfields) 
            NormRow = NormRow1{ll}(mm,:); 
            if size (StatsToT.(Statfields{ll}).(interfields{mm}),1) == (length(SLOWstatPath)+length(FASTstatPath))
               for  cnt1 = 1:size (StatsToT.(Statfields{ll}).(interfields{mm}),1)  
                  StatsToT.(Statfields{ll}).(interfields{mm})(cnt1,:) = RelComp( NormRow, StatsToT.(Statfields{ll}).(interfields{mm})(cnt1,:));
               end 
            elseif SLOWinc == 0  
               for  cnt1 = 1:size (StatsToT.(Statfields{ll}).(interfields{mm}),1)  
                  StatsToT.(Statfields{ll}).(interfields{mm})(cnt1,:) = RelComp(  NormRow,StatsToT.(Statfields{ll}).(interfields{mm})(cnt1,:));
               end  
            end
        end
     else
        for mm = 1:length(interfields)
            if mm<= length( fieldnames(StatsToT.DEL))
            DELRow = DELRow1(mm,:); 
            if size (StatsToT.(Statfields{ll}).(interfields{mm}),1) == (length(SLOWstatPath)+length(FASTstatPath))
               for  cnt1 = 1:size (StatsToT.(Statfields{ll}).(interfields{mm}),1)  
                  StatsToT.(Statfields{ll}).(interfields{mm})(cnt1,:) = StatsToT.(Statfields{ll}).(interfields{mm})(cnt1,:)./ DELRow;
               end
            elseif SLOWinc == 0     
               for  cnt1 = 1:size (StatsToT.(Statfields{ll}).(interfields{mm}),1)  
                  StatsToT.(Statfields{ll}).(interfields{mm})(cnt1,:) = StatsToT.(Statfields{ll}).(interfields{mm})(cnt1,:)./ DELRow;
               end  
            end
            end
        end         
     end
    
end
end

%% Plotting Statistics and DEL

if PlotStatistics == 1

if VisiblePlots == 0 
 set(0,'DefaultFigureVisible','off');    
end   

    for kk=1:length(plotStatfields)   

        figure('name',[plotStatfields{kk} ' M-M']);  
        subplot(3,1,1)
        plot1=plot (speeds,StatsToT.MEAN.(plotStatfields{kk}));
        title([plotStatfields{kk} ' Mean'])
        grid on
        legend(LegendCell,'Location','northwestoutside')    
        for i=1:length(plot1)
        set(plot1(i),'Color',col_ord(i,:),'LineStyle','none','Marker',MarkerCell{i},'markers',8);
        end
%         hold on
        subplot(3,1,2)
        plot1=plot (speeds,StatsToT.MIN.(plotStatfields{kk}));
        title([plotStatfields{kk} ' Min'])
        grid on
        for i=1:length(plot1)
        set(plot1(i),'Color',col_ord(i,:),'LineStyle','none','Marker',MarkerCell{i},'markers',8);
        end
%         hold on
        subplot(3,1,3)
        plot1=plot (speeds,StatsToT.MAX.(plotStatfields{kk}));
        title([plotStatfields{kk} ' Max'])
        grid on
        for i=1:length(plot1)
        set(plot1(i),'Color',col_ord(i,:),'LineStyle','none','Marker',MarkerCell{i},'markers',8);
        end
        set(gcf,'units','points','position',[60,60,700,450])
%         hold off              
        if SaveStatPlots == 1
              if SLOWinc == 1  
                  FigDir=[PathToSave 'Statistics\'];
              else
                  FigDir=[PathToSave 'Statistics\FASTonly\'];
              end
            if ~exist(FigDir, 'dir');    mkdir(FigDir);     end
            saveas(gcf,[FigDir plotStatfields{kk} '_MMM' '.fig'])
            saveas(gcf,[FigDir plotStatfields{kk} '_MMM' '.png'])
        end
        
        figure('name',[plotStatfields{kk} ' STD']);  
%         plot1=plot (speeds,StatsToT.MEAN.(plotStatfields{kk}));
        plot1=plot (speeds,(StatsToT.STD.(plotStatfields{kk})));
        legend(LegendCell)
                grid on
        for i=1:length(plot1)
        set(plot1(i),'Color',col_ord(i,:),'LineStyle','none','Marker',MarkerCell{i},'markers',8);
        end
        title([plotStatfields{kk} ' STD'])
        set(gcf,'units','points','position',[60,60,700,450])
        if SaveStatPlots == 1        
            if ~exist(FigDir, 'dir');    mkdir(FigDir);     end
            saveas(gcf,[FigDir plotStatfields{kk} '_STD' '.fig'])
            saveas(gcf,[FigDir plotStatfields{kk} '_STD' '.png'])
        end
        
    end
    

   for kk=1:length(plotDELfields) 
        figure('name',[plotDELfields{kk} ' DEL']);  
        plot1=plot (speeds,StatsToT.DEL.(plotDELfields{kk}));
        title([plotDELfields{kk} ' DEL'])
        grid on
        legend(LegendCell)
        for i=1:length(plot1)
            set(plot1(i),'Color',col_ord(i,:),'LineStyle','none','Marker',MarkerCell{i},'markers',8);
        end
        set(gcf,'units','points','position',[60,60,700,450])    
        if SaveStatPlots == 1    
            if ~exist(FigDir, 'dir');    mkdir(FigDir);     end
                saveas(gcf,[FigDir plotDELfields{kk} '_DEL' '.fig'])
                saveas(gcf,[FigDir plotDELfields{kk} '_DEL' '.png'])
        end             
   end
    
if VisiblePlots == 0 
 set(0,'DefaultFigureVisible','on');    
end

end

%% Plotting frequencies
if PlotFreqs == 1

if VisiblePlots == 0
 set(0,'DefaultFigureVisible','off');    
end

for kk=1:length(plotFreqfields)        
    for i=1:length(speeds)
        figure('name',['PSD ' plotFreqfields{kk} ' V ' num2str(speeds(i))]);  

        if SLOWinc == 1  
            for pp=1:(length(SLOWstatPath)+length(FASTstatPath))
              plot1=plot (FreqToT.f.(plotFreqfields{kk}){pp,i},FreqToT.S.(plotFreqfields{kk}){pp,i});    
              hold all
    %         for i=1:length(plot1)
    %         set(plot1(i),'Color',col_ord(i,:),'LineStyle','none','Marker',MarkerCell{i},'markers',8);
    %         end
    %         hold on

    %         hold off      
            end
        else
           for pp=1:(length(FASTstatPath))
              plot1=plot (FreqToT.f.(plotFreqfields{kk}){pp,i},FreqToT.S.(plotFreqfields{kk}){pp,i});    
              hold all
           end
        end
        title([plotFreqfields{kk} ' PSD V=' num2str(speeds(i))])
        hold off
        grid on
        legend(LegendCell,'Location','S')
        AXpr = gca;
        AXpr.XScale = 'log';
        AXpr.YScale = 'log';
        AXpr.YMinorGrid = 'off';        
        AXpr.XLim    = [1e-3 1e0];
        set(gcf,'units','points','position',[60,60,700,450])
        
if SaveFreqPlots==1
        if SLOWinc == 1
            FigDir2=[PathToSave 'Spectra\Vel_' num2str(speeds(i)) '\' ];
        else
            FigDir2=[PathToSave 'Spectra\FASTonly\Vel_' num2str(speeds(i)) '\' ];
        end
        if ~exist(FigDir2, 'dir');    mkdir(FigDir2);     end
        saveas(gcf,[FigDir2 plotFreqfields{kk} '_PSD_V' num2str(speeds(i)) '.fig'])
        saveas(gcf,[FigDir2 plotFreqfields{kk} '_PSD_V' num2str(speeds(i)) '.png'])
end
    
    end
end

if VisiblePlots == 0 
 set(0,'DefaultFigureVisible','on');
end

end

