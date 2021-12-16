%
% Script to identify desired trajectories for controller with different
% output levels. The input is an identified space of CP(TSR,theta) and a
% CT(TSR,theta).
%
% Desired levels of output need to be indicated by the user so that all the
% possible trajectories that follow this can be found.
%
% Vasilis Pettas 11/2020 University of Stuttgart
%
% To DOs
% -Update control logic to reflect realistic limits and switches
% -Automate so that so that some output of the identified feasible space is
%  pushed as input to the calculation of steady state for set point
%  identification

clc
clearvars
close all

%% Inputs
% Identification file with Cp,Ct,TSR,pitch (.mat file full path)
Ident_File = 'D:\Dropbox\SWE\PhD_Topic\ControllerAndSimulation\WitlisVP\trunkSC\DTU10MW_VP\PowerAndThrustCoefficientsDTU10MW';

% Turbine nominal characterstics
turbine.Ploss = 220; %(W)
turbine.el_eff = 0.94; % efficiency of the generator [-]
turbine.Pnom_el = 10*10^6; % nominal electrical power [W]
turbine.P_nom  = turbine.Ploss+turbine.Pnom_el*(1/turbine.el_eff);  % nominal aerodynamic power for 100% output (including losses and efficiency)
turbine.w_rat_nom = convangvel(9.6,'rpm','rad/s');   % rated rotor speed for 100% power in rpm
turbine.r = 89.2;                 % rotor radius [m]
turbine.min_w = convangvel(6,'rpm','rad/s'); % minimum rotor speed accepted in rpm
turbine.max_w = convangvel(9.6*1.2,'rpm','rad/s'); % maximum rotor speed accepted in rpm
turbine.V_nom = 11.4;   % rated wind speed m/s
turbine.min_theta_nom = 0;  % pitch for below rated in nominal operation in degrees
turbine.GB_ratio = 50; % GenSp/RotSp

% Get nominal values from inputs
turbine.Tq_nom = turbine.P_nom/(turbine.GB_ratio*turbine.w_rat_nom); % nominal torque (Nm)
turbine.Cp_nom = turbine.P_nom./(0.5.*1.225.*(turbine.V_nom.^3).*(pi.*turbine.r.^2)); % optimal cp required to reach nominal 100% rating (10MW)
turbine.TSR_nom = (turbine.w_rat_nom.*turbine.r)./turbine.V_nom;    % optimal TSR ar rated rotor speed for nominal 100% rating (10MW)

% Options
options.TSR_discr = 0.025; % 0.01  % discretization of resampled data for TSR (the identification was done with 0.1 steps)
options.theta_discr = 0.025; %0.0125; % discretization of resampled data for theta (the identification was done with 0.5 deg steps)
options.V     = 4:0.05:26;   % wsp vector to calculate the steady state values m/s
options.Power_factor_vec = 0.5:0.05:1.3; % multiplier of the nominal power. This applies to the whole wind pseed region
options.Cp_tolerance = 0.0002; % tolerance to search for feasible TSR and theta for a target value of Cp
options.w_min = 6; %RPM minimum allowed rotor speed for all cases
options.Vvec = 6:0.1:17; % possible rated wind speeds lower than nominal to find transition regions. Not meaningfull to change it
options.TSR_targ_vec = 5.5:0.5:10.5;  %^TSR targeted for geometrical trajectories
restriction.TSR = [5.5 11]; % [min max]
restriction.theta = [-4 12]; % [min max] RPM
options.control_names = {'minCT' 'constθ' 'constTSR' 'lin55' 'lin60' 'lin65' 'lin70' 'lin75' 'lin80' 'lin85' 'lin90' 'lin95' 'lin100' 'lin105'} ; % just for plotting
options.plot.feasible_restr_space_case = 12; % choose the index of the corresponding value in options.Power_factor_vec for CP CT space plotting
options.flag.steady_states_all_rates_case  = 5; % corresponds to the trajectory to be plotted. Their order is defined in options.control_names
options.steady_states_compare_case = 10; % choose the index of the corresponding value in options.Power_factor_vec for comparison of control strategies in steady states

% Flags for plotting
flag.plot_coefficents_space_all = 0; % Do you want to see the coefficent space
flag.plot_steady_states_all_rates = 1;  % plot the steady state curves for one of the designed controllers at all rated values
flag.plot_steady_states_compare = 1;  % plot the steady state for all the designed controllers at one rated value
flag.plot_feasible_restr_space = 0; % plotting the restricted space for a specific target Cp target value
flag.plot_traj = 0; % plot trajectories in the feasible space

%% Create data space

% Get identified values
DataIn = load(Ident_File); % load Cp-Ct (TSR,theta) identified tables from FAST
[IdDat.theta,IdDat.lambda,IdDat.C_T,IdDat.C_P,IdDat.C_M] = InterpAndClean(rad2deg(DataIn.theta),DataIn.lambda,DataIn.c_T,DataIn.c_P,options.TSR_discr,options.theta_discr);
[IdDat.theta2D,IdDat.lambda2D] = meshgrid(IdDat.theta,IdDat.lambda);

%% Calculation of feasible space
% Get all feasible space for all levels of power output requested

%Get target Cp from power factor
P_rat_el_vec = options.Power_factor_vec*turbine.Pnom_el;
P_rat_aero_vec = turbine.Ploss+P_rat_el_vec*(1/turbine.el_eff);

for ii = 1:length(P_rat_aero_vec)
    for i = 1:length(options.Vvec)
        Cp_targ_all_vec(ii,i) =  P_rat_aero_vec(ii)/(0.5.*1.225.*(options.Vvec(i).^3)*(pi*turbine.r.^2));
    end
    [~,pos] = sort(abs(Cp_targ_all_vec(ii,:)-turbine.Cp_nom*options.Power_factor_vec(ii)));
    Cp_targ(ii) = Cp_targ_all_vec(ii,pos(1)); % tricky point I take the slightly higher Cp check 1 or 2...
    if Cp_targ(ii)>turbine.Cp_nom
        Cp_targ(ii) = turbine.Cp_nom;
    end
end

for i =1:length(options.Power_factor_vec)
    %     [feasib.lambda_vec{i},feasib.theta_vec{i},feasib.CP_out{i},feasib.CT_out{i},feasib.CM_out{i}] = get_Cp_space(IdDat.C_P,round(options.Power_factor_vec(i)*turbine.Cp_nom,3)-options.Cp_tolerance,round(options.Power_factor_vec(i)*turbine.Cp_nom,3)+options.Cp_tolerance,IdDat.C_T,IdDat.C_M,IdDat.lambda,IdDat.theta);
    [feasib.lambda_vec{i},feasib.theta_vec{i},feasib.CP_out{i},feasib.CT_out{i},feasib.CM_out{i}] = get_Cp_space(IdDat.C_P,Cp_targ(i)-options.Cp_tolerance,Cp_targ(i)+options.Cp_tolerance,IdDat.C_T,IdDat.C_M,IdDat.lambda,IdDat.theta);
end

%% Choose set points from the feasible space

% Restrictions:
% TSR has to betweem 5.5 and 10.5 outside this the speeds are really weird
% Min pitch should be below 10 degrees so there is reserve before the stall

for i =1:length(options.Power_factor_vec)
    
    % reduce feasible space according to restrictions
    feasib_restr = get_restricted_space(feasib.CT_out{i},feasib.CP_out{i},feasib.CM_out{i},feasib.lambda_vec{i},feasib.theta_vec{i},restriction);
    feasib_restr2{i} = feasib_restr; % saving it for plotting. Remove later
    
    
    % min CT find the minimum CT with restrictions %%%%%%%%
    if ~isempty(feasib_restr.CT_out)
        [~,minCTpos] = min(feasib_restr.CT_out);
        Traject.minCT.CT(i)  = feasib_restr.CT_out(minCTpos);
        Traject.minCT.CP(i)  = feasib_restr.CP_out(minCTpos);
        Traject.minCT.CM(i)  = feasib_restr.CM_out(minCTpos);
        Traject.minCT.TSR(i) = feasib_restr.lambda_vec(minCTpos);
        Traject.minCT.theta(i) = feasib_restr.theta_vec(minCTpos);
        Traject.minCT.kOmega(i) = 0.5*1.225*pi*(turbine.r^5)*Traject.minCT.CP(i)/ ((Traject.minCT.TSR(i)^3)*(turbine.GB_ratio^3));
    else
        Traject.minCT.CT(i)  = [];
        Traject.minCT.CP(i)  = [];
        Traject.minCT.CM(i)  = [];
        Traject.minCT.TSR(i) = [];
        Traject.minCT.theta(i) = [];
        Traject.minCT.kOmega(i) = [];
    end
    
    
    %%%% Only TSR: Keep min pitch angle (almost) constant and change only TSR
    if ~isempty(feasib_restr.CT_out)
        pos_ind_theta = find(feasib_restr.theta_vec==turbine.min_theta_nom);
        if i==1
            [~,pos_theta] = min(abs(feasib_restr.theta_vec-turbine.min_theta_nom)) ;
        else
            if size(pos_ind_theta)==1 % if one point just keep it
                pos_theta = pos_ind_theta;
            else  % if more than one point, take the closest value. Pay attention to the manual tolerance for point definition
                pos_theta = get_disamb_val(feasib_restr.theta_vec,turbine.min_theta_nom,feasib_restr.lambda_vec,0.005,Traject.theta.theta,Traject.theta.TSR);
            end
        end
        Traject.TSR.CT(i)  = feasib_restr.CT_out(pos_theta);
        Traject.TSR.CP(i)  = feasib_restr.CP_out(pos_theta);
        Traject.TSR.CM(i)  = feasib_restr.CM_out(pos_theta);
        Traject.TSR.TSR(i) = feasib_restr.lambda_vec(pos_theta);
        Traject.TSR.theta(i) = feasib_restr.theta_vec(pos_theta);
        Traject.TSR.kOmega(i) = 0.5*1.225*pi*(turbine.r^5)*Traject.TSR.CP(i)/ ((Traject.TSR.TSR(i)^3)*(turbine.GB_ratio^3));
    else
        Traject.TSR.CT(i)  = [];
        Traject.TSR.CP(i)  = [];
        Traject.TSR.CM(i)  = [];
        Traject.TSR.TSR(i) = [];
        Traject.TSR.theta(i) = [];
        Traject.TSR.kOmega(i) = [];
    end
    
    
    %%%%%%% Only pitch:  keep TSR steady to the nominal
    if ~isempty(feasib_restr.CT_out)
        pos_ind_TSR = find(feasib_restr.lambda_vec==turbine.TSR_nom);
        if size(pos_ind_TSR)==1 % if one point just keep it
            pos_TSR = pos_ind_TSR;
        else  % if more than one point or no point found, take the closest value
            if i==1
                [~,pos_TSR] = min(abs(feasib_restr.lambda_vec-turbine.TSR_nom));
            else
                pos_TSR = get_disamb_val(feasib_restr.lambda_vec,turbine.TSR_nom,feasib_restr.theta_vec,0.005,Traject.theta.TSR,Traject.theta.theta); % check the manual tolerance if needed
            end
        end
        Traject.theta.CT(i)  = feasib_restr.CT_out(pos_TSR);
        Traject.theta.CP(i)  = feasib_restr.CP_out(pos_TSR);
        Traject.theta.CM(i)  = feasib_restr.CM_out(pos_TSR);
        Traject.theta.TSR(i) = feasib_restr.lambda_vec(pos_TSR);
        Traject.theta.theta(i) = feasib_restr.theta_vec(pos_TSR);
        Traject.theta.kOmega(i) = 0.5*1.225*pi*(turbine.r^5)*Traject.theta.CP(i)/ ((Traject.theta.TSR(i)^3)*(turbine.GB_ratio^3));
    else
        Traject.theta.CT(i)  = [];
        Traject.theta.CP(i)  = [];
        Traject.theta.CM(i)  = [];
        Traject.theta.TSR(i) = [];
        Traject.theta.theta(i) = [];
        Traject.theta.kOmega(i) = [];
    end
    
    
    %%%%% Combination: a bit of pitch and a bit of TSR based on geometry...
    if 1==~isempty(feasib_restr.lambda_vec)  % just to be able to fold it
        
        % Define 2 points and line based on target TSRs
        if i==1
            lin.P1 = [turbine.min_theta_nom turbine.TSR_nom];
            for ii=1:numel(options.TSR_targ_vec)
                [~,TSR_ind_find] = min(abs(feasib_restr.lambda_vec-options.TSR_targ_vec(ii)));
                % find second point based on target TSR of the minimum Cptarg
                lin.P2{ii} = [feasib_restr.theta_vec(TSR_ind_find) feasib_restr.lambda_vec(TSR_ind_find)];
                lin.x_lin{ii} = 0:0.01:restriction.theta(2);  % we re in the positive part of the pitch axis. This helps to make the output monotonic
                % create all coordinate pairs with a specific discretization
                lin.y_lin{ii} = ((lin.x_lin{ii}-lin.P1(1)).*(lin.P2{ii}(2)-lin.P1(2))./(lin.P2{ii}(1)-lin.P1(1)))+lin.P1(2);
                % find all Cps corresponding to the line
                for i_x = 1:numel(lin.x_lin{ii})
                    lin.CP_xy{ii}(i_x) = GetCoefficient(lin.x_lin{ii}(i_x),lin.y_lin{ii}(i_x),IdDat.theta,IdDat.lambda,IdDat.C_P)   ;
                end
            end
        end
        
        %Calculate set point crossing through the line and the cp
        for ii=1:numel(options.TSR_targ_vec)
            [~,Cp_ind] = min(abs(lin.CP_xy{ii}-Cp_targ(i)));
            %              RESPECT TO FIONA!
            %   [row,col]=find(feasib_restr.theta_vec>9.8375-0.01&X<9.8375+0.01&Y<6.87+0.01&Y>6.87-0.01&(abs(IdDat.C_P-0.23451)==min(abs(IdDat.C_P-0.23451))))
            %   IdDat.C_P(row,col)
            
            Traject.GEO{ii}.CT(i)  = GetCoefficient(lin.x_lin{ii}(Cp_ind),lin.y_lin{ii}(Cp_ind),IdDat.theta,IdDat.lambda,IdDat.C_T)  ;
            Traject.GEO{ii}.CP(i)  = GetCoefficient(lin.x_lin{ii}(Cp_ind),lin.y_lin{ii}(Cp_ind),IdDat.theta,IdDat.lambda,IdDat.C_P)  ;
            Traject.GEO{ii}.CM(i)  = GetCoefficient(lin.x_lin{ii}(Cp_ind),lin.y_lin{ii}(Cp_ind),IdDat.theta,IdDat.lambda,IdDat.C_M)  ;
            Traject.GEO{ii}.TSR(i) = lin.y_lin{ii}(Cp_ind);
            Traject.GEO{ii}.theta(i) = lin.x_lin{ii}(Cp_ind);
            Traject.GEO{ii}.kOmega(i) = 0.5*1.225*pi*(turbine.r^5)*Traject.GEO{ii}.CP(i)/ ((Traject.GEO{ii}.TSR(i)^3)*(turbine.GB_ratio^3));
        end
        
    else
        for ii=1:numel(options.TSR_targ_vec)
            Traject.GEO{ii}.CT(i)  = [];
            Traject.GEO{ii}.CP(i)  = [];
            Traject.GEO{ii}.CM(i)  = [];
            Traject.GEO{ii}.TSR(i) = [];
            Traject.GEO{ii}.theta(i) = [];
            Traject.GEO{ii}.kOmega(i) = [];
        end
    end
    
end

%% Calculation of steady state set points

% Manual input from looking at the previous values. Automate this somehow!
namesf = fieldnames(Traject);
for i_traj = 1:length(namesf)
    Traj_to_calc = Traject.(namesf{i_traj});
    if ~strcmp(namesf{i_traj},'GEO')
        TSR_val  = [turbine.TSR_nom   Traj_to_calc.TSR  ]; %7.8;
        theta_min_val = [turbine.min_theta_nom  Traj_to_calc.theta]; %0;
        Pnom_val = [1*turbine.P_nom P_rat_aero_vec]; %#ok<*NBRAK>
        
        for i =1:length(TSR_val)
            Res = GetOperChar(TSR_val(i),theta_min_val(i),Pnom_val(i),options.w_min,turbine.r,IdDat.theta,IdDat.lambda,IdDat.C_P,IdDat.C_T,options.V,turbine.GB_ratio,turbine.Ploss,turbine.el_eff);
            SS{i_traj,i} = Res;
            clear Res
        end

    else
        for i_geo = 1:length(Traj_to_calc)
            TSR_val  = [turbine.TSR_nom   Traj_to_calc{i_geo}.TSR]; %7.8;
            theta_min_val = [turbine.min_theta_nom  Traj_to_calc{i_geo}.theta]; %0;
            Pnom_val = [1*turbine.P_nom P_rat_aero_vec]; %#ok<*NBRAK>
            cnt = size(SS,1);
            for i =1:length(TSR_val)
                Res = GetOperChar(TSR_val(i),theta_min_val(i),Pnom_val(i),options.w_min,turbine.r,IdDat.theta,IdDat.lambda,IdDat.C_P,IdDat.C_T,options.V,turbine.GB_ratio,turbine.Ploss,turbine.el_eff);
                SS{cnt+1,i} = Res;
                clear Res
            end
        end
    end
end

%% Plotting

if flag.plot_coefficents_space_all 
    % C_P,CT,CM surfaces
    
    % figure,surf(IdDat.theta2D,IdDat.lambda2D,IdDat.C_P,IdDat.C_P)
    % xlabel('theta')
    % ylabel('lambda')
    % zlabel('Cp')
    % set(gca,'CLim',[0 0.5]);
    % colorbar
    
    figure,contour(IdDat.theta2D,IdDat.lambda2D,IdDat.C_P,[0:0.025:0.475 max(max(IdDat.C_P))],'ShowText','on')
    xlabel('theta')
    ylabel('lambda')
    xlim([restriction.theta(1) restriction.theta(2)])
    ylim([[restriction.TSR(1) restriction.TSR(2)]])
    title('Cp')
    
    
    % figure,surf(IdDat.theta2D,IdDat.lambda2D,C_T,IdDat.C_T)
    % xlabel('theta')
    % ylabel('lambda')
    % zlabel('Ct')
    % zlim(gca,[0 3]);
    % set(gca,'CLim',[0 3]);
    % colorbar
    
    figure,contour(IdDat.theta2D,IdDat.lambda2D,IdDat.C_T,0.1:0.05:2.0,'ShowText','on')
    xlabel('theta')
    ylabel('lambda')
    xlim([restriction.theta(1) restriction.theta(2)])
    ylim([[restriction.TSR(1) restriction.TSR(2)]])
    title('Ct')
    
    
    % figure,surf(IdDat.theta2D,IdDat.lambda2D,IdDat.C_M,IdDat.C_M)
    % xlabel('theta')
    % ylabel('lambda')
    % zlabel('Cm')
    % zlim(gca,[0 0.1]);
    % set(gca,'CLim',[0 0.1]);
    % colorbar
end

if flag.plot_feasible_restr_space
    % PLoting specific feasible space think of a smarter way to automate!
    plot_case = options.plot.feasible_restr_space_case;
    
    figure
    subplot(1,2,1)
    scatter3(feasib_restr2{plot_case}.theta_vec,feasib_restr2{plot_case}.lambda_vec,feasib_restr2{plot_case}.CP_out)
    view(0,90);
    ylabel('lambda')
    xlabel('theta')
    xlim([restriction.theta(1) restriction.theta(2)])
    ylim([[restriction.TSR(1) restriction.TSR(2)]])
    zlabel('Cp')
    grid on
    title(['Cp for Target Cp = ' num2str(round(Cp_targ(plot_case),4)) ' Pel = ' num2str(round(P_rat_el_vec(plot_case)/10^6,2)) ' MW'])
    subplot(1,2,2)
    scatter3(feasib_restr2{plot_case}.theta_vec,feasib_restr2{plot_case}.lambda_vec,feasib_restr2{plot_case}.CT_out,[],feasib_restr2{plot_case}.CT_out)
    ylabel('lambda')
    xlabel('theta')
    xlim([restriction.theta(1) restriction.theta(2)])
    ylim([[restriction.TSR(1) restriction.TSR(2)]])
    zlabel('Ct')
    grid on
    view(0,90);
    colorbar
    title(['Ct for Target Cp = ' num2str(round(Cp_targ(plot_case),4))])
end

if flag.plot_traj
    
    figure
    subplot(1,2,1)
    contour(IdDat.theta2D,IdDat.lambda2D,IdDat.C_P,[Cp_targ],'ShowText','on')
    xlabel('theta')
    ylabel('lambda')
    xlim([restriction.theta(1) restriction.theta(2)])
    ylim([[restriction.TSR(1) restriction.TSR(2)]])
    title('Cp')
    hold on
    plot(Traject.TSR.theta,Traject.TSR.TSR,'O-','MarkerSize',10)
    hold on
    plot(Traject.theta.theta,Traject.theta.TSR,'X-','MarkerSize',10)
    hold on
    plot(Traject.minCT.theta,Traject.minCT.TSR,'^-','MarkerSize',10)
    for i =1:length(Traject.GEO)
        plot(Traject.GEO{i}.theta,Traject.GEO{i}.TSR,'>-','MarkerSize',10)
    end
    subplot(1,2,2)
    contour(IdDat.theta2D,IdDat.lambda2D,IdDat.C_T,0.1:0.05:2.0,'ShowText','on')
    xlabel('theta')
    ylabel('lambda')
    xlim([restriction.theta(1) restriction.theta(2)])
    ylim([[restriction.TSR(1) restriction.TSR(2)]])
    title('Ct')
    hold on
    plot(Traject.TSR.theta,Traject.TSR.TSR,'O-','MarkerSize',10)
    hold on
    plot(Traject.theta.theta,Traject.theta.TSR,'X-','MarkerSize',10)
    hold on
    plot(Traject.minCT.theta,Traject.minCT.TSR,'^-','MarkerSize',10)
    for i =1:length(Traject.GEO)
        plot(Traject.GEO{i}.theta,Traject.GEO{i}.TSR,'>-','MarkerSize',10)
    end      
end

if flag.plot_steady_states_all_rates   
    plot_case = options.flag.steady_states_all_rates_case;
    
    figure
    subplot(2,2,1)
    for i = 1:size(SS,2)
        plot(options.V,SS{plot_case,i}.P_out,'linewidth',2),grid on,xlabel('WSP'),ylabel('Power')
        hold on
    end
    legend(arrayfun(@(a)num2str(a),[turbine.Pnom_el/10^6 P_rat_el_vec/10^6],'uni',0))
    title (['Controller ' options.control_names(plot_case)])
    
    subplot(2,2,2)
    for i = 1:size(SS,2)
        plot(options.V,SS{plot_case,i}.Cp_out,'linewidth',2),grid on,xlabel('WSP'),ylabel('Cp')
        hold on
    end
%     legend(arrayfun(@(a)num2str(a),[turbine.Pnom_el/10^6 P_rat_el_vec/10^6],'uni',0))
    
    subplot(2,2,3)
    for i=1:size(SS,2)
        plot(options.V,SS{plot_case,i}.P_out./SS{plot_case,1}.P_out,'linewidth',2),grid on,xlabel('WSP'),ylabel('P/Pnom')
        hold on
    end
%     legend(arrayfun(@(a)num2str(a),[turbine.Pnom_el/10^6 P_rat_el_vec/10^6],'uni',0))
    
    subplot(2,2,4)
    for i = 1:size(SS,2)
        plot(options.V,SS{plot_case,i}.Ct_out,'linewidth',2),grid on,xlabel('WSP'),ylabel('Ct')
        hold on
    end
%     legend(arrayfun(@(a)num2str(a),[turbine.Pnom_el/10^6 P_rat_el_vec/10^6],'uni',0))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    subplot(2,2,1)
    for i = 1:size(SS,2)
        plot(turbine.GB_ratio*SS{plot_case,i}.w_out,SS{plot_case,i}.Tq_out,'linewidth',2),grid on,xlabel('Gen Sp'),ylabel('Gen Tq')
        hold on
    end
%     legend(arrayfun(@(a)num2str(a),[turbine.Pnom_el/10^6 P_rat_el_vec/10^6],'uni',0))
    title (['Controller ' options.control_names(plot_case)])
    
    subplot(2,2,2)
    for i = 1:size(SS,2)
        plot(options.V,SS{plot_case,i}.theta_out,'linewidth',2),grid on,xlabel('WSP'),ylabel('Pitch')
        hold on
    end
%     legend(arrayfun(@(a)num2str(a),[turbine.Pnom_el/10^6 P_rat_el_vec/10^6],'uni',0))
    
    subplot(2,2,3)
    for i = 1:size(SS,2)
        plot(options.V,SS{plot_case,i}.TSR_out,'linewidth',2),grid on,xlabel('WSP'),ylabel('TSR')
        hold on
    end
%     legend(arrayfun(@(a)num2str(a),[turbine.Pnom_el/10^6 P_rat_el_vec/10^6],'uni',0))
    
    subplot(2,2,4)
    for i = 1:size(SS,2)
        plot(options.V,SS{plot_case,i}.w_out,'linewidth',2),grid on,xlabel('WSP'),ylabel('Rot speed')
        hold on
    end
    legend(arrayfun(@(a)num2str(a),[turbine.Pnom_el/10^6 P_rat_el_vec/10^6],'uni',0))
    
end

if flag.plot_steady_states_compare
    plot_case = options.steady_states_compare_case +1 ;
    
    figure   
    subplot(2,2,1)
    for i = 1:size(SS,1)
        plot(options.V,SS{i,plot_case}.P_out,'linewidth',2),grid on,xlabel('WSP'),ylabel('Power')
        hold on
    end
    title (['Pel = ' num2str(P_rat_el_vec(plot_case-1)/10^6) ' MW'])
    legend(options.control_names)
    
    subplot(2,2,2)
    for i = 1:size(SS,1)
        plot(options.V,SS{i,plot_case}.Cp_out,'linewidth',2),grid on,xlabel('WSP'),ylabel('Cp')
        hold on
    end
    
    subplot(2,2,3)
    for i=1:size(SS,1)
        plot(options.V,SS{i,plot_case}.P_out./SS{1,1}.P_out,'linewidth',2),grid on,xlabel('WSP'),ylabel('P/Pnom')
        hold on
    end
    
    subplot(2,2,4)
    for i = 1:size(SS,1)
        plot(options.V,SS{i,plot_case}.Ct_out,'linewidth',2),grid on,xlabel('WSP'),ylabel('Ct')
        hold on
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    subplot(2,2,1)
    for i = 1:size(SS,1)
        plot(turbine.GB_ratio*SS{i,plot_case}.w_out,SS{i,plot_case}.Tq_out,'linewidth',2),grid on,xlabel('Gen Sp'),ylabel('Gen Tq')
        hold on
    end
    title (['Pel = ' num2str(P_rat_el_vec(plot_case-1)/10^6) ' MW'])    
    
    subplot(2,2,2)
    for i = 1:size(SS,1)
        plot(options.V,SS{i,plot_case}.theta_out,'linewidth',2),grid on,xlabel('WSP'),ylabel('Pitch')
        hold on
    end
    
    
    subplot(2,2,3)
    for i = 1:size(SS,1)
        plot(options.V,SS{i,plot_case}.TSR_out,'linewidth',2),grid on,xlabel('WSP'),ylabel('TSR')
        hold on
    end
    
    subplot(2,2,4)
    for i = 1:size(SS,1)
        plot(options.V,SS{i,plot_case}.w_out,'linewidth',2),grid on,xlabel('WSP'),ylabel('Rot speed')
        hold on
    end
    legend(options.control_names)
    
end

%% Functions

function [theta,lambda,C_T,C_P,C_M] = InterpAndClean(theta,lambda,C_T,C_P,TSR_discr,theta_discr)

%Interpolation for higher resolution
[x,y] = ndgrid(lambda,theta);
FCp   = griddedInterpolant(x,y,C_P);
FCt   = griddedInterpolant(x,y,C_T);
[X,Y] = ndgrid(min(lambda):TSR_discr:max(lambda),min(theta):theta_discr:max(theta));
c_P_new = FCp(X,Y);
c_T_new = FCt(X,Y);
lambda  = X(:,1);
theta   = Y(1,:);
C_P = c_P_new;
C_T = c_T_new;
C_M = bsxfun(@rdivide,C_P,lambda);

% Clean up useless values:
C_P(C_P<0) = nan;
C_T(C_P<0) = nan;
C_M(C_P<0) = nan;
C_P(C_T<0) = nan;
C_T(C_T<0) = nan;
C_M(C_T<0) = nan;

end

function [lambda_vec,theta_vec,CP,CT,CM] = get_Cp_space(C_P,Cpmin,Cpmax,C_T,C_M,lambda,theta)

[lambda_ind,theta_ind] = find(C_P>=Cpmin & C_P<=Cpmax);
if ~isempty(lambda_ind)
    for i=1:length(lambda_ind)
        CP(i) = C_P(lambda_ind(i),theta_ind(i));
        CT(i) = C_T(lambda_ind(i),theta_ind(i)); %#ok<*AGROW>
        CM(i) = C_M(lambda_ind(i),theta_ind(i));
        lambda_vec(i) = lambda(lambda_ind(i));
        theta_vec(i) = theta(theta_ind(i));
    end
else
    lambda_vec = [];
    theta_vec = [];
    CP = [];
    CT = [];
    CM = [] ;
end
end

function out = GetCoefficient(val1,val2,Distr1,Distr2,DistrOut)

aa = abs(Distr1-val1);
[~,val1_ind] = min(aa);
bb=abs(Distr2-val2);
[~,val2_ind] = min(bb);
[sz1,~] = size(DistrOut);
if sz1==length(Distr1)
    out = DistrOut(val1_ind,val2_ind);
else
    out = DistrOut(val2_ind,val1_ind);
end
end

function out = GetCoefficientReverse(val1,val2,Distr1,Coeff2,DistrOut,tolerance,Evaluation)

Differ = abs(Distr1-val1);
[~,val1_ind] = min(Differ);
% Get subset of coefficient with the requested value
[sz1,~] = size(Coeff2);
if sz1==length(Distr1)
    Vals = Coeff2(val1_ind,:);
else
    Vals = Coeff2(:,val1_ind);
end
[bb,bb_ind] = sort(abs(Vals-val2));
cc = bb-bb(1);
[~,pos_cc]= find(cc<tolerance);
if numel(pos_cc)==1
    val2_ind = bb_ind(pos_cc);
else
    if isempty(Evaluation)
        [~, val2_ind] =  min(abs(Vals-val2));
    else
        %        [~, kk] = min(abs(DistrOut(bb_ind(pos_cc))-Evaluation(end))) ;
        jj = DistrOut(bb_ind(pos_cc));
        pp = find (jj<Evaluation(end));
        jj(pp) = 100; %#ok<*FNDSB>
        [~,jj_ind] = min(jj-Evaluation(end));
        val2_ind = bb_ind(jj_ind);
    end
end
out = DistrOut(val2_ind);
end

function SS = GetOperChar(TSR_val,theta_min_val,Pnom_val,w_min,r,theta,lambda,C_P,C_T,V,GB_ratio,losses,efficiency)

% ToDo:
% Add the transition region by checking rotational speed, generator torque and power. Put limits in them

Cp_val =  GetCoefficient(theta_min_val,TSR_val,theta,lambda,C_P) ;
Vnom_val  = (Pnom_val/(Cp_val*(0.5.*1.225).*(pi.*r.^2)) )^(1/3);
W_nom_val = convangvel(Vnom_val*TSR_val/r,'rad/s','rpm');

% Get steady states for the input values
for i = 1:length(V)
    
    %%% below rated %%%
%     SS.w_out(i) = convangvel(V(i)*TSR_val/r,'rad/s','rpm'); %#ok<*SAGROW>
    SS.w_out(i) = V(i)*TSR_val*(30/pi)/r; %#ok<*SAGROW>
    if SS.w_out(i) < w_min % min speed requirement
        SS.w_out(i) = w_min;
    end
%     SS.TSR_out(i)= convangvel(SS.w_out(i)*r/V(i),'rpm','rad/s');
    SS.TSR_out(i)= SS.w_out(i)*(pi/30)*r/V(i);
    SS.Cp_out(i) = GetCoefficient(theta_min_val,SS.TSR_out(i),theta,lambda,C_P) ;%Cp(theta_val,TSR_out(i));
    SS.Ct_out(i) = GetCoefficient(theta_min_val,SS.TSR_out(i),theta,lambda,C_T) ;%Ct(theta_val,TSR_out(i));
    SS.P_out(i)  = (efficiency)*SS.Cp_out(i)*(0.5.*1.225.*(V(i).^3).*(pi.*r.^2))-losses; % el power incl losses and efficiency
    SS.theta_out(i) = theta_min_val;
%     SS.Tq_out(i) = SS.Cp_out(i)*(0.5.*1.225.*(V(i).^3).*(pi.*r.^2))/convangvel(GB_ratio*SS.w_out(i),'rpm','rad/s');
    SS.Tq_out(i) = SS.Cp_out(i)*(0.5.*1.225.*(V(i).^3).*(pi.*r.^2))/(GB_ratio*SS.w_out(i)*(pi/30));
    
    %%% above rated %%%
    if  SS.Cp_out(i)*(0.5.*1.225.*(V(i).^3).*(pi.*r.^2)) >= Pnom_val
        SS.w_out(i)  = W_nom_val ;
        SS.P_out(i)  = Pnom_val;
%         SS.Tq_out(i) = SS.P_out(i)/convangvel(GB_ratio*SS.w_out(i),'rpm','rad/s');
        SS.Tq_out(i) = SS.P_out(i)/(GB_ratio*SS.w_out(i)*(pi/30));
        SS.Cp_out(i) = SS.P_out(i)/(0.5.*1.225.*(V(i).^3).*(pi.*r.^2));
        SS.P_out(i)  = (efficiency)*Pnom_val-losses; 
%         SS.TSR_out(i)= convangvel(SS.w_out(i)*r/V(i),'rpm','rad/s');
        SS.TSR_out(i) = SS.w_out(i)*(pi/30)*r/V(i);
        SS.theta_out(i) = GetCoefficientReverse(SS.TSR_out(i),SS.Cp_out(i),lambda,C_P,theta,1e-3,SS.theta_out) ;%theta(Cp,TSR);
        SS.Ct_out(i) = GetCoefficient(SS.theta_out(i),SS.TSR_out(i),theta,lambda,C_T) ;%Ct(theta_val,TSR_out(i));
        
    end
end
end

function feasib_restr = get_restricted_space( CT_out,CP_out,CM_out,lambda_vec,theta_vec,restriction)

TSR_restr   = find(lambda_vec<restriction.TSR(1) | lambda_vec>restriction.TSR(2));
theta_restr = find(theta_vec<restriction.theta(1) | theta_vec>restriction.theta(2));
ind_rem = unique([TSR_restr,theta_restr]);

feasib_restr.CT_out = CT_out;
feasib_restr.CP_out = CP_out;
feasib_restr.CM_out = CM_out;
feasib_restr.lambda_vec = lambda_vec;
feasib_restr.theta_vec = theta_vec;
% now remove the values that are out of the space
if ~isempty(ind_rem)
    feasib_restr.CT_out(ind_rem) = [];
    feasib_restr.CP_out(ind_rem) = [];
    feasib_restr.CM_out(ind_rem) = [];
    feasib_restr.lambda_vec(ind_rem) = [];
    feasib_restr.theta_vec(ind_rem) = [];
end

end

function pos_out = get_disamb_val (InpSpaceX,targetX,InpSpaceY,tolerance,EvaluationX,EvaluationY)
%choose the set point from multiple that are close based on the proximity
%to the previous point. It is implemented to avoid jumps in the points when
%a trajectory is followed. Better implementation can be thought of...

%sort the points according to procimity to target
[aa,pos_aa] = sort(abs(InpSpaceX-targetX));
% find the difference in distance between the closest point and the next closer points due to sorting the first is the closest
bb =aa-aa(1);
% extract the points within a tolerance
[~,pos_cc]= find(bb<tolerance); %be careful of the manual tolerance here!

if numel(pos_cc)==1
    pos_out = pos_aa(pos_cc); % only one point found in tolerance no test
else  % more than one point found: get the minimum euclidean distance
    Y(1) = EvaluationX(end);
    Y(2) = EvaluationY(end);
    for i =1:length(pos_cc)
        X(1) = InpSpaceX(pos_aa(pos_cc(i)));
        X(2) = InpSpaceY(pos_aa(pos_cc(i)));
        D(i) = norm(X-Y);
        %     [dd,dd_pos] = min(abs(feasib_restr.theta_vec(pos_aa(pos_cc))-EvaluationX(end)));
    end
    [~,pos_d] = min(D);
    pos_out = pos_aa(pos_d);
end
end










