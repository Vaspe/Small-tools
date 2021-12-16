clc
clearvars
close all

% Thoretical caclulation of lean angle speed, corner radius, friction coefficient and braking/cornenring. It is
% oversimplified but can be used as a first order aproxiamation
% For equations see 
% https://www.stevemunden.com/leanangle.html
% https://web.archive.org/web/20130105031937/http://www.howfastcanigo.com/howitworks.html
% https://en.wikipedia.org/wiki/Bicycle_and_motorcycle_dynamics
%
% Vasilis Pettas  19/10/2020

%% Inputs
r_cor = 13;  % radius of corner in meters  see this source https://nacto.org/publication/urban-street-design-guide/intersection-design-elements/corner-radii/
m_mot = 250; % mass of motorcycle in kg with rider and equipment
lean_ang = 30;   % lean angle in degrees
vel_mot_in = 96.6; % forward velocity in km/h
g = 9.81; % gravitational constant in SI
acc_mot = 1.1*(-9.81) ; % constant deccelerations m/s^2 1*9.81 = 1g (1-1.1g average is the maximum for best motorcycles today)
wdth_tire = 0.180 ; %m width of rear tires in meters
ht_tire = 0.1; % m height of tire
% CG: by positioning out body towards
% the inner part of the corner and into the tank we decrease it but also change the
% reduce the lean angle itself by changing the location of the CG relevant
% to y plane (x forward motion axis, z height axis, y side to side axis)
% here we dont take into consideration the body positon in the y axis
h_CG = 0.54; % m center of gravity
Whb = 1.386; % m wheelbase 
rake = 24;   % degrees caster angle

% Swipe data
lean_vec = 0:60; % vector to swipe for lean angles in degrees
speed_vec_in = 0:160; % vector to swipe for speeds in kmh 
r_vec =  0:0.5:75; %vector to swipe for corner radius in meters
breakdist_vec = 10:200;
acc_vec = [2:0.1:15]; %#ok<*NBRAK> %vector to swipe for mean negative accelerations  
steer_angle_vec = linspace(0,20,length(lean_vec)); %steering angle input


%% Calculations
speed_vec = speed_vec_in *1000/3600; %convert to m/s
vel_mot = vel_mot_in* 1000/3600; %convert to m/s
t_rear = wdth_tire/2; % we use the half of the width for the calculaiton
lean_vec_cor = lean_vec + asind(t_rear*sind(lean_vec)/(h_CG-t_rear)); % corrected vector of lean angles with tire widtht ang CG
lean_ang_cor = lean_ang + asind(t_rear*sind(lean_ang)/(h_CG-t_rear)); % corrected input lean angle with tire widtht ang CG

% Calculate velocity for a range of lean angles and a fixed radius
v_sweep = sqrt(tand(lean_vec).*g.*r_cor).*(3600/1000); % in km/h
% Calculate velocity for a range of radii and a fixed lean angle
v_sweep2 = sqrt(tand(lean_ang).*g.*r_vec).*(3600/1000); % in km/h

% Calculate lean angle for a range of radii and fixed velocity
lean_sweep_ideal = atand((vel_mot^2)./(r_vec.*g)); % lean angle in degrees
lean_sweep_cor = lean_sweep_ideal + asind(t_rear*sind(lean_sweep_ideal)/(h_CG-t_rear)); % correcting the lean angle vector accorfing to with of tire and CG height)
% Calculate lean angle for a range of velocities and a fixed radius
lean_sweep2_ideal = atand((speed_vec.^2)./(r_cor.*g)); % lean angle in degrees
lean_sweep2_cor = lean_sweep2_ideal + asind(t_rear*sind(lean_sweep2_ideal)/(h_CG-t_rear)); % correcting the lean angle vector accorfing to with of tire and CG height)


% Calculate radius of corner for a range of lean angles and a fixed velocity
r_sweep = (vel_mot.^2) ./ (g.*tand(lean_vec)) ; % radius in m
% Calculate radius of corner for a range of velocities and a fixed lean angle
r_sweep2 = (speed_vec.^2) ./ (g.*tand(lean_ang)) ; % radius in m

% calculate the turn radius decrease according to bikes geomerty, steering angle and lean angle
r_decr = Whb .* cosd(lean_vec) ./(cosd(steer_angle_vec).*cosd(rake));

%% Braking calculation to stop 
% Calculate break distance for a range of speeds with fixed ideal constant deccelerations
beakdist_sweep = (speed_vec.^2)./(2*(-acc_mot));

% Calculate break distance for a range of deccelerations with fixed speeds
beakdist_sweep2 = (vel_mot.^2)./(2*(acc_vec));


%% Plotting

% Ideal lean angle
figure
subplot(3,2,1)
plot(lean_vec,v_sweep,'linewidth',2)
xlabel('ideal lean angle in deg','FontSize', 10)
ylabel('km/h','FontSize', 12)
title(['For a corner with radius = ' num2str(r_cor) ' m'])
grid on

subplot(3,2,2)
plot(r_vec,v_sweep2,'linewidth',2)
xlabel('Corner radius in m','FontSize', 10)
ylabel('km/h','FontSize', 10)
title(['For ideal lean angle = ' num2str(lean_ang) ' deg'])
grid on


subplot(3,2,3)
plot(r_vec,lean_sweep_ideal,'linewidth',2)
xlabel('Corner radius in m','FontSize', 10)
ylabel('ideal lean angle in deg','FontSize', 10)
title(['For a speed of = ' num2str(vel_mot_in) ' km/h'])
grid on

subplot(3,2,4)
plot(speed_vec_in,lean_sweep2_ideal,'linewidth',2)
xlabel('Speed in km/h','FontSize', 10)
ylabel('ideal lean angle in deg','FontSize', 10)
title(['For a corner radius of = ' num2str(r_cor) ' m'])
grid on


subplot(3,2,5)
plot(lean_vec,r_sweep,'linewidth',2)
xlabel('ideal lean angle in deg','FontSize', 10)
ylabel('Corner radius in m','FontSize', 10)
title(['For a speed of = ' num2str(vel_mot_in) ' km/h'])
grid on

subplot(3,2,6)
plot(speed_vec_in,r_sweep2,'linewidth',2)
xlabel('Speed in km/h','FontSize', 10)
ylabel('Corner radius in m','FontSize', 10)
title(['For ideal lean angle = ' num2str(lean_ang) ' deg'])
grid on


% Corrected lean angle
figure
subplot(3,2,1)
plot(lean_vec_cor,v_sweep,'linewidth',2)
xlabel('Corrected lean angle in deg','FontSize', 10)
ylabel('km/h','FontSize', 12)
title(['For a corner with radius = ' num2str(r_cor) ' m'])
grid on

subplot(3,2,2)
plot(r_vec,v_sweep2,'linewidth',2)
xlabel('Corner radius in m','FontSize', 10)
ylabel('km/h','FontSize', 10)
title(['For corrected lean angle = ' num2str(lean_ang_cor) ' deg'])
grid on


subplot(3,2,3)
plot(r_vec,lean_sweep_cor,'linewidth',2)
xlabel('Corner radius in m','FontSize', 10)
ylabel('Corrected lean angle in deg','FontSize', 10)
title(['For a speed of = ' num2str(vel_mot_in) ' km/h'])
grid on

subplot(3,2,4)
plot(speed_vec_in,lean_sweep2_cor,'linewidth',2)
xlabel('Speed in km/h','FontSize', 10)
ylabel('Corrected lean angle in deg','FontSize', 10)
title(['For a corner radius of = ' num2str(r_cor) ' m'])
grid on


subplot(3,2,5)
plot(lean_vec_cor,r_sweep,'linewidth',2)
xlabel('Corrected lean angle in deg','FontSize', 10)
ylabel('Corner radius in m','FontSize', 10)
title(['For a speed of = ' num2str(vel_mot_in) ' km/h'])
grid on

subplot(3,2,6)
plot(speed_vec_in,r_sweep2,'linewidth',2)
xlabel('Speed in km/h','FontSize', 10)
ylabel('Corner radius in m','FontSize', 10)
title(['For corrected lean angle = ' num2str(lean_ang_cor) ' deg'])
grid on


% decceleration
figure
subplot(2,2,1)
plot(speed_vec_in,beakdist_sweep,'linewidth',2)
xlabel('Speed in km/h','FontSize', 10)
ylabel('Break to stop distance in m','FontSize', 10)
title(['For a negative acceleration = ' num2str(acc_mot) ' m/s2'])
grid on

subplot(2,2,2)
plot(acc_vec,beakdist_sweep2,'linewidth',2)
xlabel('Negative acceleration in m/s2','FontSize', 10)
ylabel('Break to stop distance in m','FontSize', 10)
title(['For initial speed = ' num2str(vel_mot_in) ' km/h'])
grid on

subplot(2,2,3)
plot(0:70,tand(0:70),'linewidth',2)
xlabel('Ideal lean angle deg')
ylabel('Coefficient of friction required')
title('Ideal friction calulation for cornering')
grid on 

subplot(2,2,4)
plot(acc_vec/g,beakdist_sweep2,'linewidth',2)
xlabel('Friction coefficient [-] max ideally 1.5 in track','FontSize', 10) % bmw s1000rr in test had 1.15 with road tires!!!!
ylabel('Break to stop distance in m','FontSize', 10)
title(['For initial speed = ' num2str(vel_mot_in) ' km/h'])
grid on

% figure
% plot3(lean_vec,steer_angle_vec,r_decr,'linewidth',2)
% xlabel('ideal lean angle in deg','FontSize', 10) % bmw s1000rr in test had 1.15 with road tires!!!!
% ylabel('steer angle in deg','FontSize', 10)
% zlabel('Decrease of radius turn','FontSize', 10)
% % title('Decrease of radius turn')
% grid on
% ancillary
% figure,plot (0:100,[0:100]*1000/3600,'linewidth',2),xlabel('km/h'),ylabel('m/s'),title('km/h to m/s'),grid on 
