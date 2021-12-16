% Script to calculate scanning points for the nacelle mounted lidar measuring patterns.
% Inputing a cartesian trajectory in a reference plane and outputing the
% trajectory in different planes along with the cartesian and spherical
% cartesians. Outputs can be written in txt file for uploading on a lidar
% device both for step&stare and continuous modes. Also calculates scanning
% times for a rough prediction of sampling frequency 
%
% Lidar accepts spherical coordinates. Azimuth is 0 looking in front and 90
% looking to the right. Elevation is 0 looking at Lidar height to the rotor plane and 90
% looking straight up. Origin is the lidar
%
% In cartesian x is perpendicular to the rotor, z is height, y is parallel to
% the diameter. Positive x means in front of the rotor, y left and z up, as seen from Lidar to blades
% Origin is the the hub center
%
% (D senvion 126m, Lidar-hub x offset ca. 4m, z offset 2m. -15deg is the 
%  minimum elevation for the Lidar)
%
% Dependencies: comet3_mod function a=for plotting
% ** mostly compatible with Octave (except in plotting)
%
%VPE 29.11.2018 v1

clc,clear all, close all %#ok<*CLALL>
run('..\AddingWitlisPaths.m')

%% Inputs in cartesian CS

planes = [100  250  500  750 1000]; % measurement planes x [m]
ref_pl = 250; % reference plane on which the patterns are calculated LOS dictates the position in the other planes [m] No offset considered here!
off_x  = 0 ;  % Lidar offset from hub center x [m] 4?
off_z  = 0 ;  % Lidar offset from hub center z [m] 2?
off_y  = 0 ;  % Lidar offset from hub center y [m]
flag_downstream = 1; %if this flag is used thetrajectory will switch and look at the downstream mirroring the upstream definition

orig_patern = {       % y = line 1, z = line 2, expressed on ref_pl / the last extra point is the stare point, mandatory for Streamline XR v14
%            [-55,-36.6,-18.3,0,18.3,36.6,55,0 ; 0,0,0,0,0,0,0,0];   % 7 point horizontal line +1 staring  
%            [0,0,0,0,0,0,0,0 ; -55,-36.6,-18.3,0,18.3,36.6,55,0];   % 7 point vertical line           
%            [0,47.6,47.6,0,-47.6,-47.6,0,0 ; 55,27.5,-27.5,-55,-27.5,27.5,0,0 ]; % 7 point circle
%            [-55,-55,-55,0,0,0,55,55,55,0 ; 55,0,-55,-55,0,55,55,0,-55,0]; % 9 point rectangular grid
%            [55,27.5,0,-27.5,-55,0,0,0,0,0 ; 0,0,0,0,0,55,27,-27,-55,0 ];% 9 point cross    
%            [46.7,23.3,0,-23.3,-46.7,46.7,22.9,-22.9,-46.7,0 ; -46.7,-23.3,0,23.3,46.7,46.7,22.9,-22.9,-46.7,0];% 9 point X
%            [-55,55,0 ; 0,0,0];   % 2 point horizontal line +1 staring  
%            [0,0,0 ; -55,55,0];   % 2 point vertical line     
%            [-55,-55,0,0,55,55,0 ; 55,-55,-55,55,55,-55,0]; % 6 point rectangular grid     
%            [55,-55,0,0,0 ; 0,0,55,-55,0 ];% 4 point cross   
%            [46.7,-46.7,46.7,-46.7,0 ; -46.7,46.7,46.7,-46.7,0];%4 point X   

           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51	-85.51];   % 41 point horizontal linefor Alpha Ventus campaign -20 elev          
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39	-81.39];   % 41 point horizontal linefor Alpha Ventus campaign -19 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25	-77.25];   % 41 point horizontal linefor Alpha Ventus campaign -18 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09	-73.09];   % 41 point horizontal linefor Alpha Ventus campaign -17 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91	-68.91];   % 41 point horizontal linefor Alpha Ventus campaign --16 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7	-64.7];   % 41 point horizontal linefor Alpha Ventus campaign -15 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48	-60.48];   % 41 point horizontal linefor Alpha Ventus campaign -14 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24	-56.24];   % 41 point horizontal linefor Alpha Ventus campaign -13 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98	-51.98];   % 41 point horizontal linefor Alpha Ventus campaign -12 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7	-47.7];   % 41 point horizontal linefor Alpha Ventus campaign -11 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41	-43.41];   % 41 point horizontal linefor Alpha Ventus campaign -10 elev           
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11	-39.11];  % 41 point horizontal linefor Alpha Ventus campaign -9 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79	-34.79];   % 41 point horizontal linefor Alpha Ventus campaign -8 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47	-30.47];  % 41 point horizontal linefor Alpha Ventus campaign -7 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13	-26.13]; % 41 point horizontal linefor Alpha Ventus campaign -6 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79	-21.79]; % 41 point horizontal linefor Alpha Ventus campaign -5 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44	-17.44]; % 41 point horizontal linefor Alpha Ventus campaign -4 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08	-13.08]; % 41 point horizontal linefor Alpha Ventus campaign -3 elev        
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72	-8.72];  % 41 point horizontal linefor Alpha Ventus campaign -2 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36	-4.36];  % 41 point horizontal linefor Alpha Ventus campaign -1 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];                                                                                                                                                                      % 41 point horizontal linefor Alpha Ventus campaign 0 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36	4.36];   % 41 point horizontal linefor Alpha Ventus campaign 1 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72	8.72];   % 41 point horizontal linefor Alpha Ventus campaign 2 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08	13.08];  % 41 point horizontal linefor Alpha Ventus campaign 3 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44	17.44];  % 41 point horizontal linefor Alpha Ventus campaign 4 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79	21.79];  % 41 point horizontal linefor Alpha Ventus campaign 5 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13	26.13];  % 41 point horizontal linefor Alpha Ventus campaign 6 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47	30.47];  % 41 point horizontal linefor Alpha Ventus campaign 7 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79	34.79];  % 41 point horizontal linefor Alpha Ventus campaign 8 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11	39.11];  % 41 point horizontal linefor Alpha Ventus campaign 9 elev
           [-160.7,-153.9,-146.95,-139.8,-132.5,-125.0,-117.37,-109.6,-101.7,-93.65,-85.5,-77.25,-68.9,-60.48,-51.98,-43.4,-34.79,-26.13,-17.44,-8.72,0,8.7,17.44,26.13,34.79,43.4,51.98,60.48,68.9,77.25,85.5,93.65,101.68,109.6,117.37,125.0,132.48,139.8,146.95,153.92,160.696902421635 ;43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41	43.41];  % 41 point horizontal linefor Alpha Ventus campaign 10 elev
           };
                      
pat_names= {                     % name_scan mode(ss/cs)_#pt_#range gate
%             'HorLine_7pt_31d5';   %1
%             'VerLine_7pt_31d5';   %2
%             'Circle_7pt_31d5'     %3
%             'GridRect_9pt_31d5'   %4
%             'Cross_9pt_31d5'      %5
%             'X_9pt_31d5'          %6
%             'HorLine_2pt_31d5';   %7
%             'VerLine_2pt_31d5';   %8
%             'GridRect_6pt_31d5'   %9
%             'Cross_4pt_31d5'      %10
%             'X_4pt_31d5'          %11
            
            
            'S_HorLine_41pt_30_m20elev';    %12  % Alpha Ventus TV scan -20 % name_scan mode(ss/cs)_#pt_#range gate_elevation
            'S_HorLine_41pt_30_m19elev';    %13  % Alpha Ventus TV scan -19
            'S_HorLine_41pt_30_m18elev';    %14  % Alpha Ventus TV scan -18
            'S_HorLine_41pt_30_m17elev';    %15  % Alpha Ventus TV scan -17
            'S_HorLine_41pt_30_m16elev';    %16  % Alpha Ventus TV scan -16
            'S_HorLine_41pt_30_m15elev';    %17  % Alpha Ventus TV scan -15
            'S_HorLine_41pt_30_m14elev';    %18  % Alpha Ventus TV scan -14
            'S_HorLine_41pt_30_m13elev';    %19  % Alpha Ventus TV scan -13
            'S_HorLine_41pt_30_m12elev';    %20  % Alpha Ventus TV scan -12
            'S_HorLine_41pt_30_m11elev';    %21  % Alpha Ventus TV scan -11
            'S_HorLine_41pt_30_m10elev';    %22  % Alpha Ventus TV scan -10
            'S_HorLine_41pt_30_m9elev';    %23  % Alpha Ventus TV scan -9
            'S_HorLine_41pt_30_m8elev';    %24  % Alpha Ventus TV scan -8
            'S_HorLine_41pt_30_m7elev';    %25  % Alpha Ventus TV scan -7
            'S_HorLine_41pt_30_m6elev';    %26  % Alpha Ventus TV scan -6
            'S_HorLine_41pt_30_m5elev';    %27  % Alpha Ventus TV scan -5        
            'S_HorLine_41pt_30_m4elev';    %28  % Alpha Ventus TV scan -4           
            'S_HorLine_41pt_30_m3elev';    %29  % Alpha Ventus TV scan -3            
            'S_HorLine_41pt_30_m2elev';    %30  % Alpha Ventus TV scan -2           
            'S_HorLine_41pt_30_m1elev';    %31  % Alpha Ventus TV scan -1            
            'S_HorLine_41pt_30_0elev';    %32  % Alpha Ventus TV scan 0            
            'S_HorLine_41pt_30_1elev';    %33 % Alpha Ventus TV scan 1            
            'S_HorLine_41pt_30_2elev';    %34  % Alpha Ventus TV scan 2            
            'S_HorLine_41pt_30_3elev';    %35  % Alpha Ventus TV scan 3                        
            'S_HorLine_41pt_30_4elev';    %36  % Alpha Ventus TV scan 4                        
            'S_HorLine_41pt_30_5elev';    %37  % Alpha Ventus TV scan 5                        
            'S_HorLine_41pt_30_6elev';    %38  % Alpha Ventus TV scan 6                        
            'S_HorLine_41pt_30_7elev';    %39  % Alpha Ventus TV scan 7                        
            'S_HorLine_41pt_30_8elev';    %40  % Alpha Ventus TV scan 8                        
            'S_HorLine_41pt_30_9elev';    %41  % Alpha Ventus TV scan 9                        
            'S_HorLine_41pt_30_10elev';    %42  % Alpha Ventus TV scan 10                                    
            };
            
Azim_steps  = 50e4;  % number of steps for a full 360 in azimuth
Elev_steps  = 25e4;  % number of steps for a full10 360 in elevation
max_mot_sp  = 5e4;   % maximum speed in step/second unit --> 5e4 steps/s --> 10s for a 360 in azi and 5s in elevation
max_acc     = 5e4;   % maximum speed in step/second^2 unit --> 5e4 steps/s^2 --> 1s needed to reach a speed of 5e4 steps/s
wait_time   = 0;     % [ms]time to wait after each point in continuous scanning mode (if 0 it will still reach 0 velocity and then accelerate again)

ss_mot_sp_azi  = 5e4; % speed in step/second unit is step&stare mode azimuth motor/ 5e4 = 36 deg/s
ss_mot_acc_azi = 3e4; % acceleration in step/second^2 unit is step&stare mode azimuth motor /3e4 = 21.6 deg/s^2
ss_mot_sp_ele  = 5e4; % speed in step/second unit is step&stare mode elevation motor /5e4 = 72 deg/s
ss_mot_acc_ele = 5e4; % acceleration in step/second^2 unit is step&stare mode azimuth motor /5e4 = 72 deg/s^2
CSM_sp         = 36;  % [deg/s] constant equal or smaller than the azimuth max speed: User defined
CSM_acc        = 36;  % [deg/s^2] constant equal or smaller than the azimuth max acc: User defined
CSM_wait       = 0;   % wait in points of CSM trajectory [ms]

max_pulse_rate = 10e3;   % pulses sent from laser per second. Increasing means increasing volume of data and accuracy in long range  
req_pulse_rate = 10e3;   % requestesd number of pulses <= max
comet_plot_pattern = 1;  % choose which pattern to be plotted with

repeat_traj = 25; % choose how many times to repeat the trajectory in the file 1 means it is only the trajectory described here

% save_folder = 'D:\Dropbox\SWE\Tasks\OWP controls\Lidar measurements\Scaning_patterns\'; N:\SWE\70_Messdaten\91_NeueKampagnenDokumentation\2019_AlphaVentus_OWPControlParkCast_StreamLineXR
save_folder = 'N:\SWE\70_Messdaten\91_NeueKampagnenDokumentation\2019_AlphaVentus_OWPControlParkCast_StreamLineXR\Initial_Trajectories\';
save_flag   = 1; 
plot_flag   = 1;
plot_flag_comet = 0;


%% Introducing offsets due to lidar misalignment with hub          
for i=1:length(orig_patern)
    patern{i}(2,:) = orig_patern{i}(2,:) - off_z ;    
    patern{i}(1,:) = orig_patern{i}(1,:) - off_y ;        
end
planes = planes + off_x;  % x dirextion offsets only
ref_pl = ref_pl + off_x;  % x dirextion offsets only         
           
%% Calculate in spherical CS 

% These are effectively the step and stare trajectory points
for i=1:length(patern)

    % define x,y,z for each pattern
    iY=patern{i}(1,:);
    iZ=patern{i}(2,:);
    iX=ref_pl*ones(1,length(iY));   
    %loop over every point and convert to spherical
    for j = 1 : length(iX)
        r(i,j)     = sqrt(iX(j).^2+iY(j).^2+iZ(j).^2); %radius not usefull here since LOS is the radius in multiple points
        phi(i,j)   = -atand (iY(j)./iX(j)) ;  %azimuth
        theta(i,j) = asind (iZ(j)./r(i,j)) ;  %#ok<*SAGROW> %elevation
    end
end

%% Calculate ideal time of trajectory for step&stare mode (including the stare point at the end)
scan_time = req_pulse_rate/max_pulse_rate; % time spent in each point of the ss trajectory

%reaching max speed and stopping, initial condtion: stand still
azi_max_v_t = ss_mot_sp_azi/ss_mot_acc_azi; %s
azi_max_v_x = 0.5*ss_mot_acc_azi*azi_max_v_t^2; %steps
ele_max_v_t = ss_mot_sp_ele/ss_mot_acc_ele; %s
ele_max_v_x = 0.5*ss_mot_acc_ele*ele_max_v_t^2; %steps  
for i=1:length(patern)
    iphi   = phi(i,1:length(patern{i}));
    itheta = theta(i,1:length(patern{i}));
    phi_steps   = [diff(iphi) iphi(1)-iphi(end)]; % steps in degrees required for each trajectory point
    theta_steps = [diff(itheta) itheta(1)-itheta(end)] ;
    phi_req     = phi_steps*Azim_steps/360;     % convert to motor steps required
    theta_req   = theta_steps*Elev_steps/360; 
    % loop over each step to calculate the time (from speed 0 to speed 0)
    for k=1:length(phi_steps) 
        if abs(phi_req(k)) <= 2*azi_max_v_x
           phi_t = 2* sqrt(2*(abs(phi_req(k))/2) / ss_mot_acc_azi) ; %acceleration and deacceleration equal for half the distance
        else 
           phi_t = 2*azi_max_v_t + (abs(phi_req(k))-2*azi_max_v_x)/ss_mot_sp_azi ;%2*acc + steady max sp    
        end
        if abs(theta_req(k)) <= 2*ele_max_v_x
           theta_t = 2* sqrt(2*(abs(theta_req(k))/2) / ss_mot_acc_ele) ; %acceleration and deacceleration equal for half the distance
        else 
           theta_t = 2*ele_max_v_t + (abs(theta_req(k))-2*ele_max_v_x)/ss_mot_sp_ele ;%2*acc + steady max sp    
        end
        t_step(i,k) = max(phi_t,theta_t);        
    end
    traj_ss_time(i) = sum(t_step(i,:))+ length(patern{i})*scan_time;  
    
    clear iphi itheta phi_steps theta_steps phi_req theta_req
end

%% Calculate continuous scanning mode inputs and ideal time

% In continuous scanning mode the Lidar reverts the azimuth CS making it
% counterclockwise positive as seen from z+ and elevation is positive
% towards the ground
phi_CSM = -phi;
theta_CSM = -theta;
% CSM_max_v_azi =CSM_sp*Azim_steps/360 ; %maximum speed in motor steps/s
% CSM_max_v_ele =CSM_sp*Elev_steps/360 ; %maximum speed in motor steps/s
% for patterns consisting only of straight lines the conversion to CSM is
% straight forward, for patterns including curves...

for i=1:length(patern)
    iphi   = phi_CSM(i,1:length(patern{i}));
    itheta = theta_CSM(i,1:length(patern{i}));
    phi_steps   = abs([diff([0,iphi]) iphi(1)-iphi(end)]); % steps in degrees required for each trajectory point
    theta_steps = abs([diff([0,itheta]) itheta(1)-itheta(end)]) ;
    phi_req     = iphi*Azim_steps/360;     % convert to motor steps required
    theta_req   = itheta*Elev_steps/360; 
    
    % since azimuth motor is slower and in order to have a smooth pattern the
    % magnitude of 2D speed vector should be constant and equal or smaller than the azimuth max speed
    % *assuming instataneous 0 speed at each scanning point (?)
    for k=1:length(phi_steps)
         if phi_steps(k)~=0
            sf_ang = atand(phi_steps(k)/theta_steps(k)); %scaling factor angle
         elseif phi_steps(k)==0 && phi_steps(k)==0
            sf_ang = 0;    
         else
            sf_ang = 90;
         end
         sf_azi = sind(sf_ang);   %scaling factor based on the length of the step that each motor has to perform
         sf_ele = cosd(sf_ang);
         v_azi_req = sf_azi*CSM_sp*Azim_steps/360; %required speed in motor steps/s
         v_ele_req = sf_ele*CSM_sp*Elev_steps/360;
         acc_azi_req = sf_azi*CSM_acc*Azim_steps/360; % required acceleration in motor steps/s^2
         acc_ele_req = sf_ele*CSM_acc*Elev_steps/360;
        
        % Calculating duaration of the trajectory (including the staring at the end) 
        % required acceleration intervals based on step's speed and acceleration 
        if acc_azi_req == 0 && v_azi_req ==0
            azi_max_v_t = 0;
        else
            azi_max_v_t = v_azi_req/acc_azi_req; %s
        end
        azi_max_v_x = 0.5*acc_azi_req*azi_max_v_t^2; %steps
        if acc_ele_req == 0 && v_ele_req ==0
            ele_max_v_t = 0;
        else
            ele_max_v_t = v_ele_req/acc_ele_req; %s
        end
        ele_max_v_x = 0.5*acc_ele_req*ele_max_v_t^2; %steps         
        if acc_azi_req==0 && v_azi_req==0
            phi_t = 0;
        elseif abs(phi_steps(k)*Azim_steps/360) <= 2*azi_max_v_x
           phi_t = 2* sqrt(2*(abs(phi_steps(k)*Azim_steps/360)/2) / acc_azi_req) ; %acceleration and deacceleration equal for half the distance
        else 
           phi_t = 2*azi_max_v_t + (abs(phi_steps(k)*Azim_steps/360)-2*azi_max_v_x)/abs(v_azi_req) ;%2*acc + steady max sp    
        end
        if acc_ele_req==0 && v_ele_req==0
            theta_t = 0;
        elseif abs(theta_steps(k)*Elev_steps/360) <= 2*ele_max_v_x
           theta_t = 2* sqrt(2*(abs(theta_steps(k)*Elev_steps/360)/2) / acc_ele_req) ; %acceleration and deacceleration equal for half the distance
        else 
           theta_t = 2*ele_max_v_t + (abs(theta_steps(k)*Elev_steps/360)-2*ele_max_v_x)/abs(v_ele_req) ;%2*acc + steady max sp    
        end
        t_step_CSM(i,k) = max(phi_t,theta_t); 
        
        V_CSM_azi(i,k)  = v_azi_req ;
        A_CSM_azi(i,k)  = acc_azi_req;
        V_CSM_ele(i,k)  = v_ele_req ;
        A_CSM_ele(i,k)  = acc_ele_req;        
        if k<=length(phi_req)
            Phi_CSM(i,k)   = (phi_req(k));
            Theta_CSM(i,k) = (theta_req(k)) ;  
        else
            Phi_CSM(i,k)   = 0;
            Theta_CSM(i,k) = 0;  
        end
    end    
    traj_CSM_time(i) = sum(t_step_CSM(i,:));  
    
    clear iphi itheta phi_steps theta_steps phi_req theta_req
end

% convert to correct ouput format for lidar CSM input file. The acc and
% speed values refer to the movement to reach this point. First set of A,S is ignored and defined by the internal software 
V_CSM_azi_exp = abs(round(V_CSM_azi/10));  % defined in 10steps/s
A_CSM_azi_exp = abs(round(A_CSM_azi/1000));  % defined in 1000steps/s
V_CSM_ele_exp = abs(round(V_CSM_ele/10));  
A_CSM_ele_exp = abs(round(A_CSM_ele/1000));  
Phi_CSM_exp   = [round(Phi_CSM),round(Phi_CSM)];
Theta_CSM_exp = [round(Theta_CSM),round(Theta_CSM)];  

%% Switch to scanning downstream if requested 

if flag_downstream==1  %changing only azimuth! elevation is the same
    
    phi   = phi+180;
    ind_el = find(phi>360); % remove nvalues higher than 360
    phi(ind_el) = phi(ind_el)-360;
    ind_mn = find(phi<0); % remove negative values
    phi(ind_mn) = phi(ind_mn)+360; 

end

%% calculate projection in planes - azimuth and elevation is the same have to find where it cuts the requested planes

for i = 1:length(patern) %loop over all paterns

  for k = 1:length(planes)   % loop over each scanning plane
      if flag_downstream==1
          kX = -planes(k);
      else
          kX = planes(k);  
      end
      for jj = 1:length(patern{i}) % loop over each point of the pattern
          rr = ( kX/(cosd(theta(i,jj))*cosd(phi(i,jj)) ) ) ;
          yy = rr*cosd(theta(i,jj))*sind(-phi(i,jj));
          zz = rr*sind(theta(i,jj));
          All_plane_int{k}(jj,:) = [kX,yy,zz];  % cartesians in x,y,z format
          All_plane_spher_int{k}(jj,:) = [rr,phi(i,jj),theta(i,jj)];        % sphericals in r,phi,theta format 
                           % LOS distance (r) is needed for calculating useful range gates in post processing
      end  
  end
  All_patern_cart{i}  =  All_plane_int;
  All_patern_spher{i} =  All_plane_spher_int;
  clear All_plane_int
end 

%% Plotting

if plot_flag==1
    % plot patterns initial
    figure,
    for i=1:length(patern)
        if flag_downstream==1
            plot3(-ref_pl*ones(1,length(patern{i}(1,:))),patern{i}(1,:),patern{i}(2,:),'--X','Linewidth',2)
        else
            plot3(ref_pl*ones(1,length(patern{i}(1,:))),patern{i}(1,:),patern{i}(2,:),'--X','Linewidth',2)
        end
        hold on
    end
    plot3(0,0,0,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],...
        'Marker','o','LineStyle','none');
    hold off
    grid on

%     Pattern in different planes 
    for ii=1:length(patern)
        figure,
        for i=1:length(planes)
            
            plot3(All_patern_cart{1,ii}{1, i}(:,1),All_patern_cart{1, ii}{1, i}(:,2),All_patern_cart{1, ii}{1, i}(:,3),'--X','Linewidth',2)
%              plot3(All_patern{1,ii}{1, i}(:,1),All_patern{1, ii}{1, i}(:,2),All_patern{1, ii}{1, i}(:,3))

            hold on
        end
        plot3(0,0,0,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],...
        'Marker','o','LineStyle','none');
        hold off
        grid on
        set(gca,'YGrid','on', 'FontSize', 14,'YMinorGrid','on','XMinorGrid','on','ZMinorGrid','on')
        title(pat_names{ii}, 'Interpreter', 'none')
    end

    % Comet plot
    if plot_flag_comet==1
        figure,
        for i=4:4
            for kkk=1:2
                comet3_mod(All_patern_cart{1,comet_plot_pattern}{1, i}(:,1),All_patern_cart{1, comet_plot_pattern}{1, i}(:,2),All_patern_cart{1, comet_plot_pattern}{1, i}(:,3),0.01,0.3)
                %          plot3(All_patern{1,ii}{1, i}(:,1),All_patern{1, ii}{1, i}(:,2),All_patern{1, ii}{1, i}(:,3))
                
                hold on
            end
        end
        hold off
        grid on
        set(gca,'YGrid','on', 'FontSize', 14,'YMinorGrid','on','XMinorGrid','on','ZMinorGrid','on')
        title(pat_names{ii}, 'Interpreter', 'none')
    end
end

%% Saving 

%Save coordinates in cartesian and spherical to .mat files and trajectories
%to .txt files appropriate for Lidar suitable for step&stare mode
if save_flag==1
    % create destination folder if it doesn't exist
    if ~exist(save_folder,'dir')
        mkdir(save_folder); 
    end
    
    for ii=1:length(patern)
        eval(['[fileID,msg] = fopen(''' save_folder pat_names{ii} '_ss.txt''' ',' '''w''' ');']);
        %Strict for streamliner device 3 digits before and 3 digits after, azimuth elevation no space in between
%         sav_matrix=[phi(ii,1:length(patern{ii})-1)',theta(ii,1:length(patern{ii})-1)']'; %excluding the starring point
        sav_matrix =[phi(ii,1:length(patern{ii}))',theta(ii,1:length(patern{ii}))']'; %keeping all points
        sav_matrix (sav_matrix==0) = 0;  % workoaround for -00 in text file
        sav_matrix =repmat(sav_matrix,1,repeat_traj);
        fprintf(fileID,'%07.3f%07.3f\r\n',sav_matrix); 
        fclose(fileID);
    end

    for ii=1:length(patern)
%         eval(['[fileID,msg] = fopen(''' save_folder pat_names{ii} '_CSM.txt''' ',' '''w''' ');']);
%         %Strict for streamliner device 
%         sav_matrix =[A_CSM_azi_exp(ii,1:length(patern{ii})-1)',V_CSM_azi_exp(ii,1:length(patern{ii})-1)',Phi_CSM_exp(ii,1:length(patern{ii})-1)',...
%             A_CSM_ele_exp(ii,1:length(patern{ii})-1)',V_CSM_ele_exp(ii,1:length(patern{ii})-1)',Theta_CSM_exp(ii,1:length(patern{ii})-1)',ones(length(patern{ii})-1,1)*wait_time]';
%         sav_matrix (sav_matrix==0) = 0;  % workoaround for -00 in text file
%         fprintf(fileID,'A.1=%g,S.1=%g,P.1=%g*A.2=%g,S.2=%g,P.2=%g\r\nW%g\r\n',sav_matrix); %excluding the staring point
% %         fclose(fileID);
    end
    
   considered_offsetsXYZ=[off_x,off_y,off_z];
    save([save_folder 'CS_traj.mat'],'orig_patern','All_patern_cart','All_patern_spher','considered_offsetsXYZ','req_pulse_rate','traj_ss_time','traj_CSM_time')  
    
    
end
