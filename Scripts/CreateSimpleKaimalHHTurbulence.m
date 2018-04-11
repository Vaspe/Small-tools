clc,clear all,close all

%%Atmospheric turbulence based on Kaimal spectrum
dt    = 0.0125;
timeLength  = 1060;
I     = 0.035; % turbulence intensity [-] typical value chosen
Vav   = 16; % average wind velocity [m/s]
l     = 340.2; % turbulence length scale for z>60m [m]
fhighcut_wind       = 35;      % Hertz
ConstantVShear      = 0.14;   % vertical shear with exponential law;
CreateWindfieldFile = 1;      % Choose hether to create XYZ file
HubHeight           = 115.63; % Hub height in m 
RotorRadius         = 89.2;   %including hub 

filedirectory = 'D:\code_witlis_FAST\witlisCurrent\trunkSC\Wind\';
name      = 'TurbwindHHKaimal_';

%%
t             = [0:dt:timeLength]';
df            = 1/timeLength; % frequency descretization because it comes from observation time in sec
fwind         = df:df:fhighcut_wind; % frequency domain
Swind         = 4*(I^2)*Vav*l./(1+6*(fwind*l/Vav)).^(5/3) ; % Kaimal spectrum p 32 pdf wind climate lecture
bp            = sqrt(2*Swind*df) ;             % bp for amlitude in order to create time series with fourrier p 32 pdf wind climate lecture
M             = length(t);
eps           = 2*pi*rand(size(fwind)) ;               %random phase [0,2pi] to create the stochastic wind
y             = M*real(ifft(pad2(bp.*exp(1i*eps),M))) ; % inverse Fourier to go to time domain!
Vts           = [Vav+y]' ;                              % time series of velocity for 10 min


plot (t,Vts)  

%% Create wind field file complied to windfield.mat files 
if CreateWindfieldFile==1

    
    VShearCol     = ConstantVShear*ones(length(t),1); % create shear for file 
    
    dy = 8;  % discretization of grid points in y direction
    dz = 8;  % discretization of grid points in z direction
    nt = length(t); % total time steps (can be seen as x direction)
    GridLength = round(RotorRadius+20); % the distance from 0 -rotor center to end of yz slice (square shape)
    GridLength = GridLength+dy-rem(GridLength,dy);
    
    z  = -GridLength:dz:GridLength;  %
    y  = -GridLength:dy:GridLength;  %     
    ny = length(y); % nymber of grid points in y direction
    nz = length(z);  %number of grid points in z direction 
%     [T, Y1, Z1] = meshgrid(t,y,z);
    for i = 1:ny
        Y(i,:)  = z;  
    end
    
    for i = 1:nz
        Z(:,i)  = y;  
    end
    
    u=zeros(nz,length(t),ny);
    %(nz x nt x ny)
    for i = 1:length(Vts)         
        Zdiruvec = Vts(i)*((z+HubHeight)./HubHeight).^ConstantVShear ;
        for ii = 1:length(y)    
            u(ii,i,:) = Zdiruvec;     
        end       
    end
        
    
    %Assign values
    windfield.Uref     = Vav;
    windfield.T_offset = 0; 
    windfield.dt       = dt; 
    windfield.ny       = ny;
    windfield.nz       = nz;         
    windfield.v        = 0;
    windfield.w        = 0;
    windfield.u   = u;    
    windfield.grid.ny  = ny;    
    windfield.grid.nz  = nz;
    windfield.grid.dt  = dt;    
    windfield.grid.dy  = dy;    
    windfield.grid.dz  = dz;        
    windfield.grid.nt  = nt;       
    windfield.grid.t   = t;
    windfield.grid.y   = y;
    windfield.grid.z   = z;  
    windfield.grid.Y   = Y;
    windfield.grid.Z   = Z;       

    
end

%% Save files

filename    = [filedirectory name num2str(Vav) ] ;    
zerocolumns = zeros(length(t),1);
finalvar = [t,Vts,zerocolumns,zerocolumns,VShearCol,zerocolumns,zerocolumns,zerocolumns];

dlmwrite([filename '.wnd'],finalvar,'precision','%.4f','delimiter','\t');
save(filename,'windfield')
