clc,clear all,close all

% testing mass damper system according to https://www.youtube.com/watch?v=l3Ig97bbxdg
t=0:0.00001:5;
C1=0.5;
C2=0.5;
m=100;
k=950;

%underdamping
m1=m;
k1=k;
d1=sqrt(4*m1*k1)-20.003 ; %4*m
omega1= sqrt((k1/m1) -((d1^(2))/(4*m1^2) ))/(2*pi); %Hz
NaturalFreq = sqrt(k1/m1)/(2*pi)
CriticalDamping=2*sqrt(k1*m1) %#ok<*NOPTS>
DampigRatio1 = d1/( 2*sqrt(k1*m1) )

%critical
m2=m;
k2=k;
d2=sqrt(4*m2*k2);
omega2= sqrt((k2/m2) -(d2^(2)/(4*m2^2) ))/(2*pi);
DampigRatio2 = d2/(2*sqrt(k2/m2));

%overdamping
d3=d2+0.00000003;
m3=m;
k3=k;
omega3= sqrt((k3/m3) -(d3^(2)/(4*m3^2) ))/(2*pi);
DampigRatio3 = d3/(2*sqrt(k3/m3));

%undamped 
d4=0;
m4=m;
k4=k;
omega4= sqrt((k4/m4) -(d4^(2)/(4*m4^2) ))/(2*pi);
DampigRatio4 = d4/(2*sqrt(k4/m4));
poles4_1=omega4;
poles4_2=-omega4;

% solution for the mass spring damper second order linear system
y1=C1.*exp(-( d1/(2*m1)-( sqrt(d1^(2)-4*m1*k1)*(2*m1) )).*t)+C2.*exp(-( d1/(2*m1)+( sqrt(d1^(2)-4*m1*k1)*(2*m1) )).*t);
y2=C1.*exp(-( d2/(2*m2)-( sqrt(d2^(2)-4*m2*k2)*(2*m2) )).*t)+C2.*exp(-( d2/(2*m2)+( sqrt(d2^(2)-4*m2*k2)*(2*m2) )).*t);
y3=C1.*exp(-( d3/(2*m3)-( sqrt(d3^(2)-4*m3*k3)*(2*m3) )).*t)+C2.*exp(-( d3/(2*m3)+( sqrt(d3^(2)-4*m3*k3)*(2*m3) )).*t);
y4=C1.*exp(-( d4/(2*m4)-( sqrt(d4^(2)-4*m4*k4)*(2*m4) )).*t)+C2.*exp(-( d4/(2*m4)+( sqrt(d4^(2)-4*m4*k4)*(2*m4) )).*t);

figure
plot(t,y1,t,y2,t,y3,t,y4,'Linewidth',2)
legend({'Underdamping ' 'Critical Damping' 'Over damping' 'Undamped'})
grid on

figure
plot([0,0],[poles4_1,poles4_2],'x')
xlabel ('Re(x)')
ylabel ('Im(x)')
title ('Poles for undamped case')



%%