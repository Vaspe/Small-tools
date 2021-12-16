%example script for optimization
%based on  https://www.inf.ethz.ch/personal/fukudak/lect/opt2011/aopt11note1.pdf
%VP 10.1.2019

clc,clear all,close all

 %% brute force
tic
x1=[0:1:5];% constraint x1>=0 and x1+x2<7/3
x2=[0:0.0001:sqrt(7/3)];
for i=1:length(x1)
   ix1=x1(i);
   for j=1:length(x2)
      ix2=x2(j);
      if ix1+ix2>7/3 %constraint
         res(i,j)=NaN;
      else
         res(i,j)=3*ix1^2+ix2^2;
      end
             
   end
end
% [row,col]=find(res<=7/3);
% feasib_sol = res(find(res<=7/3));
% [optsol,idx_sol] = max(res(row,col) );
[opt_sol,Ind]=max(res(:),[],'omitnan');
[x1optind,x2optind]=ind2sub(size(res),Ind);
t1=toc;

disp(['Optimal value with brute force is ' num2str(opt_sol) ' with values x1,x2 equal to '...
    num2str(x1(x1optind)) ' and ' num2str(x2(x2optind)) '.']);
disp(['Time brute force: ' num2str(t1) ' s']);
[X1,X2]=meshgrid(x1,x2);
figure,
mesh(x2,x1,res)


%% New problem p9 brute force
% x1 red,x2 white, x3 rose units sold / price per type K1,K2,K3, max supply
% of grapes per day S1,S2,S3

%Constraints
% PN<=4 % pinot grapes noir per day
% GM<=8% gamay grapes per day
% CH<=6 % chaselas grapes per day
% x1 requires 2PN or 1GM, x2 3CH, x3 2SM or 1CH


% maximize profit
K1=3;
K2=4;
K3=2;

S1=4;
S2=8;
S3=6;

% Constraints eqs
%x1/2<=S1
%x1/1+x3/2<=S2
%x2/3+x3/1<=S3


% objective function
%max x1*K1+x2*K2+x3*K3 

x1=[0:0.05:10];
x2=[0:0.05:2];
x3=[0:0.05:10];
count=1; 
for i=1:length(x1)
   ix1=x1(i);
   for j=1:length(x2)
      ix2=x2(j);
       for k=1:length(x3)
          ix3=x3(k);
          if ix1*2<=S1 && (ix1*1)+(ix3*2)<=S2 && (ix2*3)+(ix3*1)<=S3
              ResW(i,j,k)= ix1*K1+ix2*K2+ix3*K3; 
              xx(count)=ix1;
              yy(count)=ix2;
              zz(count)=ix3;
              count=count+1;
          else
              ResW(i,j,k)=NaN;
          end
       end
   end
end
[optsol,idx_sol] = max(ResW(:),[],'omitnan');
[x1optind,x2optind,x3optind]=ind2sub(size(ResW),idx_sol);
disp(['Optimal wine profit with brute force is ' num2str(optsol) ' with values x1,x2,x3 equal to '...
num2str(x1(x1optind)) ', ' num2str(x2(x2optind)) ' and ' num2str(x3(x3optind))]);

% figure,
% plot3(xx,yy,zz,'*')

shp = alphaShape(xx',yy',zz');
figure
plot(shp)








