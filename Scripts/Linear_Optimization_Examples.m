clearvars
close all

%% MAtlab example on simplex
tic
A = [1 1
    1 1/4
    1 -1
    -1/4 -1
    -1 -1
    -1 1];

b = [2 1 2 1 -1 2];
f = [-1 -1/3];
x_mat = linprog(f,A,b);
y_mat = f*x_mat;
T_mat = toc;


%% Vasilis implementation with brute force
% Stupid way but works...
% first find a rough space and then refine...

tic
incr= 0.005;
x1 = 0:incr:10;
x2 = 0:incr:10;
for i =1:length(x1)
    for ii =1:length(x2)        
%         y(i,ii)= -x1(i) -(1/3)*x2(ii); %#ok<*SAGROW> %just for reference    
        % apply constrains to the values selected by removing values that are out of bounds
        if x1(i)+x2(ii)>2 || x1(i)+x2(ii)/4>1 || x1(i)-x2(ii)>2 || -x1(i)/4-x2(ii)>1 || -x1(i)-x2(ii)>-1 || -x1(i)+x2(ii)>2 
            y1(i,ii)= nan;
        else
            y1(i,ii)= -x1(i) -(1/3)*x2(ii);            
        end       
    end
end
% [M,I]=min(y1,[],'omitnan','all','linear');
y_vas = min(y1,[],'all');
[row,col]=find(y1==y_vas);

x_vas = [x1(row);x2(col)];
T_vas = toc;
        

%% using fmincon form matlab can be more general!!!

tic
A1 = [1 1
    1 1/4
    1 -1
    -1/4 -1
    -1 -1
    -1 1];
b1 = [2 1 2 1 -1 2];
x0 = [0,10];

fun = @(x)-x(1)-(1/3)*x(2);
x_fmincon = fmincon(fun,x0,A1,b1) ;
y_fmincon = -x_fmincon(1)-(1/3)*x_fmincon(2)       ;

T_fmincon = toc;