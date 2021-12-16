clear all, clc, close all %#ok<*CLALL,*DUALC>

% Pareto front example

% Two objective optimization (minimization in this case), pareto front 
% calculation, in a bounded region. 


% Region bounds for optimization (can be seen as contraints too)
Minbnd = -1;
Maxbnd = 2;

% Objective functions
f1 = @(x)sqrt(1+x.^2);
f2 = @(x)4+2*sqrt(1+(x-1).^2);

% Step 1: figure out unconstrained minima of the two functions
[min1,fvalmin1,exit1] = fminbnd(f1,Minbnd,Maxbnd); % fminbnd finds minimum of single-variable function on fixed interval
[min2,fvalmin2,exit2] = fminbnd(f2,Minbnd,Maxbnd);

% Or for multivariable functions
initpoint = 0.5; % point to start search forminima
options = optimoptions('fminunc','Display','off','Algorithm','quasi-newton');
[min11,fvalmin11,exit11,out11] = fminunc(f1,initpoint,options); % fminunc finds minimum of unconstrained multivariable function
[min22,fvalmin22,exit22,out22] = fminunc(f2,initpoint,options);

% Step 2: Set as optimization goals the unconstrained minima
optgoals = [fvalmin1;fvalmin2];
% optgoals = [fvalmin11;fvalmin22];

% Step 3: Calculate pareto optimal front 
% 3.1 Manual calculation (applicable because it is a simpified problem)
     %Create the criterion space for all possible weights
N = 100; % discretization
spaceN = linspace(Minbnd,Maxbnd,N+1); % bounded design space discretization
factLim = 0.001;
weight_vec = [[0:1/(factLim*N):N/(factLim*N)];[N/(factLim*N):-1/(factLim*N):0]]';
% weight_vec = [[0:1/(factLim*N):N/(factLim*N)];zeros(101,1)']';

for i=1:length(spaceN)
    xval  = spaceN(i);
    fval1 = f1(xval);
    fval2 = f2(xval);   
    calcvalTot(i) = (fval1-optgoals(1))/optgoals(1) + (fval2-optgoals(2))/optgoals(2); 
    calcvalf1(i) = (fval1-optgoals(1))/optgoals(1) ; 
    calcvalf2(i) = (fval2-optgoals(2))/optgoals(2);
%     for ii = 1:length(weight_vec)
% %         calcval(i,ii) =  weight_vec(ii,1)*fval1/optgoals(1) + weight_vec(ii,2)*fval2/optgoals(2);    
%         calcval(i,ii) =  weight_vec(ii,1)*fval1/optgoals(1) + weight_vec(ii,2)*fval2/optgoals(2);    
%     end
%     optval(i,:) = min(calcval(i,:));  
%     indexN(i,:) =find(calcval(i,:)==min(calcval(i,:)));   
%     f1valPar = 
end

% 3.2 Using fgoalattain native function
options = optimoptions('fgoalattain','Display','off');
% f=@(x)deal(f1(x),f2(x));
fun = @simple_mult;
% function f = simple_mult(x)
% f(:,1) = sqrt(1+x.^2);
% f(:,2) = 4 + 2*sqrt(1+(x-1).^2);
for i=1:length(spaceN)
     weight = weight_vec(i,:);
%     weight = [2,2];
    [x_att(i,:),f_att(i,:),attainfactor(i,:),exitflag_att(i,:)] = fgoalattain(fun,initpoint,optgoals,weight,...
        [],[],[],[],[],[],[],options);       
end
% Step 4: visualize

% Design space
figure
hold on 
plot(spaceN,f1(spaceN),spaceN,f2(spaceN),'LineWidth',2);
plot([0,0],[0,8],'g--','LineWidth',2);
plot([1,1],[0,8],'g--','LineWidth',2);
hold off
xlabel('Variable')
ylabel('Objective')
legend({'f1' 'f2'})
title('Design space')
grid on


%Criterion space
figure
plot(f_att(:,1),f_att(:,2),'k.',f1(spaceN),f2(spaceN),'r.');
xlabel('f_1')
ylabel('f_2')
legend({'Using fgoalattain' 'Manual implementation'})
title('Pareto front in criterion space')
grid on

figure
hold on
plot(spaceN,calcvalTot,spaceN,calcvalf1,spaceN,calcvalf2,'LineWidth',2),grid on
legend ({'Total Cost' 'f1' 'f2'})
plot([0,0],[0,1.5],'g--','LineWidth',2);
plot([1,1],[0,1.5],'g--','LineWidth',2);
xlabel('design space')
ylabel('cost')
title('Pareto front calculated manually')

