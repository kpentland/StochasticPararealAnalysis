%% Figure 4

%% State-dependent (sampling rules 1, 2, 3, and 4) error bounds vs. numerics (B < 1)
% Here we plot the errors from SParareal with state-dependent
% perturbations against the (linear) theoretical error
% bounds as well as the 'numerical' error bounds. 

clear; close all; clc

%Inputs:
n = 100;                                     % dimension of the ODE system
Q = -2*eye(n,n) + (1/n)*(-1 + 2*rand(n,n));  % generate the matrix Q (i.e. du/dt = Qu)
u0 = -5 + 10*rand(1,n);                      % initial condition
f = @(t,u)( Q*u );                           %function handle for ODE
tspan = [0,2];                               %time interval
N = 20;                                      %no. of time sub-intervals step
Ng = N;                                      %no. of coarse steps
Nf = Ng;                                     %no. of fine steps
epsilon = 10^(-10000);                       %error tolerance
G = 'RK1';                                   %coarse solver
F = 'ExactLinearSystem';                     %fine solver
sims = 500;                                  %no. of independent sims


%%% SAMPLING RULES 2 and 4
sample_rule = 2;           
sample_rule2 = 4;

%solve with stochastic parareal
[~,U_para,~,~,~,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule,epsilon,sims,0);
[~,U_para2,~,~,~,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule2,epsilon,sims,0);


%solve using the fine solver serially (to calculate errors)
dt = (tspan(2)-tspan(1))/Nf;    t_fine = (tspan(1):dt:tspan(2));
[~,u_fine] = RK(t_fine,u0,f,F);
U_fine = u_fine(1:(Nf/N):end,:);



% NUMERICAL ERRORS (first sampling rule)
tempU = cell(sims,1);
UU_fine = repmat(U_fine,1,N+1);
for j = 1:sims
    temp = U_para{j,1} - UU_fine;
    
    temp1 = zeros(N+1,N+1);
    for i = 1:N+1
        ind = (n*(i-1)+1:n*i);
        temp1(:,i) = ( vecnorm(temp(:,ind),inf,2) ).^2; 
    end
    tempU{j,1} = temp1;
end
temppp = mean(cat(3,tempU{:}),3);
exper_bounds = max(temppp);



% NUMERICAL ERRORS (second sampling rule)
tempU2 = cell(sims,1);
UU_fine = repmat(U_fine,1,N+1);
for j = 1:sims
    temp = U_para2{j,1} - UU_fine;
    
    temp1 = zeros(N+1,N+1);
    for i = 1:N+1
        ind = (n*(i-1)+1:n*i);
        temp1(:,i) = ( vecnorm(temp(:,ind),inf,2) ).^2; 
    end
    tempU2{j,1} = temp1;
end
temppp2 = mean(cat(3,tempU2{:}),3);
exper_bounds2 = max(temppp2);



% THEORETICAL ERRORS

%stability function for coarse solver (for RK1 solver)
DT = (tspan(2)-tspan(1))/N;                      %time slice width
Lg = max(sum(abs(eye(n,n) + DT*Q),2));                   %stability function
p = 1;                                           %GTE of the RK method

% k = 0 max error (choose to be equal to max numerical error at k = 0)
e0 = exper_bounds(1,1);

%constants
C1 = max(sum(abs(expm(Q*DT) - (eye(n,n) + DT*Q)),2)/(DT^(p+1)));
C2 = 1;

% bound parameters
ee1 = DT;
ee2 = 1;
ee3 = 1/DT;
ee4 = 1/DT;
A = (C1^2)*(DT^((2*p)+2))*(1 + (1/ee1) + (1/ee2));
B = (Lg^2)*(1 + ee1 + (1/ee3));
D1 = (C1^2)*(DT^((2*p)+2))*(Lg^2)*(1 + ee4);
D2 = (C1^2)*(DT^((2*p)+2))*(Lg^2)*(1 + (1/ee4));


% numerical (superlinear) bound
superlin_bound = NaN(N+1,N+1);   %n by k
superlin_bound(logical(eye(N+1))) = 0;
superlin_bound(1,2:end) = e0;
% k = 1
BB = B.^(0:N+1);
for nn = 3:N+1
    superlin_bound(2,nn) = e0*A*sum(BB(1:nn-2));
end
% k => 2
for kk = (3:N+1)
    for nn = kk+1:N+1
        superlin_bound(kk,nn) = A*superlin_bound(kk-1,nn-1) + B*superlin_bound(kk,nn-1) + D1*superlin_bound(kk-1,nn-2) + D2*superlin_bound(kk-2,nn-2);
    end
end
superlin_bound = max(superlin_bound,[],2);  superlin_bound(superlin_bound==0) = 10^(-100);


% Linear bound calculations
const = (A + D1 + sqrt((A + D1)^2 + 4*D2*(1-B)))/(2*(1-B));

lin_bound = zeros(N+1,1);
% k = 0
lin_bound(1) = e0;

% k = 1
BB = B.^(0:N+1);
temp = zeros(N+1,1);
for i = 3:N+1
    temp(i,1) = e0*A*sum(BB(1:i-2));
end
lin_bound(2) = max(temp);

% k > 1
lin_bound(3:end,1) = e0*const.^(2:N);


figure(1)
hold on
plot((0:N),lin_bound,'-*r','LineWidth',1.2)
plot((0:N),superlin_bound,'-xb','LineWidth',1.2)
plot((0:N),exper_bounds,'-','color',[0.6350 0.0780 0.1840],'marker','s','LineWidth',1.2,'MarkerSize',9)
plot((0:N),exper_bounds2,'--','color','g','marker','+','LineWidth',1.2,'MarkerSize',5)
% plot((0:N),lin_bound2,'-*y','LineWidth',1.2)
hold off
xlabel('$k$','interpreter','latex');
ylabel('$\hat{e}^k$','interpreter','latex');
box on; grid on;
xlim([0 18])
xticks((0:2:21))
ylim([10^(-15) 10^(2)]); 
yticks(10.^[-25,-20,-15,-10,-5,0,5])
set(gca,'yscale','log')
legend({'Linear bound',' ``Numerical" bound','Numerical error (SR2)','Numerical error (SR4)'},'interpreter','latex','location','northeast')
set(gca,'FontSize',12)


%%% SAMPLING RULES 1 and 3
sample_rule = 1;
sample_rule2 = 3;


%solve with stochastic parareal
[~,U_para,~,~,~,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule,epsilon,sims,0);
[~,U_para2,~,~,~,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule2,epsilon,sims,0);


%solve using the fine solver serially (to calculate errors)
dt = (tspan(2)-tspan(1))/Nf;    t_fine = (tspan(1):dt:tspan(2));
[~,u_fine] = RK(t_fine,u0,f,F);
U_fine = u_fine(1:(Nf/N):end,:);



% NUMERICAL ERRORS (first sampling rule)
tempU = cell(sims,1);
UU_fine = repmat(U_fine,1,N+1);
for j = 1:sims
    temp = U_para{j,1} - UU_fine;
    
    temp1 = zeros(N+1,N+1);
    for i = 1:N+1
        ind = (n*(i-1)+1:n*i);
        temp1(:,i) = ( vecnorm(temp(:,ind),inf,2) ).^2; 
    end
    tempU{j,1} = temp1;
end
temppp = mean(cat(3,tempU{:}),3);
exper_bounds = max(temppp);



% NUMERICAL ERRORS (second sampling rule)
tempU2 = cell(sims,1);
UU_fine = repmat(U_fine,1,N+1);
for j = 1:sims
    temp = U_para2{j,1} - UU_fine;
    
    temp1 = zeros(N+1,N+1);
    for i = 1:N+1
        ind = (n*(i-1)+1:n*i);
        temp1(:,i) = ( vecnorm(temp(:,ind),inf,2) ).^2; 
    end
    tempU2{j,1} = temp1;
end
temppp2 = mean(cat(3,tempU2{:}),3);
exper_bounds2 = max(temppp2);



% THEORETICAL ERRORS

%stability function for coarse solver (for RK1 solver)
DT = (tspan(2)-tspan(1))/N;                      %time slice width
Lg = max(sum(abs(eye(n,n) + DT*Q),2));                   %stability function
p = 1;                                           %GTE of the RK method
Lf = max(sum(abs(Q),2));

% k = 0 max error (choose to be equal to max numerical error at k = 0)
e0 = exper_bounds(1,1);

%constants
C1 = max(sum(abs(expm(Q*DT) - (eye(n,n) + DT*Q)),2)/(DT^(p+1)));
C2 = 1;

% bound parameters
ee1 = DT;
ee2 = 1;
ee3 = 1/DT;
ee4 = 1;
ee5 = 1;
A = (C1^2)*(DT^((2*p)+2))*(1 + (1/ee1) + (1/ee2));
B = (Lg^2)*(1 + ee1 + (1/ee3));
D1 = 2*(C1^2)*(DT^((2*p)+2))*(Lg^2)*(1 + (1/DT));
D2 = 2*(C1^2)*(DT^((2*p)+2))*((Lg^2)*(1 + DT) + 2*(Lf^2));
D3 = 4*(C1^2)*(DT^((2*p)+2));

% numerical (superlinear) bound
superlin_bound = NaN(N+1,N+1);   %n by k
superlin_bound(logical(eye(N+1))) = 0;
superlin_bound(1,2:end) = e0;
% k = 1
BB = B.^(0:N+1);
for nn = 3:N+1
    superlin_bound(2,nn) = e0*A*sum(BB(1:nn-2));
end
% k => 2
for kk = (3:N+1)
    for nn = kk+1:N+1
        superlin_bound(kk,nn) = (A+D3)*superlin_bound(kk-1,nn-1) + B*superlin_bound(kk,nn-1) + D1*superlin_bound(kk-1,nn-2) + D2*superlin_bound(kk-2,nn-2);
    end
end
superlin_bound = max(superlin_bound,[],2);  superlin_bound(superlin_bound==0) = 10^(-100);



% Linear bound calculations
const = (A + D1 + D3 + sqrt((A + D1 + D3)^2 + 4*D2*(1-B)))/(2*(1-B));

lin_bound = zeros(N+1,1);
% k = 0
lin_bound(1) = e0;

% k = 1
BB = B.^(0:N+1);
temp = zeros(N+1,1);
for i = 3:N+1
    temp(i,1) = e0*A*sum(BB(1:i-2));
end
lin_bound(2) = max(temp);

% k > 1
lin_bound(3:end,1) = e0*const.^(2:N);



figure(2)
hold on
plot((0:N),lin_bound,'-*r','LineWidth',1.2)
plot((0:N),superlin_bound,'-xb','LineWidth',1.2)
plot((0:N),exper_bounds,'-','color',[0 0.4470 0.7410],'marker','o','LineWidth',1.2,'MarkerSize',9)
plot((0:N),exper_bounds2,'--','color','m','marker',"pentagram",'LineWidth',1.2,'MarkerSize',5)
hold off
xlabel('$k$','interpreter','latex');
ylabel('$\hat{e}^k$','interpreter','latex');
box on; grid on;
xlim([0 18])
xticks((0:2:21))
ylim([10^(-15) 10^(2)]); 
yticks(10.^[-25,-20,-15,-10,-5,0,5])
set(gca,'yscale','log')
legend({'Linear bound',' ``Numerical" bound','Numerical error (SR1)','Numerical error (SR3)'},'interpreter','latex','location','northeast')
set(gca,'FontSize',12)