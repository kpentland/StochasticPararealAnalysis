%% Figure 1

%% State-independent error bounds vs. numerics (B < 1)
% Here we plot the errors from SParareal with state-independent
% perturbations against the (linear and superlinear) theoretical error
% bounds. We increase the parameter q to examine the effect of increasing
% the abslute 'size' of the perturbations.

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
sample_rule = 0;                             %sampling rule to employ (Gaussian perturbations in SParareal)
sims = 500;                                  %no. of independent sims



%%% FIRST CASE
% set desired level of noise
q = 0; 

%solve with stochastic parareal
[~,U_para,~,k,~,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule,epsilon,sims,q);

%solve using the fine solver serially (to calculate errors)
dt = (tspan(2)-tspan(1))/Nf;    t_fine = (tspan(1):dt:tspan(2));
[~,u_fine] = RK(t_fine,u0,f,F);
U_fine = u_fine(1:(Nf/N):end,:);


%Fig 1: Convergence Plot (max error between SParareal vs. fine solution) at each iteration

% NUMERICAL ERRORS
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



% THEORETICAL ERRORS

%stability function for coarse solver (for RK1 solver)
DT = (tspan(2)-tspan(1))/N;                      %time slice width
Lg = max(sum(abs(eye(n,n) + DT*Q),2));           %Lipschitz constant for G
p = 1;                                           %GTE of the RK method

% k = 0 max error (choose to be equal to max numerical error at k = 0)
e0 = exper_bounds(1,1);

%constants
C1 = max(sum(abs(expm(Q*DT) - (eye(n,n) + DT*Q)),2)/(DT^(p+1)));
C2 = 1;

% bound parameters
e1 = DT;
e2 = 1;
e3 = 1/DT;
A = (1 + (1/e1) + (1/e2))*(C1^2)*(DT^((2*p)+2));
B = (1 + e1 + (1/e3))*(Lg^2);
L = (1 + e2 + e3)*(C2^2)*(DT^((2*q)+1));
D = A*e0;


% SUPERLINEAR BOUND
tempkN = zeros(N+1,N+1);

% k = 0
tempkN(:,1) = e0;

% k = 1
tempk1 = zeros(1,N+1);
for j = 2:N
    tempk1(j+1) = tempk1(j) + B^(j-2);
end
tempkN(:,2) = max(e0)*A*tempk1';

% k => 2
for nn = (3:N)
    for kk = 2:(nn-1)
        temp1 = 0;
        for l = 0:(nn-kk)
            temp1 = temp1 + nchoosek(l+kk-1,l)*(B^l);
        end
        
        temp2 = 0;
        for j = 0:(kk-2)
            for l = 0:nn-(j+1)
                temp2 = temp2 + nchoosek(l+j,l)*(A^j)*(B^l);
            end
        end
        
        tempkN(nn+1,kk+1) = D*(A^(kk-1))*temp1 + L*temp2;
    end
end

% take maximum over n to get max error at each k
superlin_bound = max(tempkN);  superlin_bound(superlin_bound==0) = 10^(-100);





% LINEAR BOUND
lin_bound = zeros(1,N+1);
lin_bound(1) = superlin_bound(1);
lin_bound(2) = superlin_bound(2);

AA = (A/(1-B)).^(0:N+1);
for i = 3:N+1
    lin_bound(i) = lin_bound(2)*AA(i-1) + (L/(1-B))*sum(AA(1:(i-1)));
end
lin_bound(lin_bound==0) = 10^(-100);


% Plot theory vs. numerics
figure(1)
hold on
plot((0:N),lin_bound,'-*r','LineWidth',1.2)
plot((0:N),superlin_bound,'-xb','LineWidth',1.2)
plot((0:N),exper_bounds,'-ok','LineWidth',1.2)
plot([0,N],DT^((2*q)+1)*[1,1],'--k','LineWidth',1)
hold off
box on; grid on;

xlabel('$k$','interpreter','latex');
xlim([0 18])
xticks((0:2:21))

ylabel('$\hat{e}^k$','interpreter','latex');
ylim([10^(-25) 10^(5)]); 
yticks(10.^[-25,-20,-15,-10,-5,0,5])
set(gca,'yscale','log')
set(gca,'FontSize',14)





%%% SECOND CASE
% set desired level of noise
q = 5;

%solve with stochastic parareal
[~,U_para,~,k,~,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule,epsilon,sims,q);


%solve using the fine solver serially (to calculate errors)
n = length(u0);
dt = (tspan(2)-tspan(1))/Nf;    t_fine = (tspan(1):dt:tspan(2));
[~,u_fine] = RK(t_fine,u0,f,F);
U_fine = u_fine(1:(Nf/N):end,:);


%Fig 1: Convergence Plot (max error between SParareal vs. fine solution) at each iteration

% NUMERICAL ERRORS
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



% THEORETICAL ERRORS

%stability function for coarse solver (for RK1 solver)
DT = (tspan(2)-tspan(1))/N;                      %time slice width
Lg = max(sum(abs(eye(n,n) + DT*Q),2));           %Lipschitz constant for G
p = 1;                                           %GTE of the RK method

% k = 0 max error (choose to be equal to max numerical error at k = 0)
e0 = exper_bounds(1,1);

%constants
C1 = max(sum(abs(expm(Q*DT) - (eye(n,n) + DT*Q)),2))/(DT^(p+1));
C2 = 1;

% bound parameters
e1 = DT;
e2 = 1;
e3 = 1/DT;
A = (1 + (1/e1) + (1/e2))*(C1^2)*(DT^((2*p)+2));
B = (1 + e1 + (1/e3))*(Lg^2);
L = (1 + e2 + e3)*(C2^2)*(DT^((2*q)+1));
D = A*e0;


% SUPERLINEAR BOUND
tempkN = zeros(N+1,N+1);

% k = 0
tempkN(:,1) = e0;

% k = 1
tempk1 = zeros(1,N+1);
for j = 2:N
    tempk1(j+1) = tempk1(j) + B^(j-2);
end
tempkN(:,2) = max(e0)*A*tempk1';

% k => 2
for nn = (3:N)
    for kk = 2:(nn-1)
        temp1 = 0;
        for l = 0:(nn-kk)
            temp1 = temp1 + nchoosek(l+kk-1,l)*(B^l);
        end
        
        temp2 = 0;
        for j = 0:(kk-2)
            for l = 0:nn-(j+1)
                temp2 = temp2 + nchoosek(l+j,l)*(A^j)*(B^l);
            end
        end
        
        tempkN(nn+1,kk+1) = D*(A^(kk-1))*temp1 + L*temp2;
    end
end

% take maximum over n to get max error at each k
superlin_bound = max(tempkN);  superlin_bound(superlin_bound==0) = 10^(-100);





% LINEAR BOUND
lin_bound = zeros(1,N+1);
lin_bound(1) = superlin_bound(1);
lin_bound(2) = superlin_bound(2);

AA = (A/(1-B)).^(0:N+1);
for i = 3:N+1
    lin_bound(i) = lin_bound(2)*AA(i-1) + (L/(1-B))*sum(AA(1:(i-1)));
end
lin_bound(lin_bound==0) = 10^(-100);


% Plot theory vs. numerics
figure(2)
hold on
plot((0:N),lin_bound,'-*r','LineWidth',1.2)
plot((0:N),superlin_bound,'-xb','LineWidth',1.2)
plot((0:N),exper_bounds,'-ok','LineWidth',1.2)
plot([0,N],DT^((2*q)+1)*[1,1],'--k','LineWidth',1)

hold off
box on; grid on;

xlabel('$k$','interpreter','latex');
xlim([0 18])
xticks((0:2:21))

ylim([10^(-25) 10^(5)]); 
yticks(10.^[-25,-20,-15,-10,-5,0,5])
set(gca,'yticklabel',[])
set(gca,'yscale','log')
set(gca,'FontSize',14)




%%% THIRD CASE
% set desired level of noise
q = 10;

%solve with stochastic parareal
[~,U_para,~,k,~,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule,epsilon,sims,q);


%solve using the fine solver serially (to calculate errors)
dt = (tspan(2)-tspan(1))/Nf;    t_fine = (tspan(1):dt:tspan(2));
[~,u_fine] = RK(t_fine,u0,f,F);
U_fine = u_fine(1:(Nf/N):end,:);


%Fig 1: Convergence Plot (max error between SParareal vs. fine solution) at each iteration

% NUMERICAL ERRORS
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



% THEORETICAL ERRORS

%stability function for coarse solver (for RK1 solver)
DT = (tspan(2)-tspan(1))/N;                      %time slice width
Lg = max(sum(abs(eye(n,n) + DT*Q),2));           %Lipschitz constant for G
p = 1;                                           %GTE of the RK method

% k = 0 max error (choose to be equal to max numerical error at k = 0)
e0 = exper_bounds(1,1);

%constants
C1 = max(sum(abs(expm(Q*DT) - (eye(n,n) + DT*Q)),2))/(DT^(p+1));
C2 = 1;

% bound parameters
e1 = DT;
e2 = 1;
e3 = 1/DT;
A = (1 + (1/e1) + (1/e2))*(C1^2)*(DT^((2*p)+2));
B = (1 + e1 + (1/e3))*(Lg^2);
L = (1 + e2 + e3)*(C2^2)*(DT^((2*q)+1));
D = A*e0;


% SUPERLINEAR BOUND (B<1)
tempkN = zeros(N+1,N+1);

% k = 0
tempkN(:,1) = e0;

% k = 1
tempk1 = zeros(1,N+1);
for j = 2:N
    tempk1(j+1) = tempk1(j) + B^(j-2);
end
tempkN(:,2) = max(e0)*A*tempk1';

% k => 2
for nn = (3:N)
    for kk = 2:(nn-1)
        temp1 = 0;
        for l = 0:(nn-kk)
            temp1 = temp1 + nchoosek(l+kk-1,l)*(B^l);
        end
        
        temp2 = 0;
        for j = 0:(kk-2)
            for l = 0:nn-(j+1)
                temp2 = temp2 + nchoosek(l+j,l)*(A^j)*(B^l);
            end
        end
        
        tempkN(nn+1,kk+1) = D*(A^(kk-1))*temp1 + L*temp2;
    end
end

% take maximum over n to get max error at each k
superlin_bound = max(tempkN);  superlin_bound(superlin_bound==0) = 10^(-100);





% LINEAR BOUND (B<1)
lin_bound = zeros(1,N+1);
lin_bound(1) = superlin_bound(1);
lin_bound(2) = superlin_bound(2);

AA = (A/(1-B)).^(0:N+1);
for i = 3:N+1
    lin_bound(i) = lin_bound(2)*AA(i-1) + (L/(1-B))*sum(AA(1:(i-1)));
end
lin_bound(lin_bound==0) = 10^(-100);


% Plot theory vs. numerics
figure(3)
hold on
plot((0:N),lin_bound,'-*r','LineWidth',1.2)
plot((0:N),superlin_bound,'-xb','LineWidth',1.2)
plot((0:N),exper_bounds,'-ok','LineWidth',1.2)
% plot((2:N),loose(2:N),'--k','LineWidth',1.2)
% plot((2:N),0.1*max(DT.^(((2*p)+1)*(2:N)),ones(1,N-1)*DT^(2*q)),'--k','LineWidth',1.2)
plot([0,N],DT^((2*q)+1)*[1,1],'--k','LineWidth',1)
% plot((0:N),parareal_bound,'-om','LineWidth',1.2)
hold off
box on; grid on;

xlabel('$k$','interpreter','latex');
xlim([0 18])
xticks((0:2:21))

ylim([10^(-25) 10^(5)]); 
yticks(10.^[-25,-20,-15,-10,-5,0,5])
set(gca,'yticklabel',[])
set(gca,'yscale','log')

set(gca,'FontSize',14)
legend({'Linear bound','Superlinear bound','Numerical error'},'interpreter','latex','location','northeast')
