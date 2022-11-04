%% Figure 5

%% Iteration k vs. tolerance epsilon
% Here we plot the exepected number of iterations taken to reach a given
% stopping tolerance. We do this for each of the sampling rules and for
% SParareal with different levels of Gaussian noise. 

% This section takes a while to run, so instead, run the next section which
% plots the results using pre-run values stored in 'sol.mat'.

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
G = 'RK1';                                   %coarse solver
F = 'ExactLinearSystem';                     %fine solver
sample_rule = 0;                             %sampling rule to employ (Gaussian perturbations in SParareal)
sims = 500;                                  %no. of independent sims


% error tolerances to test
epsilon = 10.^((-15:1:0));

% run SParareal for different values of q and store the number of
% iterations is takes on average
q = 0;
k1 = zeros(length(epsilon),1);
for i = 1:length(epsilon)
    [~,~,~,k,~,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule,epsilon(i),sims,q);   
    k1(i) = mean(k);
end


q = 5;
k2 = zeros(length(epsilon),1);
for i = 1:length(epsilon)
    [~,~,~,k,~,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule,epsilon(i),sims,q);   
    k2(i) = mean(k);
end


q = 10;
k3 = zeros(length(epsilon),1);
for i = 1:length(epsilon)
    [~,~,~,k,~,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule,epsilon(i),sims,q);   
    k3(i) = mean(k);
end

q = 25;
k4 = zeros(length(epsilon),1);
for i = 1:length(epsilon)
    [~,~,~,k,~,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule,epsilon(i),sims,q);   
    k4(i) = mean(k);
end

% do the same thing but run using one of the sampling rules 
q = 0;
sample_rule = 1;                        %sampling rule to employ
k5 = zeros(length(epsilon),1);
for i = 1:length(epsilon)
    [~,~,~,k,~,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule,epsilon(i),sims,q);   
    k5(i) = mean(k);
end

sample_rule = 2;                        %sampling rule to employ
k6 = zeros(length(epsilon),1);
for i = 1:length(epsilon)
    [~,~,~,k,~,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule,epsilon(i),sims,q);   
    k6(i) = mean(k);
end

sample_rule = 3;                        %sampling rule to employ
k7 = zeros(length(epsilon),1);
for i = 1:length(epsilon)
    [~,~,~,k,~,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule,epsilon(i),sims,q);   
    k7(i) = mean(k);
end

sample_rule = 4;                        %sampling rule to employ
k8 = zeros(length(epsilon),1);
for i = 1:length(epsilon)
    [~,~,~,k,~,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule,epsilon(i),sims,q);   
    k8(i) = mean(k);
end

% save('sol.mat')


%% Figure 5
% Plot the results from the section above. 


clear; close all; clc;

load('sol.mat')

% Plot theory vs. numerics
figure(1)
hold on
plot(k1,epsilon,'--xk','LineWidth',1,'HandleVisibility','off')
plot(k2,epsilon,'--xk','LineWidth',1,'HandleVisibility','off')
plot(k3,epsilon,'--xk','LineWidth',1,'HandleVisibility','off')
plot(k4,epsilon,'--xk','LineWidth',1,'HandleVisibility','off')

plot(k5,epsilon,'-','color',[0 0.4470 0.7410],'marker','o','LineWidth',1.2,'MarkerSize',9)
plot(k6,epsilon,'-','color',[0.6350 0.0780 0.1840],'marker','s','LineWidth',1.2,'MarkerSize',9)
plot(k7,epsilon,'--','color','m','marker',"pentagram",'LineWidth',1.2,'MarkerSize',5)
plot(k8,epsilon,'--','color','g','marker','+','LineWidth',1.2,'MarkerSize',5)

% plot(kp,epsilon,'-*m','LineWidth',1.2)

annotation('textbox',[.75 .57 .3 .3],'string','$q=0$','edgecolor','none','FitBoxToText','on','interpreter','latex')
annotation('textbox',[.75 .31 .3 .3],'string','$q=5$','edgecolor','none','FitBoxToText','on','interpreter','latex')
annotation('textbox',[.75 .04 .3 .3],'string','$q=10$','edgecolor','none','FitBoxToText','on','interpreter','latex')
annotation('textbox',[.54 .0 .3 .25],'string','$q=25$','edgecolor','none','FitBoxToText','on','interpreter','latex')


hold off
box on; grid on;

xlabel('$\textbf{E}[k]$','interpreter','latex');
xlim([2 20])
xticks((0:2:21))

ylabel('Stopping tolerance $\epsilon$','interpreter','latex');
ylim([10^(-15) 10^(0)]); 
%yticks(10.^[-25,-20,-15,-10,-5,0,5])
set(gca,'yscale','log')

legend({'SR1','SR2','SR3','SR4'},'interpreter','latex','location','southwest')
set(gca,'FontSize',12)
