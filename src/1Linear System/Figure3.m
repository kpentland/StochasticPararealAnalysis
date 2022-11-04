%% Figure 3

%% Second moments vs. k 
% Here we plot the second moments of the state-dependent perturbations in
% SParareal. 

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
sample_rule1 = 1;                            %sampling rule 1
sample_rule2 = 2;                            %sampling rule 2
sample_rule3 = 3;                            %sampling rule 3
sample_rule4 = 4;                            %sampling rule 4
sims = 500;                                  %no. of independent sims
q = 0;

%solve with stochastic parareal
[~,U_para1,~,k1,UG1,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule1,epsilon,sims,q);
[~,U_para2,~,k2,UG2,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule2,epsilon,sims,q);
[~,U_para3,~,k3,UG3,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule3,epsilon,sims,q);
[~,U_para4,~,k4,UG4,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule4,epsilon,sims,q);


% extract the maximal (over the dimension of the system) seond moment at each
% time and iteration (and simulation)
tempUG1 = zeros(sims,1);
tempUG2 = zeros(sims,1);
tempUG3 = zeros(sims,1);
tempUG4 = zeros(sims,1);
for j = 1:sims
    temp1 = UG1{j,1};
    temp2 = UG2{j,1};
    temp3 = UG3{j,1};
    temp4 = UG4{j,1};
    for i = 1:N
        ind = (n*(i-1)+1:n*i);
        ind_next = ((n*i)+1:n*(i+1));
        tempUG1(j,i) = max( ( vecnorm(temp1(:,ind_next) - temp1(:,ind),inf,2) ).^2 );
        tempUG2(j,i) = max( ( vecnorm(temp2(:,ind_next) - temp2(:,ind),inf,2) ).^2 );
        tempUG3(j,i) = max( ( vecnorm(temp3(:,ind_next) - temp3(:,ind),inf,2) ).^2 );
        tempUG4(j,i) = max( ( vecnorm(temp4(:,ind_next) - temp4(:,ind),inf,2) ).^2 );
    end
end
% take the mean over all simulations (at each time)
temppp1 = mean(tempUG1,1); temppp1(temppp1 == 0) = 10^(-100);
temppp2 = mean(tempUG2,1); temppp2(temppp2 == 0) = 10^(-100);
temppp3 = mean(tempUG3,1); temppp3(temppp3 == 0) = 10^(-100);
temppp4 = mean(tempUG4,1); temppp4(temppp4 == 0) = 10^(-100);



% plot the maximal moments at each iteration
figure(1)
DT = (tspan(2)-tspan(1))/N;                      %time slice width
hold on
plot((0:N),(DT^((2*0) + 1))*ones(1,N+1),'--k','LineWidth',1,'HandleVisibility','off')
plot((0:N),(DT^((2*5) + 1))*ones(1,N+1),'--k','LineWidth',1,'HandleVisibility','off')
plot((0:N),(DT^((2*10) + 1))*ones(1,N+1),'--k','LineWidth',1,'HandleVisibility','off')

plot((1:N),temppp1,'-','color',[0 0.4470 0.7410],'marker','o','LineWidth',1.2,'MarkerSize',9)
plot((1:N),temppp2,'-','color',[0.6350 0.0780 0.1840],'marker','s','LineWidth',1.2,'MarkerSize',9)
plot((1:N),temppp3,'--','color','m','marker',"pentagram",'LineWidth',1.2,'MarkerSize',5)
plot((1:N),temppp4,'--','color','g','marker','+','LineWidth',1.2,'MarkerSize',5)


annotation('textbox',[.78 .58 .3 .3],'string','$q=0$','edgecolor','none','FitBoxToText','on','interpreter','latex')
annotation('textbox',[.78 .4 .3 .3],'string','$q=5$','edgecolor','none','FitBoxToText','on','interpreter','latex')
annotation('textbox',[.78 .22 .3 .3],'string','$q=10$','edgecolor','none','FitBoxToText','on','interpreter','latex')

hold off
xlabel('$k$','interpreter','latex');
ylabel('Largest second moment of $\xi^k_n(U^k_n)$','interpreter','latex');
box on; grid on;
xlim([0 20])
ylim([10^(-40) 10^(5)]);
xticks([0,5,10,15,20])
% yticks(10.^[-15,-10,-5,0])
set(gca,'yscale','log')
legend({'SR1','SR2','SR3','SR4'},'interpreter','latex','location','southwest')
set(gca,'FontSize',12)
