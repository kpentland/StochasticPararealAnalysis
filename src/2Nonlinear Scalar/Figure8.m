%% Figure 3

%% Second moments vs. k 
% Here we plot the second moments of the state-dependent perturbations in
% SParareal. 

clear; close all; clc

%Inputs:
f = @(t,u)( sqrt(u^2 + 2) );            %function handle for ODE
tspan = [-1,1];                         %time interval
u0 = 5;                                 %intial conditions
N = 20;                                 %no. of time sub-intervals step
Ng = N;                                 %no. of coarse steps (in each sub-interval)
Nf = Ng;                                %no. of fine steps (in each sub-interval)
epsilon = 10^(-10000);                  %error tolerance
G = 'RK1';                              %coarse solver
F = 'ExactLinear';                      %fine solver
sample_rule1 = 1;                       %sampling rule to employ
sample_rule2 = 2;                       %sampling rule to employ
sample_rule3 = 3;                       %sampling rule to employ
sample_rule4 = 4;                       %sampling rule to employ
sims = 500;                             %no. of independent sims
q = 0;

%solve with stochastic parareal
[~,U_para1,~,k1,UG1,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule1,epsilon,sims,q);
[~,U_para2,~,k2,UG2,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule2,epsilon,sims,q);
[~,U_para3,~,k3,UG3,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule3,epsilon,sims,q);
[~,U_para4,~,k4,UG4,~] = SParareal(f,tspan,u0,N,Ng,Nf,F,G,sample_rule4,epsilon,sims,q);


% extract the maximal seond moment at each time and iteration (and simulation)
tempUG1 = cell(sims,1);
tempUG2 = cell(sims,1);
tempUG3 = cell(sims,1);
tempUG4 = cell(sims,1);
for j = 1:sims
    temp1 = UG1{j,1};
    temp2 = UG2{j,1};
    temp3 = UG3{j,1};
    temp4 = UG4{j,1};
    tempUG1{j,1} = (temp1(:,2:end) - temp1(:,1:end-1) ).^2;
    tempUG2{j,1} = (temp2(:,2:end) - temp2(:,1:end-1) ).^2;
    tempUG3{j,1} = (temp3(:,2:end) - temp3(:,1:end-1) ).^2;
    tempUG4{j,1} = (temp4(:,2:end) - temp4(:,1:end-1) ).^2;
end
% take the mean over all simulations (at each time)
temppp1 = mean(cat(3,tempUG1{:}),3); temppp1(temppp1 == 0) = 10^(-100);
temppp2 = mean(cat(3,tempUG2{:}),3); temppp2(temppp2 == 0) = 10^(-100);
temppp3 = mean(cat(3,tempUG3{:}),3); temppp3(temppp3 == 0) = 10^(-100);
temppp4 = mean(cat(3,tempUG4{:}),3); temppp4(temppp4 == 0) = 10^(-100);



% plot the maximal moments at each iteration
figure(1)
DT = (tspan(2)-tspan(1))/N;                      %time slice width
hold on
plot((0:N),(DT^((2*0) + 1))*ones(1,N+1),'--k','LineWidth',1,'HandleVisibility','off')
plot((0:N),(DT^((2*5) + 1))*ones(1,N+1),'--k','LineWidth',1,'HandleVisibility','off')
plot((0:N),(DT^((2*10) + 1))*ones(1,N+1),'--k','LineWidth',1,'HandleVisibility','off')

plot((1:N),max(temppp1),'-','color',[0 0.4470 0.7410],'marker','o','LineWidth',1.2,'MarkerSize',9)
plot((1:N),max(temppp2),'-','color',[0.6350 0.0780 0.1840],'marker','s','LineWidth',1.2,'MarkerSize',9)
plot((1:N),max(temppp3),'--','color','m','marker',"pentagram",'LineWidth',1.2,'MarkerSize',5)
plot((1:N),max(temppp4),'--','color','g','marker','+','LineWidth',1.2,'MarkerSize',5)


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
