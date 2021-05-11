%% Simulating No-environmental feedback Bimatrix Model
% For model details, see main text section 2.1 and appendix A

%% Adding functions
addpath('../matlab_functions')

%% First case we will discuss, iron starved conditions
% We assume cooperator hosts have a unilateral advantage
% Setting parameters
hmat=[4,3;2,1]; % host payoff matrix
vmat=[2,1;1,3]; % virus payoff matrix
model_pars.hmat=hmat; %Assigning parameters for integration
model_pars.vmat=vmat;
t=[0 20]; % Running for 20 generations
ics=[0.1,0.1]; % Initial fractions cooperator hosts, ferrojan viruses
options=odeset('reltol',1e-7,'MaxStep',1e-3); % Numerical integration setting

% Solving model
[t1_out,xy1_out]=ode45(@bimat_no_feedback,t,ics,options,model_pars);

% Solving variant of this case where alpha<beta
vmat2=[1,2;1,3];
model_pars.vmat=vmat2;
[t1b_out,xy1b_out]=ode45(@bimat_no_feedback,t,ics,options,model_pars);

%% Now iron replete conditions
% We assume cooperator hosts have a unilateral disadvantage
% Setting parameters
hmat=[1,2;3,4];
vmat=[2,1;3,1];
model_pars.hmat=hmat;
model_pars.vmat=vmat;

% Simulating model
[t2_out,xy2_out]=ode45(@bimat_no_feedback,t,ics,options,model_pars);

% Simulating for gamma<delta
vmat2=[2,1;1,3];
model_pars.vmat=vmat2;
[t2b_out,xy2b_out]=ode45(@bimat_no_feedback,t,ics,options,model_pars);

%% Plotting
f=figure;

f(1)=subplot(2,2,1);
hold on
plot(xy1_out(:,1),xy1_out(:,2),'LineWidth',1.5)
scatter(xy1_out(1,1),xy1_out(1,2),120,'ro')
scatter(xy1_out(end,1),xy1_out(end,2),120,'r','filled')
xlim([-0.1 1.1])
ylim([-0.1 1.1])
title('Iron-Starved, \alpha>\beta')
xlabel('Proportion Producer Hosts (P)')
ylabel('Proportion Ferrojan Viruses (F)')
hold off

f(2)=subplot(2,2,2);
plot(xy2_out(:,1),xy2_out(:,2),'LineWidth',1.5)
hold on
scatter(xy2_out(1,1),xy2_out(1,2),120,'ro')
scatter(xy2_out(end,1),xy2_out(end,2),120,'r','filled')
xlim([-0.1 1.1])
title('Iron-Replete, \gamma>\delta')
ylim([-0.1 1.1])
xlabel('Proportion Producer Hosts (P)')
ylabel('Proportion Ferrojan Viruses (F)')
hold off

f(3)=subplot(2,2,3);
hold on
plot(t1_out,xy1_out(:,1),'LineWidth',1.5,'color','green')
plot(t1_out,xy1_out(:,2),'LineWidth',1.5,'color','magenta')
scatter(t1_out(1),xy1_out(1,1),120,'ro')
scatter(t1_out(1),xy1_out(1,2),120,'ro')
scatter(t1_out(end),xy1_out(end,1),120,'r','filled')
scatter(t1_out(end),xy1_out(end,2),120,'r','filled')
ylim([-0.1 1.1])
xlabel('Generations')
ylabel('Prevalence of Strategy')
legend('Producer Host (P)','Ferrojan Virus (F)','Location','east')
legend boxoff
hold off

f(4)=subplot(2,2,4);
hold on
plot(t2_out,xy2_out(:,1),'LineWidth',1.5,'color','green')
plot(t2_out,xy2_out(:,2),'LineWidth',1.5,'color','magenta')
scatter(t2_out(1),xy2_out(1,1),120,'ro')
scatter(t2_out(1),xy2_out(1,2),120,'ro')
scatter(t2_out(end),xy2_out(end,1),120,'r','filled')
scatter(t2_out(end),xy2_out(end,2),120,'r','filled')
ylim([-0.1 1.1])
xlabel('Generations')
ylabel('Prevalence of Strategy')
legend('Producer Host (P)','Ferrojan Virus (F)','Location','east')
legend boxoff
hold off

set(f,'FontSize',18)
fullplot=gcf;
fullplot.PaperUnits='inches';
fullplot.PaperPosition=[0,0,12,8];
print(fullplot,'../figures/figure_2.png','-dpng','-r600')