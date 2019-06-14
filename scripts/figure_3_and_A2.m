%% Simulating Environmental Feedback Bimatrix Model
% For model details, see main text section 2.2 and appendix B

%% Adding filepath
addpath('../matlab_functions/')

%% We also set the exterior equilibria for stability analysis
exterior_x=[1,1,1,1,0,0,0,0]';
exterior_y=[0,0,1,1,0,0,1,1]';
exterior_n=[0,1,0,1,0,1,0,1]';

%% First case we want to discuss, trivial disaster case
% see [main text results] for details
% Setting parameters
hmat_n0=[3,2;5,4]'; %Host payoffs for n=0 case
hmat_n1=[3,2;5,4]'; %Host payoffs for n=1 case
vmat_n0=[2.5,2;4.5,4]; %Virus payoffs for n=0 case
vmat_n1=[2.5,2;4.5,4]; %Virus payoffs for n=1 case
thetay=1; %Environmental depletion per frequency ferrojan virus
thetax=1; %Environmental restoration per frequency cooperator host

%% Find internal fixed points
rc_nstar=solve_nstar(hmat_n0,hmat_n1,vmat_n0,vmat_n1,thetax,thetay);
rc_xystar=solve_xystar(rc_nstar,hmat_n0,hmat_n1,thetax,thetay);

% Perform stability analysis on all fixed points
rc_fpx=vertcat(rc_xystar(:,1),exterior_x);
rc_fpy=vertcat(rc_xystar(:,2),exterior_y);
rc_fpn=vertcat(rc_nstar,exterior_n);
rc_stability=eval_jac(rc_fpx,rc_fpy,rc_fpn,hmat_n0,hmat_n1,vmat_n0,vmat_n1,thetax,thetay);

% Numerical simulation of dynamics
resource_crash=run_feedback(hmat_n0,vmat_n0,hmat_n1,vmat_n1,thetax,thetay,false);
%% Next case, heteroclinic network

% Setting parameters as above
hmat_n0=[4,4;3,2];
vmat_n0=[3,4;1,4];
hmat_n1=[1,1;2,2.1];
vmat_n1=[2,1;3,1];
thetay=1;
thetax=4;

%% Solving interior fixed points
hc_nstar=solve_nstar(hmat_n0,hmat_n1,vmat_n0,vmat_n1,thetax,thetay);
hc_xystar=solve_xystar(hc_nstar,hmat_n0,hmat_n1,thetax,thetay);
hc_fpx=vertcat(hc_xystar(:,1),exterior_x);
hc_fpy=vertcat(hc_xystar(:,2),exterior_y);
hc_fpn=vertcat(hc_nstar,exterior_n);

% Evaluating stability interior fixed points
hc_stability=eval_jac(hc_fpx,hc_fpy,hc_fpn,hmat_n0,hmat_n1,vmat_n0,vmat_n1,thetax,thetay);

% Simulate dynamics
hc_net=run_feedback(hmat_n0,vmat_n0,hmat_n1,vmat_n1,thetax,thetay,false);

%% Next case, neutral orbits

% Let's say in this case viral iron incorporation is dominating
hmat_n0=[4,4;3,2];
hmat_n1=[1,1;2,2];
vmat_n0=[3,2;2,1];
vmat_n1=[4,2;3,1];
thetax=2.5;
thetay=1;

%% Solve interior fixed points
o2d_nstar=solve_nstar(hmat_n0,hmat_n1,vmat_n0,vmat_n1,thetax,thetay);
o2d_xystar=solve_xystar(o2d_nstar,hmat_n0,hmat_n1,thetax,thetay);
o2d_fpx=vertcat(o2d_xystar(:,1),exterior_x);
o2d_fpy=vertcat(o2d_xystar(:,2),exterior_y);
o2d_fpn=vertcat(o2d_nstar,exterior_n);

% Evaluate stability all fixed points
o2d_stability=eval_jac(o2d_fpx,o2d_fpy,o2d_fpn,hmat_n0,hmat_n1,vmat_n0,vmat_n1,thetax,thetay);

% Run dynamics
orbit_2d=run_feedback(hmat_n0,vmat_n0,hmat_n1,vmat_n1,thetax,thetay,false);
%% Last Case, 3d neutral orbit

% Setting up parameters
vmat_n0=[3,2;3,2];
vmat_n1=[2,3;2,3];
hmat_n0=[4,4;3.5,3.5];
hmat_n1=[1.5,1.5;2,2];
thetax=2.5;
thetay=0.5;

%% Finding interior fixed points
o3d_nstar=solve_nstar(hmat_n0,hmat_n1,vmat_n0,vmat_n1,thetax,thetay);
o3d_xystar=solve_xystar(o3d_nstar,hmat_n0,hmat_n1,thetax,thetay);
o3d_fpx=vertcat(exterior_x);
o3d_fpy=vertcat(exterior_y);
o3d_fpn=vertcat(exterior_n);

% Evaluating stability fixed points
o3d_stability=eval_jac(o3d_fpx,o3d_fpy,o3d_fpn,hmat_n0,hmat_n1,vmat_n0,vmat_n1,thetax,thetay);

% Simulating dynamics
orbit_3d=run_feedback(hmat_n0,vmat_n0,hmat_n1,vmat_n1,thetax,thetay,false);

%% GENERATING FIGURE [XX]
% Making phase plane demonstration for supplemental figure X [supplemental]

% Setting up grid for phase plane
int_ystars=0:0.1:1;
int_xstars=(thetay.*int_ystars+1)./thetax;
int_nstars=ones(size(int_ystars)).*0.5;
int_stab=eval_jac(int_xstars,int_ystars,int_nstars,hmat_n0,hmat_n1,vmat_n0,vmat_n1,thetax,thetay);
plane_y=0:0.01:1;
plane_n=0:0.01:1;
[bigy,bign]=meshgrid(plane_y,plane_n);
plane_x=(1+thetay.*bigy)./thetax;
new_bign=ones(101).*0.5;

% Graphics
figure;
hold on
phase=quiver3(orbit_3d.N,orbit_3d.Y,orbit_3d.X,orbit_3d.dndt,orbit_3d.dydt,orbit_3d.dxdt);
set(phase,'AutoScale','on','AutoScaleFactor',2)
plot3(orbit_3d.outs(:,3),orbit_3d.outs(:,2),orbit_3d.outs(:,1),'LineWidth',1.5)
xlabel('n','FontSize',18)
ylabel('y','FontSize',18)
zlabel('x','FontSize',18)
axis([-0.1 1.1,-0.1 1.1,-0.1 1.1])
view([-12.7 0.5])
surf(plane_x,bigy,bign)
surf(bigy,bign,new_bign)
hold off
phase_fig=gcf;
print(phase_fig,'../figures/figure_A2.png','-dpng','-r600')

%% Setting up plotting for main text figure 3
fig=[];
figure;

fig(1)=subplot(4,2,1);
hold on
%quiver3(resource_crash.X,resource_crash.Y,resource_crash.N,resource_crash.dxdt,resource_crash.dydt,resource_crash.dndt)
plot3(resource_crash.outs(:,1),resource_crash.outs(:,2),resource_crash.outs(:,3),'LineWidth',1.5)
scatter3(resource_crash.outs(1,1),resource_crash.outs(1,2),resource_crash.outs(1,3),120,'ro')
scatter3(resource_crash.outs(end,1),resource_crash.outs(end,2),resource_crash.outs(end,3),120,'r','filled')
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('n','FontSize',18)
title('Resource Crash','FontSize',18)
axis([-0.1 1.1,-0.1 1.1,-0.1 1.1])
view(3);
hold off

fig(4)=subplot(4,2,7);
hold on
%quiver3(hc_net.X,hc_net.Y,hc_net.N,hc_net.dxdt,hc_net.dydt,hc_net.dndt)
plot3(hc_net.outs(:,1),hc_net.outs(:,2),hc_net.outs(:,3),'LineWidth',1.5)
scatter3(hc_net.outs(1,1),hc_net.outs(1,2),hc_net.outs(1,3),120,'ro')
scatter3(hc_net.outs(end,1),hc_net.outs(end,2),hc_net.outs(end,3),120,'r','filled')
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('n','FontSize',18)
title('Phenotypic Jumping','FontSize',18)
axis([-0.1 1.1,-0.1 1.1,-0.1 1.1])
view(3);
hold off

fig(2)=subplot(4,2,3);
hold on
%quiver3(orbit_2d.X,orbit_2d.Y,orbit_2d.N,orbit_2d.dxdt,orbit_2d.dydt,orbit_2d.dndt)
plot3(orbit_2d.outs(:,1),orbit_2d.outs(:,2),orbit_2d.outs(:,3),'LineWidth',1.5)
scatter3(orbit_2d.outs(1,1),orbit_2d.outs(1,2),orbit_2d.outs(1,3),120,'ro')
scatter3(orbit_2d.outs(end,1),orbit_2d.outs(end,2),orbit_2d.outs(end,3),120,'r','filled')
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('n','FontSize',18)
title('Dominating Phenotypes','FontSize',18)
axis([-0.1 1.1,-0.1 1.1,-0.1 1.1])
view(3);
hold off

fig(3)=subplot(4,2,5);
hold on
%quiver3(orbit_3d.X,orbit_3d.Y,orbit_3d.N,orbit_3d.dxdt,orbit_3d.dydt,orbit_3d.dndt)
plot3(orbit_3d.outs(:,1),orbit_3d.outs(:,2),orbit_3d.outs(:,3),'LineWidth',1.5)
scatter3(orbit_3d.outs(1,1),orbit_3d.outs(1,2),orbit_3d.outs(1,3),120,'ro')
scatter3(orbit_3d.outs(end,1),orbit_3d.outs(end,2),orbit_3d.outs(end,3),120,'r','filled')
%scatter3(o3d_xystar(1,1),o3d_xystar(1,2),o3d_nstar(1))
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('n','FontSize',18)
title('Phenotypic Oscillation','FontSize',18)
axis([-0.1 1.1,-0.1 1.1,-0.1 1.1])
view(3);
hold off

fig(5)=subplot(4,2,2);
plot(resource_crash.ts,resource_crash.outs(:,1),'g','LineWidth',1.5)
hold on
plot(resource_crash.ts,resource_crash.outs(:,2),'m','LineWidth',1.5)
plot(resource_crash.ts,resource_crash.outs(:,3),'b','LineWidth',1.5)
ylim([-0.1 1.1])
xlim([0 50])
legend('Hosts (x)','Viruses (y)','Fe State (n)','Location','east')
legend boxoff
hold off

fig(8)=subplot(4,2,8);
plot(hc_net.ts,hc_net.outs(:,1),'g','LineWidth',1.5)
hold on
plot(hc_net.ts,hc_net.outs(:,2),'m','LineWidth',1.5)
plot(hc_net.ts,hc_net.outs(:,3),'b','LineWidth',1.5)
ylim([-0.1 1.1])
xlim([0 50])
xlabel('Generations','FontSize',18)
ylabel('Prevalence of Strategy','FontSize',18)
hold off

fig(6)=subplot(4,2,4);
plot(orbit_2d.ts,orbit_2d.outs(:,1),'g','LineWidth',1.5)
hold on
plot(orbit_2d.ts,orbit_2d.outs(:,2),'m','LineWidth',1.5)
plot(orbit_2d.ts,orbit_2d.outs(:,3),'b','LineWidth',1.5)
ylim([-0.1 1.1])
xlim([0 50])
hold off

fig(7)=subplot(4,2,6);
plot(orbit_3d.ts,orbit_3d.outs(:,1),'g','LineWidth',1.5)
hold on
plot(orbit_3d.ts,orbit_3d.outs(:,2),'m','LineWidth',1.5)
plot(orbit_3d.ts,orbit_3d.outs(:,3),'b','LineWidth',1.5)
ylim([-0.1 1.1])
xlim([0 50])
hold off

set(fig,'FontSize',16)
full_figure=gcf;
full_figure.PaperUnits='inches';
full_figure.PaperPosition=[0,0,8,14];

%% Saving figure to output
print(full_figure,'../figures/figure_3.png','-dpng','-r600')