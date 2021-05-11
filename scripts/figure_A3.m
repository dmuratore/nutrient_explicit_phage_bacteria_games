%%% Generating Figure Supplemental for Heteroclinic Network Stability
%% Establishing path
addpath('../matlab_functions/')

%% Executing Code
% For details about interpretation, see appendix B.0.4
% Establishing parameters
hmat_n0=[2,2;1,1]; %Host payoff matrix for n=0 conditions
vmat_n0=[4,3;4,3]; %Virus payoff matrix for n=0 conditions
hmat_n1=[3,2;4,4.1]; %Host payoff matrix for n=1 conditions
vmat_n1=[1,3;1,3]; %Virus payoff matrix for n=1 conditions
thetay=1; %Environmental depletion for ferrojan virus frequency
thetax=4; %Environmental restoration for cooperator host frequency
% Model attributes
% Evaluate internal n fixed points
hc_nstar=solve_nstar(hmat_n0,hmat_n1,vmat_n0,vmat_n1,thetax,thetay);
% Evaluate x+y fixed points corresponding to the above values of n
hc_xystar=solve_xystar(hc_nstar,hmat_n0,hmat_n1,thetax,thetay);

% Solving the same things to show interior fixed points are not feasible
hc2_nstar=solve_nstar(hmat_n0,hmat_n1+[0,0;0,-0.5],vmat_n0,vmat_n1,thetax,thetay);
hc2_xystar=solve_xystar(hc2_nstar,hmat_n0,hmat_n1+[0,0;0,-0.5],thetax,thetay);

% Generating random trajectories for interior fixed point feasible/not
% feasible conditions
rng(517803705)
for i=1:10
hc_net=run_feedback(hmat_n0,vmat_n0,hmat_n1,vmat_n1,thetax,thetay,true);
hc_2=run_feedback(hmat_n0,vmat_n0,hmat_n1+[0,0;0,-0.5],vmat_n1,thetax,thetay,true);
hc_x{i}=hc_net.outs;
nhc_x{i}=hc_2.outs;
end

% Generating figure
hc_comp_fig=[];
figure;
hc_comp_fig(1)=subplot(1,2,1);
hold on
quiver3(hc_net.X,hc_net.Y,hc_net.N,hc_net.dxdt,hc_net.dydt,hc_net.dndt)
for i=1:10
plot3(hc_x{1,i}(:,1),hc_x{1,i}(:,2),hc_x{1,i}(:,3),'LineWidth',1.5)
scatter3(hc_x{1,i}(1,1),hc_x{1,i}(1,2),hc_x{1,i}(1,3),120,'o')
scatter3(hc_x{1,i}(end,1),hc_x{1,i}(end,2),hc_x{1,i}(end,3),120,'filled')
end
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('n','FontSize',18)
axis([-0.1 1.1,-0.1 1.1,-0.1 1.1])
view([-79.9 15.6]);
hold off
hc_comp_fig(2)=subplot(1,2,2);
hold on
quiver3(hc_2.X,hc_2.Y,hc_2.N,hc_2.dxdt,hc_2.dydt,hc_2.dndt)
for i=1:10
plot3(nhc_x{1,i}(:,1),nhc_x{1,i}(:,2),nhc_x{1,i}(:,3),'LineWidth',1.5)
scatter3(nhc_x{1,i}(1,1),nhc_x{1,i}(1,2),nhc_x{1,i}(1,3),120,'o')
scatter3(nhc_x{1,i}(end,1),nhc_x{1,i}(end,2),nhc_x{1,i}(end,3),120,'filled')
end
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('n','FontSize',18)
axis([-0.1 1.1,-0.1 1.1,-0.1 1.1])
view([-67.8 13])
hold off
print_fig=gcf;
%% Figure
% Saving figure to output
print(print_fig,'../figures/figure_A3.png','-dpng','-r600')