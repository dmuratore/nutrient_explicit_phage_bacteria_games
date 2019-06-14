function [ xy_eqs ] = solve_xystar( neqs,hmat_n0,hmat_n1,thetax,thetay )
% This function is designed to solve the x+y equilibria of the
% environmental feedback coupled model given solutions for n (can be found
% using solve_nstar.m)
% Parameters: neqs n equilbrium values (found with solve_nstar.m)
% hmat_n0 Host payoff matrix for n=0 conditions
% hmat_n1 Host payoff matrix for n=1 conditions
% thetax environmental restoration rate for cooperators x
% thetay environmental destruction rate for ferrojan viruses y

% Initializing some parameters for easier typing
q2=hmat_n0(1,2)-hmat_n0(2,2);
q2p=hmat_n1(1,2)-hmat_n1(2,2);
q1=hmat_n0(1,1)-hmat_n0(2,1);
q1p=hmat_n1(1,1)-hmat_n1(2,1);
del2=q2p-q2;
del1=q1p-q1;

% Solving x and y nullclines in terms of nstar
yeqs=(q2+del2.*neqs)./((q2+del2.*neqs)-(q1+del1.*neqs));
xeqs=(1+thetay.*yeqs)./thetax;

% Saving to output
xy_eqs=[xeqs,yeqs];

end

