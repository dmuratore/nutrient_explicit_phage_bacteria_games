function [ neqs ] = solve_nstar( hmat_n0,hmat_n1,vmat_n0,vmat_n1,thetax,thetay )
% This function is designed to solve for equilibrium values of n (nstar)
% for any given parameterization of bimatrix replicator dynamics model with
% coupled environmental feedback. Solutions for n-star can then be plugged
% into solve_xystar.m to get x and y equilibrium coordinations
% Parameters:
% hmat_n0 Host payoff matrix for n=0 condition
% hmat_n1 Host payoff matrix for n=1 condition
% vmat_n0 Virus payoff matrix for n=0 condition
% vmat_n1 Virus payoff matrix for n=1 condition
% thetax Environmental restoration rate per cooperator host frequency
% thetay Environmental depletion rate per ferrojan virus frequency

% Defining important parameters for evaluation of nstars
q1=hmat_n0(1,1)-hmat_n0(2,1);
q2=hmat_n0(1,2)-hmat_n0(2,2);
q3=vmat_n0(1,1)-vmat_n0(1,2);
q4=vmat_n0(2,1)-vmat_n0(2,2);
del1=(hmat_n1(1,1)-hmat_n1(2,1))-q1;
del2=(hmat_n1(1,2)-hmat_n1(2,2))-q2;
del3=(vmat_n1(1,1)-vmat_n1(1,2))-q3;
del4=(vmat_n1(2,1)-vmat_n1(2,2))-q4;

% Defining some multiplicative factors based on thetas
f1=1+thetay;
f2=-1+thetax+thetay;
f3=-1+thetax;

% Solving for parameters of polynomial for nstar
a=((f1*del2*del3)+(f2*del2*del4)-(del1*del3)-(f3*del4*del1));
b=((f1*((q2*del3)+(q3*del2)))+(f2*((q2*del4)+(q4*del2)))-((q3*del1)+(q1*del3))-(f3*((q4*del1)+(del4*q1))));
c=((f1*q2*q3)+(f2*q2*q4)-(q3*q1)-(f3*q4*q1));

% Output
neqs=roots([a b c]);



end

