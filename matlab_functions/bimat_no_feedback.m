function [ results ] = bimat_no_feedback( time,x0,pars )
%This function will integrate bimatrix resource-independent replicator
%dynamics as described in main text section 2.1.1
% Parameters
x=x0(1); %Initial fraction cooperator hosts
y=x0(2); %Initial fraction ferrojan viruses
% Establishing some combinations of payoff parameters for more succint
% dynamics functions
q1=pars.hmat(1,1)-pars.hmat(2,1);
q2=pars.hmat(1,2)-pars.hmat(2,2);
q3=pars.vmat(1,1)-pars.vmat(1,2);
q4=pars.vmat(2,1)-pars.vmat(2,2);
% Dynamics
xdot=x*(1-x)*(q1*y+q2*(1-y));
ydot=y*(1-y)*(q3*x+q4*(1-x));
% Output
results=[xdot;ydot];

end
