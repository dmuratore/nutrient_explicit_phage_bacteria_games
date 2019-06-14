function results = env_feedback( time,ics,pars )
% The purpose of this function is to run the environmentally-couple
% bimatrix replicator dynamical model as described in main text section 2.2.1
% The time
% variable is for numerical integration using ode45.
% Parameters
x=ics(1); %Fraction of hosts cooperating
y=ics(2); %Fraction of viruses ferrojan
n=ics(3); %Environmental condition
% Get these parameters from the get_pars.m function
% They are combinations of the payoff matrices
% hmat_n0,hmat_n1,vmat_n0,vmat_n1
q1=pars.q1;
q1p=pars.q1p;
q2=pars.q2;
q2p=pars.q2p;
q3=pars.q3;
q3p=pars.q3p;
q4=pars.q4;
q4p=pars.q4p;
% Environmental restoration parameters
thetax=pars.thetax;
thetay=pars.thetay;
% Dynamical System
dxdt=x*(1-x)*((q1+(q1p-q1)*n)*y+(q2+(q2p-q2)*n)*(1-y));
dydt=y*(1-y)*((q3+(q3p-q3)*n)*x+(q4+(q4p-q4)*n)*(1-x));
dndt=n*(1-n)*(-1+thetax*x-thetay*y);
% Outputs
results=[dxdt;dydt;dndt];

end

