%% Running Bimatrix Model Coupled to Environmental Feedback
function results=run_feedback(mat1,mat2,mat3,mat4,tx,ty,icrand)
%% Establishing Parameters
% Describe payoff matrices for hosts and viruses in n=0 case
hmat_n0=mat1;
vmat_n0=mat2;
% Describe payoff matrices for hosts and viruses in n=1 case
hmat_n1=mat3;
vmat_n1=mat4;
% Getting model parameters from these matrices (for brevity in
% implementation)
% These parameters are just the differences between payoffs for either host
% or virus strategy a la a-c, b-d, alpha-beta, gamma-delta
n0_qs=get_pars(hmat_n0,vmat_n0);
n1_qs=get_pars(hmat_n1,vmat_n1);
% Setting environmental restoration/degradation parameters
thetax=tx;
thetay=ty;
% Establishing parameters as parameter object
pars.q1=n0_qs(1);
pars.q2=n0_qs(2);
pars.q3=n0_qs(3);
pars.q4=n0_qs(4);
pars.q1p=n1_qs(1);
pars.q2p=n1_qs(2);
pars.q3p=n1_qs(3);
pars.q4p=n1_qs(4);
pars.thetax=thetax;
pars.thetay=thetay;
%% Setting up integration functional parameters
t=[0,200]; %Number of generations to run
% Setting up grid for phase diagrams
x=0:0.1:1;
y=0:0.1:1;
n=0:0.1:1;
[X,Y,N]=meshgrid(x,y,n);
%% Filling out arrows for phase diagram
dxdt=X.*(1-X).*((pars.q1+(pars.q1p-pars.q1).*N).*Y+(pars.q2+(pars.q2p-pars.q2).*N).*(1-Y));
dydt=Y.*(1-Y).*((pars.q3+(pars.q3p-pars.q3).*N).*X+(pars.q4+(pars.q4p-pars.q4).*N).*(1-X));
dndt=N.*(1-N).*(-1+pars.thetax.*X-pars.thetay.*Y);
%% Evaluating integration of model
% Random initial conditions?
if icrand==true
    ics=[rand rand rand];
else
    ics=[0.1 0.1 0.1];
end
% Adjusting Runge-Kutta 4th order integration parameters
options=odeset('reltol',1e-7,'MaxStep',1e-3);
% Evaluating numerical integration
[ts,outs]=ode45(@env_feedback,t,ics,options,pars);
%% Saving outputs
results.X=X;
results.Y=Y;
results.N=N;
results.dxdt=dxdt;
results.dydt=dydt;
results.dndt=dndt;
results.ts=ts;
results.outs=outs;
end