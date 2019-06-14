%% Function to generate parameters from payoff matrices
function results=get_pars(hmat,vmat)
% For details of interpretation of model parameters, see main text section
% 2.1.1
% Initializing parameters
% hmat - host payoff matrix
% vmat - virus payoff matrix
a=hmat(1,1);
b=hmat(1,2);
c=hmat(2,1);
d=hmat(2,2);
alpha=vmat(1,1);
beta=vmat(1,2);
gamma=vmat(2,1);
delta=vmat(2,2);
% Defining linear combinations of payoff convenient for later computations
q1=a-c;
q2=b-d;
q3=alpha-beta;
q4=gamma-delta;
% Output
results=[q1,q2,q3,q4];


end