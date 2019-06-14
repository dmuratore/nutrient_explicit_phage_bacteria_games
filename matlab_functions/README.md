# matlab_functions:
Contains all functions necessary to run Matlab codes, files are as follows:

## eval_jac.m 
Evaluates eigenvalues of Jacobian matrix for fixed point in bimatrix replicator dynamical model coupled to environmental feedback
*INPUTS: Fixed point (x,y,n); model parameters (host payoffs n=0, host payoffs n=1, virus payoffs n=0, virus payoffs n=1, theta_x, theta_y)
*OUTPUTS: Matrix with fixed points (x,y,n) and associated eigenvalues

## get_pars.m
Reformats host and virus payoff parameters for more succinct use in model simulation
*INPUTS: Host payoff matrix, virus payoff matrix
*OUTPUTS: list of quantities (qs) referring to payoff combinations a-c,b-d,alpha-beta,gamma-delta (for details see main text section 2.1.1)


## bimat_no_feedback.m 
Solves dynamical system for bimatrix game with no environmental feedback
*INPUTS: Host and virus payoff matrices, initial frequency cooperator hosts and ferrojan viruses
*OUTPUTS: dx/dt, dy/dt


## env_feedback.m 
Solves dynamical system for bimatrix game coupled to environmental feedback
*INPUTS: Host and virus payoff parameters (output of get_pars.m), theta_x and theta_y parameters in MATLAB structure object, initial frequency cooperator hosts, ferrojan viruses, environmental state
*OUTPUTS: dx/dt, dy/dt, dn/dt


## run_feedback.m 
Wraps several functions involved in numerical simulation of environmentally-coupled bimatrix model. Performs numerical integration and generates 3D phase portrait. Can run on a random initial conditions mode or uses default initial conditions 0.1,0.1,0.1. Model set to run for 200 generations.
*INPUTS: Host payoff matrix for n=0, virus payoff matrix for n=0, host payoff matrix for n=1, virus payoff matrix for n=1, theta_x, theta_y parameters, then icrand parameter - 'true' for random initial conditions, 'false' for 0.1,0.1,0.1 initial conditions.
*OUTPUTS: MATLAB structure with 3D grid for phase portrait plotting (X,Y,N); flow corresponding to the 3D grid (dxdt,dydt,dndt), times for temporal dynamics (ts), and integrated values for x,y,n as an x by 3 matrix (outs). 


## solve_nstar.m 
Solves for interior fixed point n coordinate using quadratic equation identified in Appendix section B.0.3. 
*INPUTS: Host and virus payoff matrices for both n=0 and n=1 conditions, theta_x and theta_y parameters
*OUTPUTS: List of n-coordinates of fixed points for model


## solve_xystar.m 
Solves for x and y coordinates of the fixed points corresponding to some particular fixed point value of n. 
*INPUTS: n-coordinates (from solve_nstar.m), host payoff matrices for n=0 and n=1 conditions, theta_x and theta_y parameters
*OUTPUTS: List of x+y-coordinates of fixed points for model