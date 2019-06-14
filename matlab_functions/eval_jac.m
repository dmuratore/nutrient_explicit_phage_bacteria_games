function [ eigtab ] = eval_jac( xstar,ystar,nstar,hmat_n0,hmat_n1,vmat_n0,vmat_n1,thetax,thetay )
% This fucntion is designed to evaluate the linearization of the bimatrix
% replicator dynamical model coupled to environmental feedback about a
% given fixed point, for model details see appendix B.0.4
% Input parameters
% xstar - equilibrium cooperator host frequency
% ystar - equilibrium ferrojan virus frequency
% nstar - equilibrium environmental condition
% hmat_n0 - Host payoffs for n=0 condition
% hmat_n1 - Host payoffs for n=1 condition
% vmat_n0 - Virus payoffs for n=0 condition
% vmat_n1 - Virus payoffs for n=1 condition
% thetax - Environmental restoration rate per frequency cooperator host
% thetay - Environmental degradation rate per frequency ferrojan virus

% Setting up some helpful quantities
qs=get_pars(hmat_n0,vmat_n0);
qps=get_pars(hmat_n1,vmat_n1);
del1=qps(1)-qs(1);
del2=qps(2)-qs(2);
del3=qps(3)-qs(3);
del4=qps(4)-qs(4);

% Generating table of eigenvalues
eigtab=[];
% For each equilibrium provided
for i=1:length(xstar)
    % Calculate partial derivatives
    dxdx=(1-2*xstar(i))*((qs(1)+del1*nstar(i))*ystar(i)+(qs(2)+del2*nstar(i))*(1-ystar(i)));
    dxdy=xstar(i)*(1-xstar(i))*((qs(1)+del1*nstar(i))-(qs(2)+del2*nstar(i)));
    dxdn=xstar(i)*(1-xstar(i))*(del1*ystar(i)+del2*(1-ystar(i)));
    dydx=ystar(i)*(1-ystar(i))*((qs(3)+del3*nstar(i))-(qs(4)+del4*nstar(i)));
    dydy=(1-2*ystar(i))*((qs(3)+del3*nstar(i))*xstar(i)+(qs(4)-del4*nstar(i))*(1-xstar(i)));
    dydn=ystar(i)*(1-ystar(i))*(del3*xstar(i)+del4*(1-xstar(i)));
    dndx=nstar(i)*(1-nstar(i))*thetax;
    dndy=-nstar(i)*(1-nstar(i))*thetay;
    dndn=(1-2*nstar(i))*(-1+thetax*xstar(i)-thetay*ystar(i));
    % Set up Jacobian
    J=[[dxdx;dxdy;dxdn],[dydx;dydy;dydn],[dndx;dndy;dndn]];
    % Evaluate Eigenvalues of Jacobian
    Jeigs=eig(J);
    % Tack on to output
    eigtab(i,:)=Jeigs;
end


end

