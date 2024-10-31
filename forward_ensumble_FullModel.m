function [c] = forward_ensumble_FullModel(Omega_1, Omega_2, lambda, timeseq, I, C0)
%%% -----------------------------------------------------------------
% solve forward ODE system as a forward problem for the full model,
% assuming bot hlinear and nonlinear parts.
% Input:     lambda = lambda conversion factor
%            omega_1 and omega_2 = constitutive relations
%            beta = parameters of the linear model
% Output:    c1(t) = solid state concentration evolution in time
%            c2(t) = side reaction concentration evolution in time




global c1grid c2grid beta

% Solve ODE system
opts = odeset('RelTol',1e-13,'AbsTol',1e-2);
[t,c] = ode45(@(t,c) odefun(t,c,Omega_1,Omega_2,lambda,I,timeseq,beta),timeseq,C0,opts);



function dcdt = odefun(t,c,Omega_1,Omega_2,lambda,I,timeseq,beta)
    % function for ODE solver
    i = interp1(timeseq, I, t,'spline'); % interpolate current
    om1 = interp1(c1grid, Omega_1, c(1),'spline');
    om2 = interp1(c2grid, Omega_2, c(2),'spline');
    dcdt(1) = beta(1) + beta(2)*c(1) + 0.006*beta(3)*c(2) + lambda*(1 - 0.006*om1*om2)*i;
    dcdt(2) = beta(4)*c(1) - beta(3)*c(2)                 + lambda*(om1*om2)*i;
    dcdt = transpose(dcdt);
end


end


