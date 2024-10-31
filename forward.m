function [c] = forward(Omega_1, Omega_2, lambda)
%%% -----------------------------------------------------------------
% solve forward ODE system as a forward problem only for the nonlinear
% part, assuming there is no linear part to the model
% Input:     I = applied current [A]
%            c0 = initial condition for concentrations c1 and c2 = initial [c1, c2]
%            t0 = initial time [s]
%            tf = final time [s]
%            Omega function as a function of c1 and c2!
%            t = time sequence of NMR data!
% Output:    c1(t) = solid state concentration evolution in time
%            c2(t) = side reaction concentration evolution in time
%            This data is stored in c as a n*2 matrix.
% constants: rho = constant that controls the competition between
% intercalation and side reaction;


global c1grid c2grid timeseq I rho_1 rho_2 sf C0 

% Solve ODE system
opts = odeset('RelTol',1e-13,'AbsTol',1e-2);
[t,c] = ode45(@(t,c) odefun(t,c,rho_1,rho_2,Omega_1,Omega_2,I,timeseq,sf,lambda),timeseq,C0,opts);



function dcdt = odefun(t,c,rho_1,rho_2,Omega_1,Omega_2,I,timeseq,sf,lambda)
    % function for ODE solver
    i = interp1(timeseq, I, t,'spline'); % interpolate current
    om1 = interp1(c1grid, Omega_1, c(1),'spline');
    om2 = interp1(c2grid, Omega_2, c(2),'spline');
    dcdt(1) = lambda*sf*rho_1(i)*om1*om2;
    dcdt(2) = lambda*sf*rho_2(i)*(1-om1*om2);
    dcdt = transpose(dcdt);
end


end


