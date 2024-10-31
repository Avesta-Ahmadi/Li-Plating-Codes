function [cstar] = adjoint_FullModel(Omega_1, Omega_2, cforward, lambda)
%%% -----------------------------------------------------------------
% solve adjoint ODE system as a terminl value problem. The only trick is
% in reversing the time for the solver. This code solves adjoint model for
% the full model, when both linear and nonlinear parts are present in the
% ODE system. Beta parameters need to be defined for this. Also, this code
% is used for when fitting to one cycle only.
% Input:     omega_1 and omega_2 as functions of c1 and c2 in time
%            cforward(t) as the solution of forward problem
%            Omega function as a function of c1 and c2!
%            t = time sequence of NMR data!
% Output:    cstar(t) = solid state concentration evolution in time for the adjoint system
%            This data is stored in c as a n*2 matrix.
% constants: rho = constant that controls the competition between
% intercalation and side reaction;


global c1tilde c2tilde c1grid c2grid timeseq I Ngrid beta weight

% Pre-compute derivatives of omega wrt state variable at all points
h1 = (c1grid(end) - c1grid(1))/(Ngrid-1);
h2 = (c2grid(end) - c2grid(1))/(Ngrid-1);
dom1 = differentiation(Omega_1, h1);
dom2 = differentiation(Omega_2, h2);


% Solve ODE system
opts = odeset('RelTol',1e-2,'AbsTol',1e-2);
cstar0 = [0, 0];
time = flip(timeseq);

[t,cstar] = ode45(@(t,c) odefun(t,c,Omega_1,Omega_2,lambda,I,timeseq,cforward,dom1,dom2,beta),time,cstar0,opts);
cstar = flip(cstar);

function dcdt = odefun(t,c,Omega_1,Omega_2,lambda,I,timeseq,cforward,dom1,dom2,beta)
    % function for ODE solver
    i = interp1(timeseq, I, t, 'spline');                    % interpolate current
    c1for = interp1(timeseq, cforward(:,1), t, 'spline');    % interpolate c1(t) concentration
    c2for = interp1(timeseq, cforward(:,2), t, 'spline');    % interpolate c2(t) concentration
    c1til = interp1(timeseq, c1tilde, t, 'spline');          % interpolate experimental concentration
    c2til = interp1(timeseq, c2tilde, t, 'spline');          % interpolate experimental concentration
    dd1 = interp1(c1grid, dom1, c1for,'spline');             % interpolate derivative of function
    dd2 = interp1(c2grid, dom2, c2for,'spline');             % interpolate derivative of function    
    om1 = interp1(c1grid, Omega_1, c1for,'spline');          % interpolate omega at equilibrium state
    om2 = interp1(c2grid, Omega_2, c2for,'spline');          % interpolate omega at equilibrium state
    dcdt(1) =        (c1for - c1til) - beta(2)*c(1) - beta(4)*c(2)         + lambda*0.006*i*om2*dd1*c(1) - lambda*i*om2*dd1*c(2);
    dcdt(2) = weight*(c2for - c2til) + beta(3)*c(2) - 0.006*beta(3)*c(1)   + lambda*0.006*i*om1*dd2*c(1) + lambda*i*om1*dd2*c(2); 
    dcdt = transpose(dcdt);
end


end


