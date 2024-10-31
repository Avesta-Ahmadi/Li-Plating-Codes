function [cstar] = adjoint_ensumble_OCV(beta, cforward, c1tilde, c2tilde, timeseq)
%%% -----------------------------------------------------------------
% solve adjoint ODE system as a terminl value problem. The only trick is
% in reversing the time for the solver. 
% Input:     omega_1 and omega_2 as functions of c1 and c2 in time
%            cforward(t) as the solution of forward problem
%            Omega function as a function of c1 and c2!
%            t = time sequence of NMR data!
% Output:    cstar(t) = solid state concentration evolution in time for the adjoint system
%            This data is stored in c as a n*2 matrix.
% constants: rho = constant that controls the competition between
% intercalation and side reaction;
global weight

opts = odeset('RelTol',1e-4,'AbsTol',1e-8);
cstar0 = [0, 0];
time = flip(timeseq);

[t,cstar] = ode45(@(t,c) odefun(t,c,beta,timeseq,cforward,c1tilde,c2tilde),time,cstar0,opts);
cstar = flip(cstar);

function dcdt = odefun(t,c,beta,timeseq,cforward,c1tilde,c2tilde)
    % function for ODE solver
    c1for = interp1(timeseq, cforward(:,1), t, 'spline');    % interpolate c1(t) concentration
    c2for = interp1(timeseq, cforward(:,2), t, 'spline');    % interpolate c2(t) concentration
    c1til = interp1(timeseq, c1tilde, t, 'spline');          % interpolate experimental concentration
    c2til = interp1(timeseq, c2tilde, t, 'spline');          % interpolate experimental concentration
    dcdt(1) =         c1for - c1til  - beta(2)*c(1)       - beta(4)*c(2);
    dcdt(2) = weight*(c2for - c2til) - 0.006*beta(3)*c(1) + beta(3)*c(2);
    dcdt = transpose(dcdt);
end


end


