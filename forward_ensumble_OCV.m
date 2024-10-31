function [c] = forward_ensumble_OCV(beta, timeseq, C0)
%%% -----------------------------------------------------------------
% solve forward ODE system as a forward problem for the OCV part, this does
% not solve the nonlinear part of the model.
% Input:     c0 = initial condition for concentrations c1 and c2 = initial [c1, c2]
%            beta = parameters of the linear model. 
%            timeseq = time
% Output:    c1(t) = solid state concentration evolution in time
%            c2(t) = side reaction concentration evolution in time
%            This data is stored in c as a n*2 matrix.
% constants: rho = constant that controls the competition between
% intercalation and side reaction;

% Solve ODE system
opts = odeset('RelTol',1e-13,'AbsTol',1e-2);
[t,c] = ode45(@(t,c) odefun(t,c,beta,timeseq),timeseq,C0,opts);



function dcdt = odefun(t,c,beta,timeseq)
    % function for ODE solver
    dcdt(1) =  beta(1) + beta(2)*c(1) + 0.006*beta(3)*c(2);
    dcdt(2) =  beta(4)*c(1) - beta(3)*c(2);
    dcdt = transpose(dcdt);
end


end


