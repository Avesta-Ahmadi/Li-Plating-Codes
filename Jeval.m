function [jj] = Jeval(tau, omega_1, omega_2, lambda, DD, which)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function evaluates the cost functional given a set of constitutive relations 
% for the RHS of the ODE system. The result of evaluation is saved and the evaluation 
% counter is incremented. 
% note that the problem is formulated in a decoupled manner. Meanning that
% the walk in the steepest descent direction is split into two parts, once
% for omega_1 and once for omega_2. So, the steepest descent direction will
% be zero for one of the constitutive relations, and nonzero for the other
% one --> then omega functions are updated (essentially for one of them)
% and cost function is computed.
% INPUT:
%       tau        - step length of the descent direction for one of
%       constitutive relations
%       omega_1    - the constitutive relation omega_1 at current iteration
%       (reference function)
%       omega_2    - the constitutive relation omega_2 at current iteration
%       (reference function)
%       DD        - Descent Direction for omega (either omega_1 or omega_2)
%       which     - An integer (1 or 2) specifying which constitutive relation is
%       getting updated and which one is remaining constant! which = 1
%       represents an update on omega_1 and which = 2 represents an update
%       for omega_2. 
% OUTPUT:
%       jj         -  cost function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global timeseq c1tilde c2tilde Feval LineMin1 LineMin2

if which == 1
    omega_1_n = omega_1 + tau * DD;
    omega_2_n = omega_2;
    lambda_n = lambda;
elseif which == 2
    omega_1_n = omega_1;
    omega_2_n = omega_2 + tau * DD;
    lambda_n = lambda;
else
    omega_1_n = omega_1;
    omega_2_n = omega_2;
    lambda_n = lambda + tau * DD;
end

c = forward(omega_1_n, omega_2_n, lambda_n);
jj = J(c,c1tilde,c2tilde,timeseq);

Feval = Feval + 1;

if which == 1
    LineMin1 = [ LineMin1; [ tau, jj ] ];
else
    LineMin2 = [ LineMin2; [ tau, jj ] ];
end
  
end
