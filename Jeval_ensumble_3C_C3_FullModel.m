function [jj] = Jeval_ensumble__FullModel(tau, omega_1, omega_2, lambda, DD, which)
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

global Feval LineMin1 LineMin2 c1tilde c2tilde I C0 timeseq

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


c1 = forward_ensumble_FullModel(omega_1_n, omega_2_n, lambda_n, timeseq.R2, I.R2, C0.R2);
jj1 = J(c1,c1tilde.R2,c2tilde.R2,timeseq.R2);

% c2 = forward_ensumble_FullModel(omega_1_n, omega_2_n, lambda_n, timeseq.R3, I.R3, C0.R3);
% jj2 = J(c2,c1tilde.R3,c2tilde.R3,timeseq.R3);
% 
% c3 = forward_ensumble_FullModel(omega_1_n, omega_2_n, lambda_n, timeseq.R4, I.R4, C0.R4);
% jj3 = J(c3,c1tilde.R4,c2tilde.R4,timeseq.R4);
% 
% c4 = forward_ensumble_FullModel(omega_1_n, omega_2_n, lambda_n, timeseq.R5, I.R5, C0.R5);
% jj4 = J(c4,c1tilde.R5,c2tilde.R5,timeseq.R5);

c5 = forward_ensumble_FullModel(omega_1_n, omega_2_n, lambda_n, timeseq.R6, I.R6, C0.R6);
jj5 = J(c5,c1tilde.R6,c2tilde.R6,timeseq.R6);

jj = (jj1 + jj5);

Feval = Feval + 1;

if which == 1
    LineMin1 = [ LineMin1; [ tau, jj ] ];
else
    LineMin2 = [ LineMin2; [ tau, jj ] ];
end
  
end
