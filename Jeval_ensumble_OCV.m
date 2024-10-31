function [jj] = Jeval_ensumble_OCV(tau, beta, DD, which)
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
%       beta      - containing the coefficients
%       DD        - Descent Direction 
%       which     - An integer (1 or 2 or 3 or 4) specifying which parameter is
%       getting updated.
% OUTPUT:
%       jj         -  cost function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Feval LineMin1 LineMin2 c1tilde c2tilde C0 timeseq
beta(which) = beta(which) + tau * DD;

c1 = forward_ensumble_OCV(beta, timeseq.R2, C0.R2);
jj1 = J(c1,c1tilde.R2,c2tilde.R2,timeseq.R2);

c2 = forward_ensumble_OCV(beta, timeseq.R3, C0.R3);
jj2 = J(c2,c1tilde.R3,c2tilde.R3,timeseq.R3);

c3 = forward_ensumble_OCV(beta, timeseq.R4, C0.R4);
jj3 = J(c3,c1tilde.R4,c2tilde.R4,timeseq.R4);

c4 = forward_ensumble_OCV(beta, timeseq.R5, C0.R5);
jj4 = J(c4,c1tilde.R5,c2tilde.R5,timeseq.R5);

c5 = forward_ensumble_OCV(beta, timeseq.R6, C0.R6);
jj5 = J(c5,c1tilde.R6,c2tilde.R6,timeseq.R6);

jj = (jj1 + jj2 + jj3 + jj4 + jj5);

Feval = Feval + 1;

if which == 1
    LineMin1 = [ LineMin1; [ tau, jj ] ];
else
    LineMin2 = [ LineMin2; [ tau, jj ] ];
end
  
end
