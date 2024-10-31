%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function evaluates the L2 inner product of the two arguments;
% INPUT:
%  z1,z2    - the two arguments;   
%  which    - determines which state variable is used for computing the L2
%               inner product. The gradients are in respect to one of the state
%               variables --> 1 for C_1 and 2 for C_2
% OUTPUT
%  ip       -  value of the inner product;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ip] = prodL2(z1, z2, which)
  
global dC1 dC2

%ip = sum( z1(2:end-1) .* z2(2:end-1) ) + 0.5 * ( z1(1) * z2(1) + z1(end) * z2(end) );
ip = sum( z1 .* z2 );


if which == 1
    ip = ip * dC1;
else
    ip = ip * dC2;
end

  
return
