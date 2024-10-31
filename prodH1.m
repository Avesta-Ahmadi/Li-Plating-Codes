%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function evaluates the H1 inner product of the two arguments
% weighted with the length-scale l of H1;
% INPUT:
%  z1,z2    - the two arguments;
% OUTPUT
%  ip       -  value of the inner product;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ip] = prodH1(z1, z2, which, l)

global dC1 dC2

if which == 1
    Dz1 = diff(z1) / dC1;
    Dz2 = diff(z2) / dC1;
    ip = prodL2(z1, z2, which) + l^2 * sum( Dz1 .* Dz2 ) * dC1;
else
    Dz1 = diff(z1) / dC2;
    Dz2 = diff(z2) / dC2;
    ip = prodL2(z1, z2, which) + l^2 * sum( Dz1 .* Dz2 ) * dC2;
end

return
