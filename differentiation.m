function [dom] = differentiation(omega, h)
%%% -----------------------------------------------------------------
% to compute the derivative of omega with respect to c at a specific point
% using finite differencing techniques
% 4-point central differencing is used for middle points!
% 2-point forward or backward differencing is used for end points
% the derivatives are computed on grid points of state variables grid
% and the precomputed vector of derivaties at grid points are passed.
% derivatives at any other points need to be interpolated.
% A differentiation matrix is used for this purpose.
% Inputs: omega --> function values at descritized points
% Ouputs: derivative at descretized points, dom

global  Ngrid


% Occupy differentiation matrix
Dx = zeros(Ngrid,Ngrid);
Dx(1,1) = -12;
Dx(end,end) = 12;
v = ones(1,Ngrid-1)*8;
v(1) = 12; v(2) = 6; v(end) = 6;
Dx = Dx + diag(v,1);
v = ones(1,Ngrid-2)*(-1);
v(1) = 0; v(2) = 0;
Dx = Dx + diag(v,2);
v = ones(1,Ngrid-1)*(-8);
v(1) = -6; v(end-1) = -6; v(end) = -12;
Dx = Dx + diag(v,-1);
v = ones(1,Ngrid-2);
v(end) = 0; v(end-1) = 0;
Dx = Dx + diag(v,-2);
Dx = (1/(12*h)) * Dx;

% Multiplly differentiation matrix by function values

dom = Dx * transpose(omega);


end


