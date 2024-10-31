function [dom] = difference(cgrid, omega, cx)
%%% -----------------------------------------------------------------
% to compute the derivative of omega with respect to c at a specific point
% using finite differencing techniques
% 4-point central differencing is used for middle points!
% 2-point forward or backward differencing is used for exterior points (exterior of cgrid)
% Inputs: cgrid and omega are discretized, omega values at cgird are present
%         cx: point where we want to compute derivative
% Ouputs: derivative at a specific point, dom

global  Ngrid

% compute function value at cx  by interpolation: 
omcx = interp1(cgrid, omega, cx, 'spline');
h = (cgrid(end) - cgrid(1))/Ngrid;
n = sum(cgrid<cx);

if n>1 && n<size(cgrid,2)-1
    % 4-point central differencing
    dom = (omega(n-1) - 8*omega(n) + 8*omega(n+1) - omega(n+2) )/(12*h);
elseif n == 1 
    % 2-point central differencing
    dom = ( omega(n+1) - omega(n) )/(2*h);
elseif n == size(cgrid,2)-1
    % 2-point central differencing
    dom = ( omega(n+1) - omega(n) )/(2*h);
elseif n < 1
    % 2-point forward differencing
    dom = ( omega(n+1) - omcx )/( cgrid(1) - cx);
else
    % 2-point backward differencing
    dom = ( omcx - omega(n) )/( cx - cgrid(end) );
end

end


