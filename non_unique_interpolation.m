function [g] = non_unique_interpolation(c, grad, cgrid)
%%% -----------------------------------------------------------------
% This funciton perform interpolation when the mapping is not unique and
% one-to-one! in the case of computing the L2 gradients, one is required to
% follow time trajectory and compute gradient from there at different
% points in time (corresponding to different points in concentration)!
% Unlike time, concentration is not monotonic. We will have a map from
% concentration to gradient, where concentration is not monotonic. Hence,
% when interpolating gradient on the concentrtion grid, there might be more
% than one interval in the forward concentrtion data that could be used for
% interpolation. gradient needs to be interpolated based on each one of
% them and then added.
% Inputs: c is the concentration trajectory in time (which might be non-monotonic)
%         grad is the vector of gradients computed at for each concentration
%         cgrid is the grid point of concentration where the gradients need to be interpolated.
% Output: gradient at a particular concentration.
% Adjust lim size --> lim=5 means that 5 points prior and after to required
% point is used for interpolation. The mapping must be one-to-one in this
% interval.



% check to see if the point cgrid fits in more than one interval or not:
a = cgrid >= c(1:end-1);
b = cgrid < c(2:end);
d = a+b; % if 0 or 2 appears in d, then cgrid lies inside that interval.
% sometimes the concentrations fluctuates a bit in a small interval of
% time, in that case many 2 or 0 will appear in this vector in sequence. In
% such situations, the gradient will be interpolated and added to each
% other for each interval. This is wrong! To prevent, make sure when we see
% a 0 or 2 in the vector, we don't see any other 0 or 2 in the vector for
% some time. This prevents adding gradients many times by mistake.
for i=1:max(size(d))
    if d(i) == 0 || d(i) == 2
        d(i+1:i+5) = 1;
    end
end

indices2 = find(d==2);
indices0 = find(d==0);
g = 0;
lim = 2;
elim = 4; % 4 points for the interpolation close to the ends points

% Interpolating for the grid point and adding results if it is visited more than once!
if ~isempty(indices2)
    for i=1:max(size(indices2))
        if indices2(i)<= lim
            g = g + interp1(c(1:elim), grad(1:elim), cgrid,'spline');
        elseif indices2(i)>=(max(size(c))-lim)
            g = g + interp1(c(end-elim:end), grad(end-elim:end), cgrid,'spline');
        else
            g = g + interp1(c(indices2(i)-lim:indices2(i)+lim), grad(indices2(i)-lim:indices2(i)+lim), cgrid,'spline');
        end
    end
end
if ~isempty(indices0)
    for i=1:max(size(indices0))
        if indices0(i)<=lim
            g = g + interp1(c(1:elim), grad(1:elim), cgrid,'spline');
        elseif indices0(i)>=(max(size(c))-lim)
            g = g + interp1(c(end-elim:end), grad(end-elim:end), cgrid,'spline');
        else
            g = g + interp1(c(indices0(i)-lim:indices0(i)+lim), grad(indices0(i)-lim:indices0(i)+lim), cgrid,'spline');
        end
    end
end
if isempty(indices0) && isempty(indices2)
    g = 0;
end


end