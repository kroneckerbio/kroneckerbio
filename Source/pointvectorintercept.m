function [p scale] = pointvectorintercept(p1, v1, p2, v2)
%POINTVECTORINTERCEPT takes two points and two vectors eminating from them 
%   and computes the point of intercept and the scaling factors of each
%   vector needed to reach that point
%
%   [p scale] = pointvectorintercept(p1, v1, p2, v2)
%
%   Given two point-vector pairs, p1 and v1, and p2 and v2, this function
%   computes p, the point at which the two vectors eminating from their
%   respective points will cross. The function also returns the scaling of
%   each vector required to reach that point. The scaling could be
%   negative.
%
%   Multiple problems can be solved with one call: all the inputs must be n
%   by 2, where n is the number of independent problems to be solved. The
%   outputs will be n by 2 as well.
%
%   Point-vector pairs that are colinear return NaN for all values. Pairs
%   that are parallel but not colinear return inf for for the scales and
%   NaN for the point. A zero vector exhibits unpredictable behavior.
%
%   Example:
%
%   [p scale] = pointvectorintercept([0 -2; 0 1; 0 1], ...
%               [1 0; 1 1; -1 1], [1 0; 1 0; 1 0], [0 1; 1 1; -1 1])
%   p =
%        1     2
%      Inf   Inf
%      NaN   NaN
%   scale =
%       -1     2
%      Inf   Inf
%      NaN   NaN

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Get number of iterations
n = size(p1,1);

% Loop over problems problems
p = zeros(n,2);
scale = zeros(n,2);
for i = 1:n
    % For ease of reading equations, define these variables
    a = p1(i,1);
    b = p1(i,2);
    c = v1(i,1);
    d = v1(i,2);
    e = p2(i,1);
    f = p2(i,2);
    g = v2(i,1);
    h = v2(i,2);
    
    % Solve problem
    alpha = ( (e-a)*h + (b-f)*g ) / (c*h-d*g);
    if isfinite(alpha)
        point = p1(i,:) + alpha * v1(i,:);
        %beta = nanmean((point - p2(i,:)) ./ v2(i,:)); % slow
        if g ~= 0
            beta = (point(1) - e) / g;
        else%g == 0
            beta = (point(2) - f) / h;
        end
    else
        point = [NaN NaN];
        beta = alpha;
    end
    
    % Store results
    p(i,:) = point;
    scale(i,:) = [alpha beta];
end
