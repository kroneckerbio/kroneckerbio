function val = piecewiselinear(t, tvalues, pvalues)
% PIECEWISELINEAR interpolates within a series of tvalues and
%   corresponding values.
%
%   val = PiecewiseLinear(t, tvalues, values)
%
%   Inputs
%       t          - A vector at which the value will be queried
%       tvalues    - A vector containing the tvalues
%       values     - A matrix with length(tvalues) rows and an arbitrary
%                    number of columns
%
%   Outputs
%       val - A array length(t) by size(values,2) containing the linear 
%             interpolation of the piecewise function. Values of t that are
%             outside the range of tvalues are equal to the nearest
%             value of values.
%
%   Examples:
%       a = piecewiselinear(0.5, [0 1], [2 1])
%       a = 1.5
%       b = piecewiselinear(0, [1 2], [2 1]) %Outside range
%       b = 2
%       c = piecewiselinear([0 0.1 0.2], [0 1], [0 2]) %Multiple queries
%       c = [0; 0.2; 0.4]
%       d = piecewiselinear(1, 0, 3) %Single data point
%       d = 3
%       e = piecewiselinear(0.5, [0;1], [0,1;0,2]) %Multiple outputs
%       e = [0 1.5];

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

count = numel(t);
val = zeros(count,size(pvalues,2));

for ind = 1:count
    k  = find(t(ind) >= tvalues, 1, 'last');
    if isempty(k)
        % Not bigger than any of them, so take the first
        val(ind,:) = pvalues(1,:);
        % Is greater than or equal to the largest one, so take the last
    elseif k == length(tvalues);
        val(ind,:) = pvalues(k,:);
    else
        % Is between tvalues
        dt = t(ind) - tvalues(k);
        T  = tvalues(k+1) - tvalues(k);
        
        val(ind,:) = pvalues(k,:)*(1-dt/T) + pvalues(k+1,:)*(dt/T);
    end
end
