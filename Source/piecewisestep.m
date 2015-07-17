function val = piecewisestep(t, tvalues, pvalues)
% PIECEWISESTEP returns the appropriate pvalues to a stepwise-defined
%   function.
%
%   val = piecewisestep(t, tvalues, pvalues)
%
%   Inputs
%       t          - A vector at which the value will be queried
%       tvalues    - A vector containing the tvalues
%       pvalues    - A matrix with length(tvalues) rows and an arbitrary
%                    number of columns
%
%   Outputs
%       val - An array of size(t) containing the stepwise interpretation of
%             the function at t; that is, the value of the timepoint
%             directly before. Values of t that are outside the range of
%             tvalues are equal to the nearest value in pvalues.
%
%   Examples:
%       a = piecewisestep(0.5, [0 1], [2 1])
%       a = 2
%       b = piecewisestep(0, [1 2], [2 1]) %Outside range
%       b = 2
%       c = piecewisestep([0 0.1 0.2], [0 1], [0 2]) %Multiple queries
%       c = [0, 0, 0]
%       d = piecewisestep(1, 0, 3) %Single data point
%       d = 3
%       e = piecewiselinear(0.5, [0;1], [0,1;0,2]) %Multiple outputs
%       e = [0 0];

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

count = numel(t);
val = zeros(count,size(pvalues,2));

for ind = 1:count
    k  = find(t(ind) >= tvalues, 1, 'last');
    
    %not bigger than any of them, so take the first
    if isempty(k)
        val(ind,:) = pvalues(1,:);
        
    %is greater than one of them, take its value
    else
        val(ind,:) = pvalues(k,:);
    end
end
