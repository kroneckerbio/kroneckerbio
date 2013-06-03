function com = compareData(outputs, times)
%COMPAREDATA Create a compare structure that simply samples the simulations
%   at the given times and outputs
%
%   com = compareData(outputs, times)
%
%   Inputs
%       outputs - Vector of ouput indexes where the measurements will be
%                 taken. There should be one entry for each measurement.
%       times   - Vector of time point at which the measurements will be
%                 taken. The length must be the same as outputs.
%
%   Outputs
%       com - Kroncker Bio simulation comparison structure

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Constants
n = numel(outputs);

% Find unique times
discreteTimes = vec(unique(times)).';

com.empty = @empty;
com.compare = @compare;
com.n = n;
com.continuous = false;
com.complex = false;
com.discreteTimes = discreteTimes;

    function val = empty(N)
        val = cell(N,2);
    end

    function result = compare(sol1, sol2, type)
        nX1 = size(sol1.c, 2);
        nX2 = size(sol2.c,2);
        
        % Extract data
        x1 = sol1.y(1:nX1,:); % x_t
        x2 = sol2.y(1:nX2,:); % x_t
        
        % Convert to outputs
        y1 = zeros(n,1);
        y2 = zeros(n,1);
        for i = 1:n
            y1(i,1) = sol1.c(outputs(i),:) * x1(:,discreteTimes == times(i));
            y2(i,1) = sol2.c(outputs(i),:) * x2(:,discreteTimes == times(i));
        end
        
        result = {y1, y2};
    end

end