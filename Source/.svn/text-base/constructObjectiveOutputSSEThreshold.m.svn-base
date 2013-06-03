function obj = constructObjectiveOutputSSEThreshold(m, outputIndex, data, threshold, tss)
% 
% 
% function obj = constructObjectiveOutputSSE(m, outputIndex, data)
% 
% This function constructs an objective function that computes the
% sum of squared error between the a model output to data points provided
% in the matrix 'data'.
% 
% Inputs
%       m               -   the model from which output is collected
%       outputIndex     -   the index of the output for which the SSE will
%                           be computed
%       data            -   a structure with two components:
%                           data.t are the data timepoints
%                           data.y are the data concentrations
% 
% Outputs
%       obj             -   an objective function structure containing
%                           three functions:
%                           obj.G is the objective function value
%                           obj.dGdx is the derivative with respect to species
%                           obj.dGdp is the derivative with respect to parameters
% 
% See also:
%   constructObjectiveOutputValue
%   constructObjectiveOutputISE
%   constructObjectiveSpeciesValue
%   constructObjectiveSpeciesISE
%   constructObjectiveSpeciesSSE
%   constructObjectiveTimingSpeciesValue
%   constructObjectiveTimingSpeciesSSE
%   constructObjectiveTimingOutputValue
%   constructObjectiveTimingOutputSSE

% This work is licensed under the Creative Commons Attribution-Noncommercial-No Derivative 
% Works 3.0 United States License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-nd/3.0/us/ or send a letter to Creative Commons,
% 171 Second Street, Suite 300, San Francisco, California, 94105, USA.

% argument processing and error checking
if nargin<2 || isempty(outputIndex)
    outputIndex = 1;
end

assert(isnumeric(outputIndex), '', 'The second argument to constructObjectiveOutputISE must be an output index.');
assert(size(data.y, 1) == length(outputIndex), '', 'The length of outputIndex must be the same as the number of outputs in data.y.');
assert(outputIndex <= m.nY, '', 'The output index %d is greater than the total number of outputs %d.', outputIndex, m.nY);

% steady state times
ssTimes = linspace(tss/2, tss, length(data.t));

obj.update = @update;
obj.G = @G;
obj.dGdx = @dGdx;
obj.dGdp = @dGdp;

%% update function
    function objNew = update(mNew)
        objNew = constructObjectiveOutputSSEThreshold(mNew, outputIndex, data, threshold, tss);
    end
%% The objective function
    function [val stopTimes] = G(t,xSol,u,nSim)
        % evaluate the ODE solver's xSol structure to get x(time)
        x = deval(xSol(nSim), data.t+tss, 1:m.nX);
        % compute the difference between x and the corresponding data
        
        val = 0;
        for i = 1:length(data.t)
            tmp       = (m.c(outputIndex,:)*x(:,i) - data.y(:,i)).*(m.c(outputIndex,:)*x(:,i) > threshold | data.y(:,i) > 1.1*threshold);
            val       = val + tmp.'*tmp;
        end
        
        % get steady state trajectory values
        xss = deval(xSol(nSim), ssTimes, 1:m.nX);
        
        for i = 1:length(data.t)
            tmp       = (m.c(outputIndex,:)*xss(:,i) - data.y(:,1)).*(m.c(outputIndex,:)*xss(:,i) > threshold | data.y(:,1) > 1.1*threshold);
            val       = val + tmp.'*tmp;
        end
        
        stopTimes = [data.t+tss ssTimes];
        
        figure(2)
        plot(data.t, m.c(outputIndex,:)*x, data.t, data.y)
        figure(1)
    end

    function val = dGdx(t,xSol,u,nSim)
        currentTime = find(t == data.t+tss);
        if ~isempty(currentTime)
            % evaluate the ODE solver's xSol structure to get x(time)
            x = deval(xSol(nSim), t, 1:m.nX);
            val = 2*((m.c(outputIndex, :)*x - data.y(:,currentTime)).*(m.c(outputIndex,:)*x > threshold | data.y(:,currentTime) > 1.1*threshold)).'*m.c(outputIndex,:);
        else
            val = zeros(1, m.nX);
        end
        
        ssTime = find(t == ssTimes);
        if ~isempty(ssTime)
            x = deval(xSol(nSim), t, 1:m.nX);
            val = 2*((m.c(outputIndex, :)*x - data.y(:,1)).*(m.c(outputIndex,:)*x > threshold | data.y(:,1) > 1.1*threshold)).'*m.c(outputIndex,:);
        end
    end

    function val = dGdp(t,xSol,u,nSim)
        val = zeros(1, m.nP);
    end
end
