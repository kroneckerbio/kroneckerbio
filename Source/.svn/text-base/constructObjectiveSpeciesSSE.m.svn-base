function obj = constructObjectiveSpeciesSSE(m, speciesIndex, data)
% 
% 
% function obj = constructObjectiveSpeciesSSE(m, speciesIndex, data)
% 
% This function constructs an objective function that computes the
% sum of squared error between the a model species to data points provided
% in the matrix 'data'.
% 
% Inputs
%       m               -   the model from which species is collected
%       speciesIndex    -   the indices of the species for which the SSE will
%                           be computed
%       data            -   a structure with two components:
%                           data.t are the data timepoints
%                           data.x are the data concentrations
% 
% Speciess
%       obj             -   an objective function structure containing
%                           three functions:
%                           obj.G is the objective function value
%                           obj.dGdx is the derivative with respect to species
%                           obj.dGdp is the derivative with respect to parameters
% 
% See also:
%   constructObjectiveOutputValue
%   constructObjectiveOutputISE
%   constructObjectiveOutputSSE
%   constructObjectiveSpeciesValue
%   constructObjectiveSpeciesISE
%   constructObjectiveTimingSpeciesValue
%   constructObjectiveTimingSpeciesSSE
%   constructObjectiveTimingOutputValue
%   constructObjectiveTimingOutputSSE

% This work is licensed under the Creative Commons Attribution-Noncommercial-No Derivative 
% Works 3.0 United States License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-nd/3.0/us/ or send a letter to Creative Commons,
% 171 Second Street, Suite 300, San Francisco, California, 94105, USA.


% argument processing and error checking
if nargin<2 || isempty(speciesIndex)
    speciesIndex = 1;
end

assert(isnumeric(speciesIndex), '', 'The second argument to constructObjectiveSpeciesSSE must be an species index.');
assert(size(data.x, 1) == length(speciesIndex), '', 'The length of speciesIndex must be the same as the number of species in data.x.');
assert(speciesIndex <= m.nY, '', 'The species index %d is greater than the total number of speciess %d.', speciesIndex, m.nY);

% Compute the matrix that selects out the interesting species.
C = eye(m.nX);
C = C(speciesIndex,:);

obj.G = @G;
obj.dGdx = @dGdx;
obj.dGdp = @dGdp;

%% The objective function
    function [val stopTimes] = G(t,xSol,u,nSim)
        % evaluate the ODE solver's xSol structure to get x(time)
        x = deval(xSol(nSim), data.t, 1:m.nX);
        % compute the difference between x and the corresponding data
        val = 0;
        for i = 1:length(data.t)
            tmp       = (C*x(:,i) - data.x(:,i));
            val       = val + tmp'*tmp;
        end
        stopTimes = data.t;
    end

    function val = dGdx(t,xSol,u,nSim)
        currentTime = find(t == data.t);
        if ~isempty(currentTime)
            % evaluate the ODE solver's xSol structure to get x(time)
            x = deval(xSol(nSim), t, 1:m.nX);
            val = 2*(C*x - data.x(currentTime))'*C;
        else
            val = zeros(1, m.nX);
        end
    end

    function val = dGdp(t,xSol,u,nSim)
        val = zeros(1, m.nP);
    end
end