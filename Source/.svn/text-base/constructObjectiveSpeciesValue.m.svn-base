function obj = constructObjectiveSpeciesValue(m, speciesIndex, time)
% 
% 
% function obj = constructObjectiveSpeciesValue(m, speciesIndex, time)
% 
% This function constructs an objective function that returns the value of
% an output at a specified time.
% 
% Inputs
%       m               -   the model from which output is collected.
%       speciesIndex     -   the index of the output of interest.
%       time            -   the time at which to compute the value of the
%                           output.
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
%   constructObjectiveOutputSSE
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
if nargin<2 || isempty(speciesIndex)
    speciesIndex = 1;
end

assert(isnumeric(speciesIndex), '', 'The second argument to constructEventSpeciesValue must be an species index.');
assert(length(speciesIndex) == 1, '', 'constructObjectiveSpeciesValue only accepts a scalar species index.');
assert(length(time) == 1, '', 'constructObjectiveSpeciesValue only accepts a single timepoint.');
assert(speciesIndex <= m.nY, '', 'The output index %d is greater than the total number of outputs %d.', speciesIndex, m.nY);

% Compute the matrix that selects out the interesting species.
C = eye(m.nX);
C = C(speciesIndex,:);

obj.update = @update;
obj.G = @G;
obj.dGdx = @dGdx;
obj.dGdp = @dGdp;

    function objNew = update(mNew)
        objNew = constructObjectiveSpeciesValue(mNew, speciesIndex, time);
    end

    function [val stopTimes] = G(t,xSol,u,nSim)
        % evaluate the ODE solver's xSol structure to get x(time)
        x = deval(xSol(nSim), time, 1:m.nX);
        
        % compute y(time)
        val       = C*x;
        stopTimes = time;
    end

    function val = dGdx(t,xSol,u,nSim)
        if t == time
            val = C;
        else
            val = zeros(1, m.nX);
        end
    end

    function val = dGdp(t,xSol,u,nSim)
        val = zeros(1, m.nP);
    end
end