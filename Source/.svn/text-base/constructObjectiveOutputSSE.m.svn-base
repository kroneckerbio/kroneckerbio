function obj = constructObjectiveOutputSSE(m, outputIndex, data)
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

obj.update = @update;
obj.G = @G;
obj.dGdx = @dGdx;
obj.dGdp = @dGdp;
obj.d2Gdx2 = @d2Gdx2;
obj.d2Gdp2 = @d2Gdp2;
obj.d2Gdpdx = @d2Gdpdx;
obj.d2Gdxdp = @d2Gdxdp;

%% update function
    function objNew = update(mNew)
        objNew = constructObjectiveOutputSSE(mNew, outputIndex, data);
    end

%% The objective function
    function [val stopTimes] = G(t,xSol,u,nSim)
        val = 0;
        stopTimes = [];
        if isempty(data.t), return, end
        
        % evaluate the ODE solver's xSol structure to get x(time)
        x = deval(xSol(nSim), data.t, 1:m.nX);
        % compute the difference between x and the corresponding data
        
        val = 0;
        for i = 1:length(data.t)
            tmp       = (m.c(outputIndex,:)*x(:,i) - data.y(:,i));
            val       = val + tmp.'*tmp;
        end
        
        stopTimes = data.t;
    end

    function val = dGdx(t,xSol,u,nSim)
        val = zeros(1, m.nX);
        if isempty(data.t), return, end
        
        currentTime = find(t == data.t);
        if ~isempty(currentTime)
            % evaluate the ODE solver's xSol structure to get x(time)
            x = deval(xSol(nSim), t, 1:m.nX);
            val = 2*(m.c(outputIndex, :)*x - data.y(:,currentTime)).'*m.c(outputIndex,:);
        else
            val = zeros(1, m.nX);
        end
    end

    function val = dGdp(t,xSol,u,nSim)
        val = zeros(1, m.nP);
    end

    function val = d2Gdx2(t,xSol,u,nSim)
        val = zeros(m.nX, m.nX);
        if isempty(data.t), return, end
        
        currentTime = find(t == data.t);
        if ~isempty(currentTime)
            % evaluate the ODE solver's xSol structure to get x(time)
            x = deval(xSol(nSim), t, 1:m.nX);
            val = 2*m.c(outputIndex, :).'*m.c(outputIndex, :);
        else
            val = zeros(m.nX, m.nX);
        end
    end

    function val = d2Gdp2(t,xSol,u,nSim)
        val = zeros(m.nP, m.nP);
    end

end
