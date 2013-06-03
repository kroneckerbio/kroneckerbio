function obj = constructObjectiveOutputISEWeight(m, outputIndex, data, weight)
% 
% 
% function obj = constructObjectiveOutputISE(m, outputIndex, data)
% 
% This function constructs an objective function that computes the
% sum of squared error between a set of model outputs and a continuous data
% function provided as a function handle.
% 
% Inputs
%       m               -   the model from which output is collected
%       outputIndex     -   the index of the output for which the ISE will
%                           be computed
%       data            -   a function of the form yData = data(t), where
%                           yData is a vector of the same size as
%                           outputIndex and represents the y data at time
%                           t.
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
%   constructObjectiveOutputSSE
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

assert(isnumeric(outputIndex), '', 'The second argument to constructObjectiveOutputSSE must be an output index.');
assert(outputIndex <= m.nY, '', 'The output index %d is greater than the total number of outputs %d.', outputIndex, m.nY);

obj.update      = @update;
obj.g           = @g;
obj.dgdx        = @dgdx;
obj.dgdp        = @dgdp;
obj.d2gdx2      = @d2gdx2;
obj.d2gdp2      = @d2gdp2;
obj.d2gdpdx     = @d2gdpdx;
obj.d2gdxdp     = @d2gdxdp;

%% update function
    function objNew = update(mNew)
        objNew = constructObjectiveOutputISEWeight(mNew, outputIndex, data, weight);
    end

    function val = g(t,x,u,nSim)
        val = m.c(outputIndex, :)*x - data(t);
        val = weight*val.'*val;
    end
    function val = dgdx(t,x,u,nSim)
        val = weight*2*(m.c(outputIndex, :)*x - data(t)).'*m.c(outputIndex, :);
    end
    function val = dgdp(t,x,u,nSim)
        val = zeros(1, m.nP);
    end
    function val = d2gdx2(t,x,u,nSim)
        val = weight*2*m.c(outputIndex, :).'*m.c(outputIndex, :);
    end
    function val = d2gdp2(t,x,u,nSim)
        val = zeros(m.nP, m.nP);
    end
    function val = d2gdpdx(t,x,u,nSim)
        val = zeros(m.nX, m.nP);
    end
    function val = d2gdxdp(t,x,u,nSim)
        val = zeros(m.nP, m.nX);
    end
end
