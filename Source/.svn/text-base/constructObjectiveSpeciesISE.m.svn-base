function obj = constructObjectiveSpeciesISE(m, speciesIndex, data)
% 
% 
% function obj = constructObjectiveSpeciesISE(m, speciesIndex, data)
% 
% This function constructs an objective function that computes the
% sum of squared error between a set of model species and a continuous data
% function provided as a function handle.
% 
% Inputs
%       m               -   the model from which output is collected
%       speciesIndex     -   the index of the output for which the ISE will
%                           be computed
%       data            -   a function of the form xData = data(t), where
%                           xData is a vector of the same size as
%                           speciesIndex and represents the x data at time
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
%   constructObjectiveOutputISE
%   constructObjectiveOutputSSE
%   constructObjectiveSpeciesValue
%   constructObjectiveSpeciesSSE
%   constructObjectiveTimingSpeciesValue
%   constructObjectiveTimingSpeciesSSE
%   constructObjectiveTimingOutputValue
%   constructObjectiveTimingOutputSSE

% This work is licensed under the Creative Commons Attribution-Noncommercial-No Derivative 
% Works 3.0 United States License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-nd/3.0/us/ or send a letter to Creative Commons,
% 171 Second Street, Suite 300, San Francisco, California, 94105, USA.

assert(isnumeric(speciesIndex), '', 'The second argument to constructObjectiveSpeciesISE must be an species index.');
assert(any(speciesIndex <= m.nX), '', 'The species index %d is greater than the total number of species %d.', speciesIndex, m.nX);

% Compute the matrix that selects out the interesting species.
C = eye(m.nX);
C = C(speciesIndex,:);

obj.update      = @update;
obj.g           = @g;
obj.dgdx        = @dgdx;
obj.dgdp        = @dgdp;
obj.d2gdpdx     = @d2gdpdx;
obj.d2gdxdp     = @d2gdxdp;
obj.d2gdx2      = @d2gdx2;
obj.d2gdp2      = @d2gdp2;

%% update function
    function objNew = update(mNew)
        objNew = constructObjectiveSpeciesISE(mNew, speciesIndex, data);
    end
    
    function val = g(t,x,u,nSim)
        val = C*x - data(t);
        val = val.'*val;
    end
    function val = dgdx(t,x,u,nSim)
        val = 2*(C*x - data(t)).'*C;
    end
    function val = dgdp(t,x,u,nSim)
        val = zeros(1, m.nP);
    end
    function val = d2gdx2(t,x,u,nSim)
        val = 2*C.'*C;
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