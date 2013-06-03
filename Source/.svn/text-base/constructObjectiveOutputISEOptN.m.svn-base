function obj = constructObjectiveOutputISEOptN(m, outputIndex, data)
% 
% 
% function obj = constructObjectiveOutputISEOptN(m, outputIndex, data)
% 
% This function constructs an objective function that computes the
% sum of squared error between a set of model outputs and a continuous data
% function provided as a function handle.
% 
% Inputs
%       m               -   the model from which output is collected
%       outputIndex     -   the index of the output for which the ISEOptN will
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
%   constructObjectiveSpeciesISEOptN
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

Omin = 1e-2;

obj.update      = @update;
obj.g           = @g;
obj.dgdx        = @dgdx;
obj.dgdp        = @dgdp;
obj.dgdO1       = @dgdO1;
obj.G           = @G;
obj.dGdx        = @dGdx;
obj.dGdp        = @dGdp;
obj.dGdO1       = @dGdO1;

%% update function
    function objNew = update(mNew)
        objNew = constructObjectiveOutputISEOptN(mNew, outputIndex, data);
    end
    
    function val = g(t,x,u,nSim,O1)
        val = diag(1./(O1+Omin))*m.c(outputIndex, :)*x - data(t);
        val = val.'*val;
    end
    function val = dgdx(t,x,u,nSim,O1)
        val = 2*(diag(1./(O1+Omin))*m.c(outputIndex, :)*x - data(t)).'*(diag(1./(O1+Omin))*m.c(outputIndex, :));
    end
    function val = dgdp(t,x,u,nSim,O1)
        val = zeros(1, m.nP);
    end
    function val = dgdO1(t,x,u,nSim,O1)
        val = 2*(diag(1./(O1+Omin))*m.c(outputIndex, :)*x - data(t)).'*(-diag(1./(O1+Omin).^2)*diag(m.c(outputIndex, :)*x));
    end

    function [val deltaTimes] = G(t,x,u,nSim,O1)
        val = 0;
        deltaTimes= [];
    end
    function val = dGdx(t,x,u,nSim,O1)
        val = zeros(1,m.nX);
    end
    function val = dGdp(t,x,u,nSim,O1)
        val = zeros(1,m.nP);
    end
    function val = dGdO1(t,x,u,nSim,O1)
        val = zeros(1,length(O1));
    end
end
