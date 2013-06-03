function obj = constructObjectiveOutputISE(m, outputIndex, data)
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

% normalize to the integral of each data trajectory
N = zeros(length(outputIndex), 1);
for i = 1:length(N)
    N(i) = quadgk(@(t) index(data(t)', (i-1)*length(t)+(1:length(t))), 0, 20);
end
N = 10*diag(1./N);
N(2,2) = 0.1;

obj.g           = @g;
obj.dgdx        = @dgdx;
obj.dgdp        = @dgdp;

    function val = g(t,x,u,nSim)
        val = N*(m.c(outputIndex, :)*x - data(t));
        val = conj(val')*val;
    end
    function val = dgdx(t,x,u,nSim)
        val = 2*conj((N*(m.c(outputIndex, :)*x - data(t)))')*(N*m.c(outputIndex, :));
    end
    function val = dgdp(t,x,u,nSim)
        val = zeros(1, m.nP);
    end
end