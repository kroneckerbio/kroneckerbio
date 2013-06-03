function m = AddInput(m, name, compartment, value, parameters, units)
%AddInput Add an input species to a KroneckerBio model
%
%   m = AddInput(m, name, compartment, value, parameters, units)
%
%   Inputs
%   m: [ model struct scalar ]
%       The model to which the input species will be added
%   name: [ string ]
%       A name for the input. This is the name by which reactions will
%       refer to it.
%   compartment: [ string ]
%       The name of the compartment it will be added
%   value: [ nonnegative scalar | string | 
%            function handle @(t,q) returns nonnegative scalar {0} ]
%       Optional
%       A scalar indicates that the input has a constant concentration, a
%       string is turned into a function handle of t and q, and the
%       function handle returns the value of the function at all time.
%   parameters: [ vector nq ]
%       This vector of parameters is passed to value(t,q) when determining
%       the value of the input at a particular time t.
%   units: [ '' ]
%       This is currently not implemented. The canonical representaion of
%       species values is in amount.
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new input added.

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean-up inputs
if nargin < 6
    units = [];
    if nargin <5
        parameters = [];
        if nargin < 4
            value = [];
        end
    end
end

% Increment counter
nxu = m.add.nxu + 1;
m.add.nxu = nxu;
m.add.Species = growSpecies(m.add.Species, nxu);

% Add item
m.add.Species(nxu).Name = fixName(name);
m.add.Species(nxu).Compartment = compartment;
m.add.Species(nxu).IsInput = true;
m.add.Species(nxu).Value.Function = fixInputValueFunction(value);
m.add.Species(nxu).Value.Parameters = fixInputParameters(parameters);
m.add.Species(nxu).Units = fixSpeciesUnits(units);

m.Ready = false;