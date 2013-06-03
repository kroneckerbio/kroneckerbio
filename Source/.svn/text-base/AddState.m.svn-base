function m = AddState(m, name, compartment, initialValue, units)
%AddState Add a state species to a KroneckerBio model
%
%   m = AddState(m, name, compartment, initialValue, units)
%
%   Inputs
%   m: [ model struct scalar ]
%       The model to which the state will be added
%   name: [ string ]
%       A name for the state. This is the name by which reactions will
%       refer to it.
%   compartment: [ string ]
%       The name of the compartment it will be added
%   initialValue: [ nonnegative scalar {0} ]
%       Optional
%       This is the initial value of the state.
%   units: [ '' ]
%       This is currently not implemented. The canonical representaion of
%       species values is in amount.
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new state added.

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean-up inputs
if nargin < 5
    units = [];
    if nargin < 4
        initialValue = [];
    end
end

if isempty(initialValue)
    initialValue = 0;
end
if isempty(units)
    units = '';
end

% Increment counter
nxu = m.add.nxu + 1;
m.add.nxu = nxu;
m.add.Species = growSpecies(m.add.Species, nxu);

% Add item
m.add.Species(nxu).Name = fixName(name);
m.add.Species(nxu).Compartment = compartment;
m.add.Species(nxu).IsInput = false;
m.add.Species(nxu).Value = initialValue;
m.add.Species(nxu).Units = fixSpeciesUnits(units);

m.Ready = false;