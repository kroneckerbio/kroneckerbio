function m = AddParameter(m, name, value, units)
%AddParameter Add a rate parameter to a KroneckerBio model
%
%   m = AddParameter(m, name, value, units)
%
%   Rate parameters determine the rate of reactions to which they are
%   associated.
%
%   Inputs
%   m: [ model struct scalar ]
%       The model to which the parameter will be added
%   name: [ string ]
%       A name for the parameter. This is the name by which reactions will
%       refer to it.
%   value: [ nonnegative scalar ]
%       The numeric value of the rate parameter.
%   units: [ '' ]
%       This is currently not implemented. The canonical representation of
%       rate parameters is in units of concentration.
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new compartment added.

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean-up inputs
if nargin < 4
    units = [];
end

% Increment counter
nk = m.add.nk + 1;
m.add.nk = nk;
m.add.Parameters = growParameters(m.add.Parameters, nk);

% Add item
m.add.Parameters(nk).Name = fixName(name);
m.add.Parameters(nk).Value = fixParameterValue(value, fixSpeciesUnits(units));

m.Ready = false;