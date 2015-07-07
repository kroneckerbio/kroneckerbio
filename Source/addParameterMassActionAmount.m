function m = addParameterMassActionAmount(m, name, value)
%AddParameter Add a rate parameter to a KroneckerBio model
%
%   m = AddParameter(m, name, value)
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
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new parameter added

% (c) 2013 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Increment counter
nk = m.add.nk + 1;
m.add.nk = nk;
m.add.Parameters = growParameters(m.add.Parameters, nk);

% Add item
m.add.Parameters(nk).Name = fixParameterName(name);
m.add.Parameters(nk).Value = fixParameterValue(value);

m.Ready = false;
