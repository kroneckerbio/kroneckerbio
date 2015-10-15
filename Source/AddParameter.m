function m = AddParameter(m, name, value)
%AddParameter Add a rate parameter to a model. The parameter must be a
%nonnegative scalar.
%
%   m = addParameter(m, name, value)
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

% Increment counter
nk = m.nk + 1;
m.nk = nk;
m.Parameters = growParameters(m.Parameters, nk);

% Add item
m.Parameters(nk).Name  = fixParameterName(name);
m.Parameters(nk).Value = fixParameterValue(value);

m.Ready = false;
