function AddParameter(this, name, value)
%AddParameter Add a rate parameter to a model. The parameter must be a
%nonnegative scalar.
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
% TODO: rate parameter must be a nonnegative scalar for massaction models but
% not for analytic models - generalize

% Increment counter
nk = this.nk + 1;
this.nk = nk;
this.growParameters;

% Add item
this.Parameters(nk).Name  = FieldValidator.ParameterName(name);
this.Parameters(nk).Value = FieldValidator.ParameterValue(value);

this.Ready = false;
