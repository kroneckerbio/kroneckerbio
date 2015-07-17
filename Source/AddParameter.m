function m = AddParameter(m, varargin)
%AddParameter Add a rate parameter to a generic model
%
%   m = AddParameter(m, ...)
%
%   This is a generic function. Select the appropriate help file below
%   depending on the type of the model.
%
%   Model.MassActionAmount
%       help addParameterMassActionAmount
%
%   Model.MassActionConcentration
%       help addParameterMassActionConcentration
%
%   Model.Analytic
%       help addParameterAnalytic

% (c) 2015 David R Hagen
% This work is released under the MIT license.

if is(m, 'Model.MassActionAmount')
    m = addParameterMassActionAmount(m, varargin{:});
elseif is(m, 'Model.MassActionConcentration')
    m = addParameterMassActionConcentration(m, varargin{:});
elseif is(m, 'Model.Analytic')
    m = addParameterAnalytic(m, varargin{:});
else
    error('KroneckerBio:AddParameter:m', 'm must be a model')
end
