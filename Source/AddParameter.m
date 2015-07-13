function m = AddParameter(m, varargin)
%AddParameter Add a parameter to a generic model
%
%   m = AddParameter(m, ...)
%
%   This is a generic function. Select the appropriate help file below
%   depending on the type of the model.
%
%   Model.MassActionAmount
%       help addParameterMassAction
%
%   Model.MassActionConcentration
%       help addParameterMassAction
%
%   Model.Analytic
%       help addParameterMassAnalytic

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if is(m, 'Model.MassActionAmount') || is(m, 'Model.MassActionConcentration')
    m = addParameterMassAction(m, varargin{:});
elseif is(m, 'Model.Analytic')
    m = addParameterAnalytic(m, varargin{:});
else
    error('KroneckerBio:AddCompartment:m', 'm must be a model')
end
