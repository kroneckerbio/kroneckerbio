function m = AddState(m, varargin)
%AddState Add a state species to a generic model
%
%   m = AddState(m, ...)
%
%   This is a generic function. Select the appropriate help file below
%   depending on the type of the model.
%
%   Model.MassActionAmount
%       help addStateMassActionAmount
%
%   Model.MassActionConcentration
%       help addStateMassActionConcentration
%
%   Model.Analytic
%       help addStateAnalytic

% (c) 2015 David R Hagen
% This work is released under the MIT license.

if is(m, 'Model.MassActionAmount')
    m = addStateMassActionAmount(m, varargin{:});
elseif is(m, 'Model.MassActionConcentration')
    m = addStateMassActionConcentration(m, varargin{:});
elseif is(m, 'Model.Analytic')
    m = addStateAnalytic(m, varargin{:});
else
    error('KroneckerBio:AddState:m', 'm must be a model')
end
