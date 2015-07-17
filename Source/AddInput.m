function m = AddInput(m, varargin)
%AddInput Add an input species to a generic model
%
%   m = AddInput(m, ...)
%
%   This is a generic function. Select the appropriate help file below
%   depending on the type of the model.
%
%   Model.MassActionAmount
%       help addInputMassActionAmount
%
%   Model.MassActionConcentration
%       help addInputMassActionConcentration
%
%   Model.Analytic
%       help addInputAnalytic

% (c) 2015 David R Hagen
% This work is released under the MIT license.

if is(m, 'Model.MassActionAmount')
    m = addInputMassActionAmount(m, varargin{:});
elseif is(m, 'Model.MassActionConcentration')
    m = addInputMassActionConcentration(m, varargin{:});
elseif is(m, 'Model.Analytic')
    m = addInputAnalytic(m, varargin{:});
else
    error('KroneckerBio:AddInput:m', 'm must be a model')
end
