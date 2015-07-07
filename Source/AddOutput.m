function m = AddOutput(m, varargin)
%AddOutput Add an output to a generic model
%
%   m = AddOutput(m, ...)
%
%   This is a generic function. Select the appropriate help file below
%   depending on the type of the model.
%
%   Model.MassActionAmount
%       help addOutputMassActionAmount
%
%   Model.MassActionConcentration
%       help addOutputMassActionConcentration
%
%   Model.Analytic
%       help addOutputAnalytic

% (c) 2015 David R Hagen
% This work is released under the MIT license.

if is(m, 'Model.MassActionAmount')
    m = addOutputMassActionAmount(m, varargin{:});
elseif is(m, 'Model.MassActionConcentration')
    m = addOutputMassActionConcentration(m, varargin{:});
elseif is(m, 'Model.Analytic')
    m = addOutputAnalytic(m, varargin{:});
else
    error('KroneckerBio:AddOutput:m', 'm must be a model')
end
