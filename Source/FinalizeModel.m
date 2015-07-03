function m = FinalizeModel(m, varargin)
%FinalizeModel Update the mathematical components of the model to reflect
%   changes made to the model
%
%   m = FinalizeModel(m, ...)
%
%   This is a generic function. Select the appropriate help file below
%   depending on the type of the model.
%
%   Model.MassActionAmount
%       help finalizeModelMassActionAmount
%
%   Model.MassActionConcentration
%       help finalizeModelMassActionConcentration
%
%   Model.Analytic
%       help finalizeModelAnalytic

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if is(m, 'Model.MassActionAmount')
    m = finalizeModelMassActionAmount(m, varargin{:});
elseif is(m, 'Model.MassActionConcentration')
    m = finalizeModelMassActionConcentration(m, varargin{:});
elseif is(m, 'Model.Analytic')
    m = finalizeModelAnalytic(m, varargin{:});
else
    error('KroneckerBio:AddOutput:m', 'm must be a model')
end
