function m = AddSeed(m, varargin)
%AddSeed Add a seed to a generic model
%
%   m = AddSeed(m, ...)
%
%   This is a generic function. Select the appropriate help file below
%   depending on the type of the model.
%
%   Model.MassActionAmount
%       help addSeedMassAction
%
%   Model.MassActionConcentration
%       help addSeedMassAction
%
%   Model.Analytic
%       help addSeedAnalytic

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if is(m, 'Model.MassActionAmount') || is(m, 'Model.MassActionConcentration')
    m = addSeedMassAction(m, varargin{:});
elseif is(m, 'Model.Analytic')
    m = addSeedAnalytic(m, varargin{:});
else
    error('KroneckerBio:AddCompartment:m', 'm must be a model')
end
