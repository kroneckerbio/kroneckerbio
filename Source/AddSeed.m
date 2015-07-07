function m = AddSeed(m, varargin)
%AddSeed Add a seed parameter to a generic model
%
%   m = AddSeed(m, ...)
%
%   This is a generic function. Select the appropriate help file below
%   depending on the type of the model.
%
%   Model.MassActionAmount
%       help addSeedMassActionAmount
%
%   Model.MassActionConcentration
%       help addSeedMassActionConcentration
%
%   Model.Analytic
%       help addSeedAnalytic

% (c) 2015 David R Hagen
% This work is released under the MIT license.

if is(m, 'Model.MassActionAmount')
    m = addSeedMassActionAmount(m, varargin{:});
elseif is(m, 'Model.MassActionConcentration')
    m = addSeedMassActionConcentration(m, varargin{:});
elseif is(m, 'Model.Analytic')
    m = addSeedAnalytic(m, varargin{:});
else
    error('KroneckerBio:AddSeed:m', 'm must be a model')
end
