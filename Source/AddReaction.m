function m = AddReaction(m, varargin)
%AddReaction Add a reaction to a generic model
%
%   m = AddReaction(m, ...)
%
%   This is a generic function. Select the appropriate help file below
%   depending on the type of the model.
%
%   Model.MassActionAmount
%       help addReactionMassActionAmount
%
%   Model.MassActionConcentration
%       help addReactionMassActionConcentration
%
%   Model.Analytic
%       help addReactionAnalytic

% (c) 2015 David R Hagen
% This work is released under the MIT license.

if is(m, 'Model.MassActionAmount')
    m = addReactionMassActionAmount(m, varargin{:});
elseif is(m, 'Model.MassActionConcentration')
    m = addReactionMassActionConcentration(m, varargin{:});
elseif is(m, 'Model.Analytic')
    m = addReactionAnalytic(m, varargin{:});
else
    error('KroneckerBio:AddReaction:m', 'm must be a model')
end
