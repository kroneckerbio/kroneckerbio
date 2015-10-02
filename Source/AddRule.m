function m = AddRule(m, varargin)
%AddRule Add a rule to a generic model
%
%   m = AddRule(m, ...)
%
%   This is a generic function. Select the appropriate help file below
%   depending on the type of the model.
%
%   Model.MassActionAmount
%       help addRuleMassActionAmount
%
%   Model.MassActionConcentration
%       help addRuleMassActionConcentration
%
%   Model.Analytic
%       help addRuleAnalytic

% (c) 2015 David R Hagen
% This work is released under the MIT license.

if is(m, 'Model.MassActionAmount')
    error('KroneckerBio:AddRule:MassActionAmount', 'Rules are not implemented for massaction models')
elseif is(m, 'Model.MassActionConcentration')
    m = addRuleMassActionConcentration(m, varargin{:});
elseif is(m, 'Model.Analytic')
    m = addRuleAnalytic(m, varargin{:});
else
    error('KroneckerBio:AddRule:m', 'm must be a model')
end
