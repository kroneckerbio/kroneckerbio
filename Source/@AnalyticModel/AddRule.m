function AddRule(this, name, expression, type)
% Add a rule to a Model.Analytic. Separate rules field simplifies conversion
% between kroneckerbio and SBML formats.
%
% Inputs:
%   m [ scalar AnalyticModel ]
%       The model to which the rule will be added
%   name [ string ]
%       LHS of rule
%   expression [ string ]
%       RHS of rule
%   type [ {'repeated assignment'} | 'initial assignment' | 'rate' ]
%       Rule type, taken from SBML spec.

% Clean up inputs
if nargin < 5
    type = 'repeated assignment';
end

% Check valid rule type
validRuleTypes = {'repeated assignment', 'initial assignment', 'rate'};
assert(ismember(type, validRuleTypes), 'KroneckerBio:addRuleAnalytic:InvalidRuleType', 'Rule type %s is not a supported rule type', type)

% Increment counter
nz = this.nz + 1;
this.nz = nz;
this.growRules;

% Add item
this.Rules(nz).Name       = name;
this.Rules(nz).Expression = FieldValidator.RuleExpression(expression);
this.Rules(nz).Type       = type;

this.Ready = false;
