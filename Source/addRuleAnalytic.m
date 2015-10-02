function m = addRuleAnalytic(m, name, expression)
% Add a rule to a Model.Analytic. Separate rules field simplifies conversion
% between kroneckerbio and SBML formats.
% Inputs:
%   m [ Model.Analytic struct ]
%       The model to which the rule will be added
%   name [ string ]
%       Name of rule
%   expression [ string ]
%       RHS of rule
%
% Inputs:
%   m [ Model.Analytic struct ]
%       The model with the new rule added

% Increment counter
nz = m.nz + 1;
m.nz = nz;
m.Rules = growRules(m.Rules, m.nz);

% Add item
m.Rules(nz).Name       = name;
m.Rules(nz).Expression = fixRuleExpression(expression);

m.Ready = false;
