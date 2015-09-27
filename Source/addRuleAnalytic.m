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
nz = m.add.nz + 1;
m.add.nz = nz;
m.add.Rules = growRules(m.add.Rules, m.add.nz);

% Add item
m.add.Rules(nz).Name = name;
m.add.Rules(nz).Expression = fixRuleExpression(expression);

m.Ready = false;
