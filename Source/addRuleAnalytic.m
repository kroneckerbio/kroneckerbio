function m = addRuleAnalytic(m, name, expression)
%addRuleAnalytic Add a rule to a Model.Analytic. Rules are names associated
%   with expressions of compartments, parameters, inputs, states and other
%   rules.
%
%   m = addRuleAnalytic(m, name, expression)
%
%   Inputs:
%   m [ Model.Analytic struct ]
%       The model to which the rule will be added
%   name [ string ]
%       Name of rule
%   expression [ string ]
%       Value of rule as an expression of model components
%
%   Outputs:
%   m [ Model.Analytic struct ]
%       The model with the new rule added

% (c) 2015 Kevin Shi, David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Increment counter
nz = m.nz + 1;
m.nz = nz;
m.Rules = growRules(m.Rules, m.nz);

% Add item
m.Rules(nz).Name       = name;
m.Rules(nz).Expression = fixRuleExpression(expression);

m.Ready = false;
