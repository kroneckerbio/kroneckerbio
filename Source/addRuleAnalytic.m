function m = addRuleAnalytic(m, name, expression, id)
% Add a rule to a Model.Analytic. Separate rules field simplifies conversion
% between kroneckerbio and SBML formats.
% Inputs:
%   m [ Model.Analytic struct ]
%       The model to which the rule will be added
%   name [ string ]
%       Name of rule
%   expression [ string ]
%       RHS of rule
%   id: [ string {[]}]
%       A unique, valid variable name
% Inputs:
%   m [ Model.Analytic struct ]
%       The model with the new rule added
%
% Note/TODO: other SBML rule types aren't implemented yet
% Note/TODO: initial assignment rules allow non-seed values (like k's) to form the
%   expression for an initial condition, but subs in constants at model
%   finalization time. This will cause the initial condition to be inconsistent if
%   those k's change, i.e., during fitting.
% Note: the `AddRule` function clashes with a function in the fuzzy logic
%   toolbox

if nargin < 4
    id = [];
end

if isempty(id)
    id = '';
elseif issym(id)
    id = char(id);
end

% Increment counter
nz = m.add.nz + 1;
m.add.nz = nz;
m.add.Rules = growRulesAnalytic(m.add.Rules, m.add.nz);

% Add item
m.add.Rules(nz).Name = name;
m.add.Rules(nz).ID = id;
m.add.Rules(nz).Expression = fixRuleExpression(expression);

m.Ready = false;
