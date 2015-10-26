function expression = RuleExpression(expression)

if ischar(expression)
    % pass
else
    error('KroneckerBio:FieldValidator:RuleExpression', 'Invalid rule expression')
end