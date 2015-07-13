function expression = fixRuleExpression(expression)

if ischar(expression)
    % pass
else
    error('fixRuleExpression: invalid rule expression')
end