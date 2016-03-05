function expression = fixRuleExpression(expression)

if ischar(expression)
    warnAboutQuotedFullName(expression)
else
    error('fixRuleExpression: invalid rule expression')
end
