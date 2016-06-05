function expression = fixRuleExpression(expression)

if ischar(expression)
    warnAboutQuotedFullName(expression)
else
    error('KroneckerBio:Rule:expression', 'Invalid rule expression')
end
