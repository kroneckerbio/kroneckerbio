function expression = fixRateExpressionAnalytic(expression)

if isempty(expression) || (isnumeric(expression) && isscalar(expression) && expression == 0)
    expression = '';
elseif ischar(expression)
    warnAboutQuotedFullName(expression)
else
    error('KroneckerBio:Reaction:forward_reverse', 'Invalid rate expression')
end
