function expression = fixOutputAnalytic(expression)

if ischar(expression)
    warnAboutQuotedFullName(expression)
elseif isnumeric(expression) && isscalar(expression)
    expression = num2str(expression);
else
    error('KroneckerBio:Output:expression', 'Invalid output expression')
end
