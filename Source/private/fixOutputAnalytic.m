function expression = fixOutputAnalytic(expression)

if ischar(expression)
    % pass
elseif isnumeric(expression) && isscalar(expression)
    expression = num2str(expression);
else
    error('KroneckerBio:Output:InvalidAnalyticExpression', 'Invalid output expression')
end