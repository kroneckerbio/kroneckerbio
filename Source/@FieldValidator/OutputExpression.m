function expression = OutputExpression(expression)

if ischar(expression)
    % pass
elseif isnumeric(expression) && isscalar(expression)
    expression = num2str(expression);
else
    error('KroneckerBio:FieldValidator:OutputExpression', 'Invalid output expression')
end