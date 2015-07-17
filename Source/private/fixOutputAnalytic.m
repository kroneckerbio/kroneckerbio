function expression = fixOutputAnalytic(expression)

if ischar(expression)
    % pass
elseif isnumeric(expression) && isscalar(expression)
    expression = num2str(expression);
else
    error('fixOutputAnalytic: invalid output expression')
end