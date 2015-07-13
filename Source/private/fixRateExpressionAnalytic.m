function expr = fixRateExpressionAnalytic(expr)

if isempty(expr) || (isnumeric(expr) && isscalar(expr) && expr == 0)
    expr = [];
elseif ischar(expr)
    % pass
else
    error('fixRateExpressionAnalytic: invalid rate expression')
end
