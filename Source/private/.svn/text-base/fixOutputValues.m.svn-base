function values = fixOutputValues(nExpr, values, units)

defaultValue = 1;

if isempty(values)
    values = zeros(nExpr,1) + defaultValue;
end

assert(numel(values) == nExpr, 'KroneckerBio:fixOutputValues:WrongValuesLength', 'The number of values must be equal to the number of expressions')

for i = 1:numel(units)
    if ~isempty(units{i}); error; end
end