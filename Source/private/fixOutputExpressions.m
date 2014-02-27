function expressions = fixOutputExpressions(expressions)

if isempty(expressions)
    expressions = cell(0,2);
elseif isnumeric(expressions)
    expressions = {'', expressions(end)};
elseif ischar(expressions)
    expressions = {row(expressions), 1};
elseif iscellstr(expressions)
    expressions = [vec(expressions), num2cell(ones(numel(expressions),1))];
else
    assert(size(expressions,2) == 2 && isnumeric([expressions{:,2}]), ...
        'KroneckerBio:Output:InvalidContributors', 'Output has an invalid list of contributors')
end
