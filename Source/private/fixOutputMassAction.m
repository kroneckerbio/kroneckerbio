function expression = fixOutputMassAction(expression)

if ischar(expression)
    expression = {expression, 1};
elseif iscellstr(expression)
    expression = [vec(expression), num2cell(ones(numel(expression),1))];
else
    assert(size(expression,2) == 2 && isnumeric([expression{:,2}]), ...
        'KroneckerBio:Output:InvalidContributors', 'Output has an invalid list of contributors')
end
