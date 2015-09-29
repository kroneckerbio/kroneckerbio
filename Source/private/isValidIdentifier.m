function tf = isValidIdentifier(strings)

tf = ~cellfun(@isempty, regexp(strings, '^[A-Za-z_][A-Za-z_0-9]*$', 'once'));
