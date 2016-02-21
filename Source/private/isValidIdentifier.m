function tf = isValidIdentifier(strings)
if ~iscell(strings)
    strings = {strings};
end

tf = ~cellfun(@isempty, regexp(strings, '^[A-Za-z_][A-Za-z_0-9]*([.][A-Za-z_][A-Za-z_0-9]*)?$', 'once'));
