function out = is(obj, type)
% Returns true if obj.Type starts with type

if isfield(obj, 'Type') && all(strncmp({obj.Type}, type, numel(type)))
    out = true;
else
    out = false;
end
