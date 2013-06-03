function out = dealfield2cell(obj, field)
%DEALFIELD2CELL Move all contents of a field to a cell vector
%
%   out = dealfield2cell(obj,field)
%
%   For each member of obj, the contents of specified field are copied to a
%   cell vector.

n = numel(obj);

out = cell(n,1);
for i = 1:n
    out{i} = obj(i).(field);
end
