function obj = keepfields(obj, fields)

assert(iscellstr(fields), 'KroneckerBio:keepfields:fields', 'fields must be a cell array of strings')

fields = vec(fields);
values = vec(cellfun(@(name)obj.(name), fields, 'UniformOutput', false));

% Alternate field ame with field value
args = row([fields, values].');

% Return updated struct with only these fields and values
obj = struct(args{:});
