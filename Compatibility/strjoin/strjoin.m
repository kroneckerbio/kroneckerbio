function joined = strjoin(strings, delimiter)
%STRJOIN Join cell array of strings into one string seperated by delimiter
%
%   joined = strjoin(strings, delimiter)

strings = row(strings);

combined = [strings; repmat({delimiter}, 1,numel(strings)-1), {''}];

joined = cat(2, combined{:});
