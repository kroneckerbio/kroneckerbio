function s = emptystruct(varargin)

if isnumeric(varargin{1})
    % Dimensions are provided
    dim = varargin{1};
    
    % Correct for Matlab assuming scalar dimension means square matrix
    if isscalar(dim)
        dim = [dim, 1];
    end
    
    start = 2;
else
    % Default dimensions of 1
    dim = [1,1];
    start = 1;
end

% Number of entries in struct
nEntries = nargin - start + 1;

% Arguments for an empty struct of the appropriate size
entries = [varargin(start:end); repmat({cell(dim)}, 1,nEntries)];

% Build empty structure
s = struct(entries{:});