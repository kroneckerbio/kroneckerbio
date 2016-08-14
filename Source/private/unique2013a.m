function [ C, ia, ic ] = unique2013a( varargin )
%UNIQUE Find unique elements of vector
%   Compatibility form that reproduces the behavior of unique in R2013a and
%   later versions in all versions of MATLAB. Should be placed in the
%   private folder so that only Kronecker functions will call this version
%   of unique.

versionYear = version('-release');
versionYear = str2double(versionYear(1:4));
if versionYear > 2012
    [C, ia, ic] = unique(varargin{:});
    return
end

% Determine whether 'stable' argument was provided
isStable = strcmp(varargin, 'stable');
varargin(isStable) = [];
isStable = any(isStable);

assert(~any(strcmp(varargin, 'legacy')), 'KroneckerBio:unique:legacyNotSupported', ...
    '''legacy'' option currently not implemented in this version of ''unique''.')

[b,m,n] = unique(varargin{:}, 'first');
m = m(:);
n = n(:);

if isStable
    oldm = m;
    m = sort(m);
    [unused, mmap] = ismember(oldm, m);
    oldn = n;
    n = zeros(size(n));
    for ii = 1:numel(oldm)
        n(oldn == ii) = mmap(ii);
    end
    if isvector(varargin{1})
        b = varargin{1}(m);
    else
        b = varargin{1}(m,:);
    end
end

C = b;
ia = m;
ic = n;

end

