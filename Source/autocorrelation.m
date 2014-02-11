function acf = autocorrelation(y, lags, dim)
%AUTOCORRELATION Autocorrelation of sample
%
%   acf = autocorrelation(y, lags, dim)
%
%   The sample autocorrelation of sample y is computed for lag points lags
%   along dimension dim. If lags is missing, the autocorrelation is
%   computed for all available points. If dim is missing, the
%   autocorrelation is computed along the first non-singleton dimension.
%
%   Lags is 1-indexed so that the legal values are from 1 to size(y, dim).
%   This is different from the Matlab econometrics autocorr function which
%   uses 0 for the first lag.

% (c) 2014 David R Hagen and Bruce Tidor
% This work is released under the MIT license.

if nargin < 3
    dim = find(size(y) ~= 1, 1);
    if isempty(dim)
        dim = 1;
    end
    
    if nargin < 2
        lags = 1:size(y,dim);
    end
end

assert(size(y,dim) >= 1, 'KroneckerBio:autocorrelation:EmptyDimension', 'Dimension along which autocorrelation is computed must not be empty')

n_dims_y = max(ndims(y), dim);

% Center data
y_centered = bsxfun(@minus, y, mean(y, dim));

% Compute autocorrelation according to FFT method
FR = fft(y_centered, 2^(nextpow2(size(y,dim))+1), dim);
S = FR .* conj(FR);
acf = ifft(S, [], dim);

% Normalize to first element
first_elements = cell(n_dims_y,1);
for i = 1:n_dims_y
    if i == dim
        first_elements{i} = 1;
    else
        first_elements{i} = 1:size(y,i);
    end
end

acf = bsxfun(@rdivide, acf, acf(first_elements{:}));

% Keep only desired values
keep_elements = cell(n_dims_y,1);
for i = 1:n_dims_y
    if i == dim
        keep_elements{i} = lags;
    else
        keep_elements{i} = 1:size(y,i);
    end
end

acf = acf(keep_elements{:});