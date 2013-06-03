function [UseParams, nTk] = fixUseParams(UseParams, nk)
%FIXUSEICS A helper function for functions that support variently active
%   parameters to standardize them in a column vector of logicals.
%
%   [UseParams, nTk] = fixUseParams(UseParams, nk)
%
%   Inkuts
%       UseParams - Can be either of the following:
%           If UseModelICs = true
%               1) vector of linear indexes
%               2) vector of logical indexes length of nX
%       nk - Scalar natural number for how many total rate parameters are
%            in this model
%
%   Outputs
%       UseParams - A column vector of logicals length of nk
%       nTk - Number of active rate parameters

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if islogical(UseParams)
    UseParams = vec(UseParams);
else%if isnumeric
    temp = UseParams;
    UseParams = zeros(nk, 1);
    UseParams(temp) = 1;
    UseParams = logical(UseParams);
end
nTk = sum(UseParams);
