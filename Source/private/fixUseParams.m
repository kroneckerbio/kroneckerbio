function [UseParams, nTk] = fixUseParams(UseParams, nk)
%fixUseParams A helper function for functions that support variently active
%   kinetic parameters to standardize them in a column vector of logicals.
%
%   [UseParams, nTk] = fixUseParams(UseParams, nk)
%
%   Inputs
%   UseParams: [ logical vector nk | positive integer vector ]
%       Indicates the kinetic parameters that will be allowed to vary
%       during the optimization
%   nk: [ nonegative integer scalar ]
%       The number of kinetic parameters in the model
%
%   Outputs
%   UseParams: [ logical vector nk ]
%       Standard form of UseParams
%   nTk: [ nonnegative integer scalar ]
%       Number of active kinetic parameters

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if isnumeric(UseParams) && isscalar(UseParams) && isnan(UseParams)
    % Default
    UseParams = true(nk,1);
elseif islogical(UseParams)
    UseParams = vec(UseParams);
else%if isnumeric
    temp = UseParams;
    UseParams = zeros(nk, 1);
    UseParams(temp) = 1;
    UseParams = logical(UseParams);
end
nTk = sum(UseParams);
