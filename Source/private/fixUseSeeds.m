function [UseSeeds, nTs] = fixUseSeeds(UseSeeds, ns, nCon)
%fixUseSeeds A helper function for functions that support variant initial
%   concentrations in experiments. It converts a variety of ways to
%   specific the active seed parameters and standardizes them into
%   matrix of logicals indicating the active parameters.
%
%   [UseSeeds, nTs] = fixUseSeeds(UseSeeds, UseModelSeeds, ns, nCon)
%
%   Inputs
%   UseSeeds: [ logical matrix ns by nCon | positive integer vector ]
%       1) linear index vector into ns, assumed same for all conditions
%       2) logical index vector ns to indicate that all conditions have the
%           same active parameters
%       3) matrix of logical indexes size ns by nCon
%   ns: [ nonegative integer scalar ]
%       The number of seed parameters in the model
%   nCon: [ nonnegative integer ]
%       Number of experimental conditions
%
%   Outputs
%   UseSeeds: [ logical matrix ns by nCon ]
%       If UseModelSeeds = true, then UseSeeds will be a logical column
%       vector. Otherwise, it will be a logical matrix.
%   nTs: [ nonnegative integer ]
%       Number of active seed parameters

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if isnumeric(UseSeeds) && isscalar(UseSeeds) && isnan(UseSeeds)
    % Default
    UseSeeds = true(ns,1);
elseif isnumeric(UseSeeds) && all(floor(vec(UseSeeds)) == vec(UseSeeds)) && all(UseSeeds >= 1)
    % Linear index
    UseSeeds = vec(UseSeeds);
    assert(all(UseSeeds <= ns), 'KroneckerBio:UseSeeds:LinearIndexOutOfRange', 'UseSeeds, when a linear index, can have no value larger than m.ns. Use a logical matrix to refer to different seeds on each condition.')
    temp = false(ns,nCon);
    temp(UseSeeds,:) = true;
    UseSeeds = temp;
elseif islogical(UseSeeds)
    % Logical index
    assert(numel(UseSeeds) == ns*nCon, 'KroneckerBio:UseSeeds:InvalidLogicalSize', 'UseSeeds, when a logical index, must have a number of elements equal to ns*nCon')
    UseSeeds = reshape(UseSeeds, ns,nCon);
else
    error('KroneckerBio:UseSeeds:InvalidType', 'UseSeeds must be provided as logical or linear index into con.s')
end

nTs = nnz(UseSeeds);
