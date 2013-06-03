function [UseSeeds, nTs] = fixUseSeeds(UseSeeds, UseModelSeeds, ns, nCon)
%fixUseSeeds A helper function for functions that support variant initial
%   concentrations in experiments. It converts a variety of ways to
%   specific the active seed parameters and standardizes them into
%   either a vector (use model seeds) or matrix (use experiment seeds) of
%   logicals indicating the active parameters.
%
%   [UseSeeds, nTs] = fixUseSeeds(UseSeeds, UseModelSeeds, ns, nCon)
%
%   Inputs
%   UseSeeds: [ logical vector ns | logical matrix ns by nCon 
%               | positive integer vector ]
%       If UseModelSeeds = true
%           1) vector of linear indexes
%           2) vector of logical indexes length of ns
%       If UseModelSeeds = false
%           1) vector of linear indexes into ns, assumed same for all
%              conditions
%           2) vector of logical indexes length(ns) to indicate that
%              all conditions have the same active parameters
%           3) matrix of logical indexes size ns by nCon)
%   UseModelSeeds: [ logical scalar ]
%       Indicates that the model's seeds should be used, not the
%       experiments'.
%   nk: [ nonegative integer scalar ]
%       The number of seed parameters in the model
%       nCon - Scalar natural number for how many experiments are used in
%              this function
%
%   Outputs
%   UseSeeds: [ logical vector ns | logical matrix ns by nCon ]
%       If UseModelSeeds = true, then UseSeeds will be a logical column
%       vector. Otherwise, it will be a logical matrix.
%   nTs: [ nonnegative integer scalar ]
%       Number of active seed parameters

% (c) 2013 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

assert(islogical(UseSeeds) || isempty(UseSeeds) || (all(floor(UseSeeds) == UseSeeds) && all(UseSeeds >= 1)), 'KroneckerBio:UseSeeds:Integer', 'Some entries in UseSeeds are not natural numbers.')
assert(islogical(UseSeeds) || isempty(UseSeeds) || (max(UseSeeds) <= ns), 'KroneckerBio:UseSeeds:InvalidIndex', 'An entry in UseSeeds was larger than the number of seeds in the model.')

if UseModelSeeds
    if isnumeric(UseSeeds)
        assert(all(floor(UseSeeds) == UseSeeds) && all(UseSeeds >= 1), 'KroneckerBio:UseSeeds:InvalidValue', 'UseSeeds is an invalid linear index')
        assert(all(UseSeeds <= ns), 'KroneckerBio:UseControls:LinearIndexOutOfRange', 'If UseSeeds is provided as a linear index, then no linear index can be larger than m.ns')
        temp = false(ns, 1);
        temp(UseSeeds) = true;
        UseSeeds = temp;
    elseif islogical(UseSeeds)
        assert(numel(UseSeeds) <= ns, 'KroneckerBio:UseSeeds:InvalidLogicalLength', 'If UseSeeds is provided as a logical index, numel(UseSeeds) cannot be larger than m.ns')
        UseSeeds = vec(UseSeeds);
    else
        error('KroneckerBio:UseSeeds:InvalidType', 'UseSeeds must be provided as logical or linear index vector into m.s')
    end
else%if use seeds on conditions
    if isnumeric(UseSeeds)
        assert(all(floor(UseSeeds) == UseSeeds) && all(UseSeeds >= 1), 'KroneckerBio:UseSeeds:InvalidValue', 'UseSeeds is an invalid linear index')
        assert(all(UseSeeds <= ns), 'KroneckerBio:UseControls:LinearIndexOutOfRange', 'If UseSeeds is provided as a linear index, then no linear index can be larger than m.ns. Use a logical matrix to refer to different seeds on each condition.')
        temp = false(ns,nCon);
        temp(UseSeeds,:) = true;
        UseSeeds = temp;
    elseif islogical(UseSeeds)
        if numel(UseSeeds) <= ns
            % Parameter seeds are same for all conditions
            temp = false(ns,nCon);
            temp(UseSeeds,:) = true;
            UseSeeds = temp;
        elseif numel(UseSeeds) == ns*nCon
            % Parameter seeds are provided for all conditions. Format is already correct.
            UseSeeds = reshape(UseSeeds, ns,nCon);
        else
            error('KroneckerBio:UseSeeds:InvalidLogicalSize', 'If UseSeeds is provided as a logical matrix, numel(UseSeeds) must be less than m.ns or equal to m.ns*numel(con)')
        end
    else
        error('KroneckerBio:UseSeeds:InvalidType', 'UseSeeds must be provided as logical or linear index into con.s')
    end
end

nTs = sum(sum(UseSeeds));
