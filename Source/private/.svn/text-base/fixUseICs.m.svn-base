function [UseICs, nTx] = fixUseICs(UseICs, UseModelICs, nx, nCon)
%FIXUSEICS A helper function for functions that support variant initial
%   concentrations in experiments. It converts a variety of ways to
%   specific the active initial conditions and standardizes them into
%   either a vector (use model ICs) or matrix (use experiment ICs) of
%   logicals indicating the active parameters.
%
%   [UseICs, nTx] = fixUseICs(UseICs, UseModelICs, nx, nCon)
%
%   Inputs
%       UseICs - Can be any of the following:
%           If UseModelICs = true
%               1) vector of linear indexes
%               2) vector of logical indexes length of nx
%           If UseModelICs = false
%               1) vector of linear indexes into nx, assumed same for all
%                  conditions
%               2) vector of logical indexes length(nx) to indicate that
%                  all conditions have the same active parameters
%               3) matrix of logical indexes size nx by nCon)
%       UseModelICs - Scalar logical indicating that the model initial
%                     conditions are active and constant for all 
%                     experiments
%       nx - Scalar natural number for how many species are in this model
%       nCon - Scalar natural number for how many experiments are used in
%              this function
%
%   Outputs
%       UseICs - If UseModelICs = true, then UseICs will be a logical
%           column vector length of nx. Otherwise, it will be a logical
%           matrix sizw nx by nCon.
%       nTx - Number of active IC parameters

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

assert(islogical(UseICs) || isempty(UseICs) || (all(floor(UseICs) == UseICs) && all(UseICs >= 1)), 'KroneckerBio:UseICs:Integer', 'Some entries in UseICs are not natural numbers.')
assert(islogical(UseICs) || isempty(UseICs) || (max(UseICs) <= nx), 'KroneckerBio:UseICs:InvalidIndex', 'An entry in UseICs was larger than the number of species in the model.')

if UseModelICs
    if isnumeric(UseICs)
        assert(all(floor(UseICs) == UseICs) && all(UseICs >= 1), 'KroneckerBio:UseICs:InvalidValue', 'UseICs is an invalid linear index')
        assert(all(UseICs <= nx), 'KroneckerBio:UseControls:LinearIndexOutOfRange', 'If UseICs is provided as a linear index, then no linear index can be larger than m.nx')
        temp = false(nx, 1);
        temp(UseICs) = true;
        UseICs = temp;
    elseif islogical(UseICs)
        assert(numel(UseICs) <= nx, 'KroneckerBio:UseICs:InvalidLogicalLength', 'If UseICs is provided as a logical index, numel(UseICs) cannot be larger than m.nx')
        UseICs = vec(UseICs);
    else
        error('KroneckerBio:UseICs:InvalidType', 'UseICs must be provided as logical or linear index vector into m.x0')
    end
else%if use ICs on conditions
    if isnumeric(UseICs)
        assert(all(floor(UseICs) == UseICs) && all(UseICs >= 1), 'KroneckerBio:UseICs:InvalidValue', 'UseICs is an invalid linear index')
        assert(all(UseICs <= nx), 'KroneckerBio:UseControls:LinearIndexOutOfRange', 'If UseICs is provided as a linear index, then no linear index can be larger than m.nx. Use a logical matrix to refer to different ICs on each condition.')
        temp = false(nx,nCon);
        temp(UseICs,:) = true;
        UseICs = temp;
    elseif islogical(UseICs)
        if numel(UseICs) <= nx
            % Parameter ICs are same for all conditions
            temp = false(nx,nCon);
            temp(UseICs,:) = true;
            UseICs = temp;
        elseif numel(UseICs) == nx*nCon
            % Parameter ICs are provided for all conditions. Format is already correct.
            UseICs = reshape(UseICs, nx,nCon);
        else
            error('KroneckerBio:UseICs:InvalidLogicalSize', 'If UseICs is provided as a logical matrix, numel(UseICs) must be less than m.nx or equal to m.nx*numel(con)')
        end
    else
        error('KroneckerBio:UseICs:InvalidType', 'UseICs must be provided as logical or linear index into [con.x0]')
    end
end
nTx = sum(sum(UseICs));