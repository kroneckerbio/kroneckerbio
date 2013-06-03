function abstol = fixAbsTolSelf(AbsTol, order, integrateObj, nx, nCon, UseAdjoint, UseModelSeeds, UseModelInputs, UseParams, UseSeeds, UseControls)
%fixAbsTolSelf Standardize the presentation of AbsTol for
%   self-sensitivities
%
%   abstol = fixAbsTolSelf(absTol, order, integrateObj, ns, nCon, useAdjoint,
%   UseModelSeeds, useModelInputs, useParams, UseSeeds, useControls)

% (c) 2012 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean up inputs
if nargin < 11
    UseControls = {};
    if nargin < 10
        UseSeeds = [];
        if nargin < 9
            UseParams = [];
            if nargin < 8
                UseModelSeeds = [];
                if nargin < 7
                    UseModelInputs = [];
                    if nargin < 6
                        UseAdjoint = [];
                    end
                end
            end
        end
    end
end

% The default AbsTol
absTolDefault = 1e-9;
if isempty(AbsTol) || (~iscell(AbsTol) && ~isstruct(AbsTol) && all(vec(isnan(AbsTol))))
    AbsTol = absTolDefault;
end

% Extract struct inside cell
if iscell(AbsTol) && isstruct(AbsTol{1})
    assert(numel(AbsTol) == 1, 'KroneckerBio:AbsTol:TooManyStructsInCells', 'AbsTol was a cell with a struct inside, but there is more than one cell, which is not allowed.')
    AbsTol = AbsTol{1};
end

% Extract cell inside cell
if iscell(AbsTol) && iscell(AbsTol{1})
    assert(numel(AbsTol) == 1, 'KroneckerBio:AbsTol:TooManyCellsInCells', 'AbsTol was a cell array with a cell array inside, but there is more than one cell in the higher array, which is not allowed.')
    AbsTol = AbsTol{1};
end

% Constants
nTk = sum(UseParams);
nTs = sum(sum(UseSeeds));
nTq = sum(cat(1, UseControls{:}));
nT = nTk + nTs + nTq;

abstol = cell(nCon,1);

switch order
    case 1
        error('Self sensitivities not yet supported for first order.')
    case 2
        if ~iscell(AbsTol)
            if UseAdjoint
                error('Adjoint method with optimal AbsTol not yet supported.')
            else
                if isscalar(AbsTol)
                    % AbsTol is starting scalar
                    for i = 1:nCon
                        if UseModelSeeds
                            % Number of parameters is constant across conditions
                            inT = nT;
                        else
                            % Number of IC parameters may vary across conditions
                            inT = sum(UseSeeds(:,i)) + nTk;
                        end
                        abstol{i} = zeros(nx+integrateObj(i)+nx*inT+integrateObj(i)*inT+nx*nx+nx*inT*nx,1) + AbsTol;
                    end
                elseif all(numel(AbsTol) == nx+integrateObj+nx*nT+nT*integrateObj+nx*nx+nx*nT*nx)
                    % AbsTol is the correct length for all experiments
                    abstol = repmat({AbsTol}, nCon,1);
                else
                    error('KroneckerBio:AbsTol:InvalidAbsTolLength', 'That is not a valid length for AbsTol.')
                end
            end
        end
    otherwise
        error('An unsupported order was passed to fixAbsTol.')
end
