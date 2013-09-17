function abstol = fixAbsTol(absTol, order, integrateObj, nx, nCon, useAdjoint, UseModelSeeds, useModelInputs, useParams, UseSeeds, useControls)
%fixAbsTol Standardize the presentation of AbsTol
%
%   abstol = fixAbsTol(absTol, order, integrateObj, ns, nCon, useAdjoint,
%   UseModelSeeds, useModelInputs, useParams, UseSeeds, useControls)
%
%   There are many ways to present the absolute integration tolerance to
%   KroneckerBio. This function is the processing center for these
%   different presentations. The many inputs define the type of problem
%   that the AbsTol is needed for; essentially, the length. The standard
%   presentation is a cell array of vectors of the correct length. Each
%   cell corresponds to an experimental condition.

% (c) 2012 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean up inputs
if nargin < 11
    useControls = {};
    if nargin < 10
        UseSeeds = [];
        if nargin < 9
            useParams = [];
            if nargin < 8
                UseModelSeeds = [];
                if nargin < 7
                    useModelInputs = [];
                    if nargin < 6
                        useAdjoint = [];
                    end
                end
            end
        end
    end
end

% The default AbsTol
absTolDefault = 1e-9;
if isempty(absTol) || (~iscell(absTol) && ~isstruct(absTol) && all(vec(isnan(absTol))))
    absTol = absTolDefault;
end

% Extract struct inside cell
if iscell(absTol) && isstruct(absTol{1})
    assert(numel(absTol) == 1, 'KroneckerBio:AbsTol:TooManyStructsInCells', 'AbsTol was a cell with a struct inside, but there is more than one cell, which is not allowed.')
    absTol = absTol{1};
end

% Extract cell inside cell
if iscell(absTol) && iscell(absTol{1})
    assert(numel(absTol) == 1, 'KroneckerBio:AbsTol:TooManyCellsInCells', 'AbsTol was a cell array with a cell array inside, but there is more than one cell in the higher array, which is not allowed.')
    absTol = absTol{1};
end

% Constants
nTk = sum(useParams);
nTs = sum(sum(UseSeeds));
nTq = sum(cat(1, useControls{:}));
nT = nTk + nTs + nTq;

abstol = cell(nCon,1);

switch order
    case 1
        % System or ObjectiveContinuous integration
        if isstruct(absTol)
            % It is a struct, extract the appropriate AbsTol
            temp = absTol;
            absTol = cell(nCon,1);
            for i = 1:nCon
                if ~integrateObj(i)
                    % System
                    assert(isfield(temp, 'System'), 'KroneckerBio:AbsTol:MissingStructField', 'AbsTol is a struct but experiment %i requires a "System" field, which does not exist on the struct', i)
                    if isscalar(temp)
                        absTol{i} = temp.System;
                    else
                        absTol{i} = temp(i).System;
                    end
                else
                    % ObjectiveContinuous
                    assert(isfield(temp, 'Objective'), 'KroneckerBio:AbsTol:MissingStructField', 'AbsTol is a struct but experiment %i requires a "Objective" field, which does not exist on the struct', i)
                    if isscalar(temp)
                        absTol{i} = temp.Objective;
                    else
                        absTol{i} = temp(i).Objective;
                    end
                end
            end
        end
        
        if ~iscell(absTol)
            % It is numeric, copy it to every experiment
            absTol = repmat({absTol}, nCon,1);
        elseif numel(absTol) == 1
            % Only one cell, copy it to every experiment
            absTol = repmat(absTol, nCon,1);
        end
        
        for i = 1:nCon
            % AbsTol is the same for all conditions
            if isscalar(absTol{i})
                % Copy the value to all species and conditions
                abstol{i} = zeros(nx+integrateObj(i),1) + absTol{i};
            elseif numel(absTol{i}) == nx
                % AbsTol is not provided for continuous objective,
                % better hope it's not needed
                assert(~integrateObj(i), 'KroneckerBio:AbsTol:VectorWithContinuousObjective', 'Failed to specify AbsTol for continuous objective')
                abstol{i} = absTol{i};
            elseif numel(absTol{i}) == nx+1
                % AbsTol is provided for 1 continuous objective,
                % use it as many times as needed
                abstol{i} = [absTol{i}(1:nx);
                    repmat(absTol{i}(nx+1), integrateObj(i),1)];
            elseif numel(absTol{i}) == nx+integrateObj(i)
                % AbsTol is the correct length
                abstol = repmat({absTol(1:nx+integrateObj)}, nCon,1);
            elseif numel(absTol{i}) >= nx
                % AbsTol is not provided for continuous objective,
                % better hope it's not needed
                assert(~integrateObj(i), 'KroneckerBio:AbsTol:VectorWithContinuousObjective', 'Failed to specify AbsTol for continuous objective')
                abstol{i} = absTol{i}(1:nx);
            else
                error('KroneckerBio:AbsTol:InvalidAbsTolLength', 'That is not a valid length for AbsTol')
            end
        end
    case 2
        % Sensitivity or GradientContinuous integration
        if isstruct(absTol)
            % It is a struct, extract the appropriate AbsTol
            temp = absTol;
            absTol = cell(nCon,1);
            for i = 1:nCon
                if ~useAdjoint
                    if ~integrateObj(i)
                        % Sensitivity
                        assert(isfield(temp, 'Sensitivity'), 'KroneckerBio:AbsTol:MissingStructField', 'AbsTol is a struct but experiment %i requires a "Sensitivity" field, which does not exist on the struct', i)
                        if isscalar(temp)
                            absTol{i} = temp.Sensitivity;
                        else
                            absTol{i} = temp(i).Sensitivity;
                        end
                    else
                        % GradientContinuous
                        assert(isfield(temp, 'Gradient'), 'KroneckerBio:AbsTol:MissingStructField', 'AbsTol is a struct but experiment %i requires a "Gradient" field, which does not exist on the struct', i)
                        if isscalar(temp)
                            absTol{i} = temp.Gradient;
                        else
                            absTol{i} = temp(i).Gradient;
                        end
                    end
                else %useAdjoint
                    % Adjoint
                    assert(isfield(temp, 'Adjoint'), 'KroneckerBio:AbsTol:MissingStructField', 'AbsTol is a struct but experiment %i requires a "Adjoint" field, which does not exist on the struct', i)
                    if isscalar(temp)
                        absTol{i} = temp.Adjoint;
                    else
                        absTol{i} = temp(i).Adjoint;
                    end
                end
            end
        end
        
        if ~iscell(absTol)
            % It is numeric, copy it to every experiment
            absTol = repmat({absTol}, nCon,1);
        elseif numel(absTol) == 1
            % Only one cell, copy it to every experiment
            absTol = repmat(absTol, nCon,1);
        end
        
        % Process each vector in each cell
        assert(numel(absTol) >= nCon, 'KroneckerBio:AbsTol:CellVectorTooShort', 'AbsTol was provided as a cell vector of length %i, but the cell vector is too short for the number of experiments %i', numel(absTol), nCon)
        for i = 1:nCon
            % If opts.UseModelSeeds is false, the number of variables can change
            if UseModelSeeds
                inTs = nTs;
            else
                inTs = nnz(UseSeeds(:,i));
            end
            
            % If opts.UseModelInputs is false, the number of variables can change
            if useModelInputs
                inTq = nTq;
            else
                inTq = nnz(useControls{i});
            end
            
            inT = nTk + inTs + inTq;
            
            if useAdjoint
                if isscalar(absTol{i})
                    % Copy the value to all species and conditions
                    abstol{i} = zeros(nx+integrateObj(i)+nx+nT,1) + absTol{i};
                else %isvector
                    % Verify vector length
                    assert(numel(absTol{i}) == nx+integrateObj(i)+nx+nT, 'KroneckerBio:AbsTol:InvalidAbsTolLength', 'That is not a valid length for AbsTol')
                    abstol{i} = absTol{i};
                end
            else %~useAdjoint
                if isscalar(absTol{i})
                    % AbsTol is starting scalar
                    abstol{i} = zeros(nx+integrateObj(i)+nx*inT+integrateObj(i)*inT,1) + absTol{i};
                elseif numel(absTol{i}) == nx+nx*inT
                    % AbsTol is not provided for continuous objective,
                    % better hope it's not needed
                    assert(~integrateObj(i), 'KroneckerBio:AbsTol:VectorWithContinuousObjective', 'Failed to specify AbsTol for continuous objective and gradient')
                    abstol{i} = absTol{i};
                elseif numel(absTol{i}) == nx+1+nx*inT+inT
                    % AbsTol is provided for 1 continuous objective,
                    % use it as many times as needed
                    abstol{i} = [absTol{i}(1:nx);
                        repmat(absTol{i}(nx+1), integrateObj(i),1);
                        absTol{i}(nx+1+1:nx+1+nx*inT);
                        repmat(absTol{i}(nx+1+nx*inT+1:nx+1+nx*inT+inT), integrateObj(i),1)];
                elseif numel(absTol) == nx+integrateObj(i)+nx*inT+inT*integrateObj(i)
                    % AbsTol is the correct length
                    abstol{i} = absTol{i};
                else
                    error('KroneckerBio:AbsTol:InvalidAbsTolLength', 'That is not a valid length for AbsTol.')
                end
            end
        end
    case 3
        % Doublesensitivity or HessianContinuous integration
        if isstruct(absTol)
            % It is a struct, extract the appropriate AbsTol
            temp = absTol;
            absTol = cell(nCon,1);
            for i = 1:nCon
                if ~integrateObj(i)
                    % Sensitivity
                    assert(isfield(temp, 'Doublesensitivity'), 'KroneckerBio:AbsTol:MissingStructField', 'AbsTol is a struct but experiment %i requires a "Doublesensitivity" field, which does not exist on the struct', i)
                    if ~iscell(temp.Doublesensitivity)
                        % It is numeric, copy it to every experiment
                        absTol{i} = temp.Doublesensitivity;
                    elseif numel(temp.Doublesensitivity) == 1
                        % Only one cell, copy it to every experiment
                        absTol(i) = temp.Doublesensitivity;
                    else
                        % Multiple cells, extract the correct one
                        assert(numel(temp.Doublesensitivity) >= i, 'KroneckerBio:AbsTol:CellVectorTooShort', 'AbsTol is a struct and a cell vector is provided for "Doublesensitivity" which is required for experiment %i, but the cell vector is not long enough to provide for this experiment', i)
                        absTol(i) = temp.Doublesensitivity(i);
                    end
                else
                    % HessianContinuous
                    assert(isfield(temp, 'HessianContinuous'), 'KroneckerBio:AbsTol:MissingStructField', 'AbsTol is a struct but experiment %i requires a "HessianContinuous" field, which does not exist on the struct', i)
                    if ~iscell(temp.HessianContinuous)
                        % It is numeric, copy it to every experiment
                        absTol{i} = temp.HessianContinuous;
                    elseif numel(temp.HessianContinuous) == 1
                        % Only one cell, copy it to every experiment
                        absTol(i) = temp.HessianContinuous;
                    else
                        % Multiple cells, extract the correct one
                        assert(numel(temp.HessianContinuous) >= i, 'KroneckerBio:AbsTol:CellVectorTooShort', 'AbsTol is a struct and a cell vector is provided for "HessianContinuous" which is required for experiment %i, but the cell vector is not long enough to provide for this experiment', i)
                        absTol(i) = temp.HessianContinuous(i);
                    end
                end
            end
        end
        
        if ~iscell(absTol)
            % It is numeric, copy it to every experiment
            absTol = repmat({absTol}, nCon,1);
        elseif numel(absTol) == 1
            % Only one cell, copy it to every experiment
            absTol = repmat(absTol, nCon,1);
        end
        
        % Process each vector in each cell
        assert(numel(absTol) >= nCon, 'KroneckerBio:AbsTol:CellVectorTooShort', 'AbsTol was provided as a cell vector of length %i, but the cell vector is too short for the number of experiments %i', numel(absTol), nCon)
        for i = 1:nCon
            % If opts.UseModelICs is false, the number of variables can change
            if UseModelSeeds
                inTs = nTs;
            else
                inTs = nnz(UseSeeds(:,i));
            end
            
            % If opts.UseModelInputs is false, the number of variables can change
            if useModelInputs
                inTq = nTq;
            else
                inTq = nnz(useControls{i});
            end
            
            inT = nTk + inTs + inTq;
            
            if isscalar(absTol{i})
                % AbsTol is starting scalar
                abstol{i} = zeros(nx+integrateObj(i)+nx*inT+integrateObj(i)*inT+nx*inT*inT+integrateObj(i)*inT*inT,1) + absTol{i};
            elseif numel(absTol{i}) == nx+nx*inT+nx*inT*inT
                % AbsTol is not provided for continuous objective,
                % better hope it's not needed
                assert(~integrateObj(i), 'KroneckerBio:AbsTol:VectorWithContinuousObjective', 'Failed to specify AbsTol for continuous objective and gradient')
                abstol{i} = absTol{i};
            elseif numel(absTol{i}) == nx+1+nx*inT+inT+nx*inT*inT+inT*inT
                % AbsTol is provided for 1 continuous objective,
                % use it as many times as needed
                abstol{i} = [absTol{i}(1:nx);
                    repmat(absTol{i}(nx+1), integrateObj(i),1);
                    absTol{i}(nx+1+1:nx+1+nx*inT);
                    repmat(absTol{i}(nx+1+nx*inT+1:nx+1+nx*inT+inT), integrateObj(i),1);
                    absTol{i}(nx+1+nx*inT+inT+1:nx+1+nx*inT+inT+nx*inT*inT);
                    repmat(absTol(absTol{i}(nx+1+nx*inT+inT+nx*inT*inT+1:nx+1+nx*inT+inT+nx*inT*inT+inT*inT)), integrateObj(i),1)];
            elseif numel(absTol) == nx+integrateObj(i)+nx*nT+nT*integrateObj(i)+nx*inT*inT+inT*inT*integrateObj(i)
                % AbsTol is the correct length
                abstol{i} = absTol{i};
            else
                error('KroneckerBio:AbsTol:InvalidAbsTolLength', 'That is not a valid length for AbsTol.')
            end
        end
    otherwise
        error('An unsupported order was passed to fixAbsTol.')
end
