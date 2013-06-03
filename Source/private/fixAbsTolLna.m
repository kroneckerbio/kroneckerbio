function abstol = fixAbsTolLna(absTol, order, integrateObj, nx, nCon)
%FixAbstolLna Standardize the presentation of AbsTol for the linear noise
%   approximation
%
%   abstol = fixAbsTolLna(absTol, order, integrateObj, nx, nCon, useAdjoint,
%   useModelICs, useModelInputs, useParams, useICs, useControls, 
%   selfSensitivities)
%
%   There are many ways to present the absolute integration tolerance to
%   KroneckerBio. This function is the processing center for these
%   different presentations. The many inputs define the type of problem
%   that the AbsTol is needed for; essentially, the length. The standard
%   presentation is a cell array of vectors of the correct length. Each
%   cell corresponds to an experimental condition.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% The default AbsTol
absTolDefault = 1e-9;
if isempty(absTol) || (~iscell(absTol) && ~isstruct(absTol) && all(vec(isnan(absTol))))
    absTol = absTolDefault;
end

% Extract struct inside cell
if iscell(absTol) && isstruct(absTol{1})
    assert(numel(absTol) == 1, 'KroneckerBio:AbsTolLna:TooManyStructsInCells', 'AbsTol was a cell with a struct inside, but there is more than one cell, which is not allowed.')
    absTol = absTol{1};
end

% Extract cell inside cell
if iscell(absTol) && iscell(absTol{1})
    assert(numel(absTol) == 1, 'KroneckerBio:AbsTolLna:TooManyCellsInCells', 'AbsTol was a cell array with a cell array inside, but there is more than one cell in the higher array, which is not allowed.')
    absTol = absTol{1};
end

% Constants
nV = numel(upperInd(nx));
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
                    assert(isfield(temp, 'Lna'), 'KroneckerBio:AbsTolLna:MissingStructField', 'AbsTol is a struct but experiment %i requires a "Lna" field, which does not exist on the struct', i)
                    if ~iscell(temp.System)
                        % It is numeric, copy it to every experiment
                        absTol{i} = temp.System;
                    elseif numel(temp.System) == 1
                        % Only one cell, copy it to every experiment
                        absTol(i) = temp.System;
                    else
                        % Multiple cells, extract the correct one
                        assert(numel(temp.System) >= i, 'KroneckerBio:AbsTolLna:CellVectorTooShort', 'AbsTol is a struct and a cell vector is provided for "Lna" which is required for experiment %i, but the cell vector is not long enough to provide for this experiment', i)
                        absTol(i) = temp.System(i);
                    end
                else
                    % ObjectiveContinuous
                    assert(isfield(temp, 'LnaObjectiveContinuous'), 'KroneckerBio:AbsTolLna:MissingStructField', 'AbsTol is a struct but experiment %i requires a "LnaObjectiveContinuous" field, which does not exist on the struct', i)
                    if ~iscell(temp.ObjectiveContinuous)
                        % It is numeric, copy it to every experiment
                        absTol{i} = temp.ObjectiveContinuous;
                    elseif numel(temp.ObjectiveContinuous) == 1
                        % Only one cell, copy it to every experiment
                        absTol(i) = temp.ObjectiveContinuous;
                    else
                        % Multiple cells, extract the correct one
                        assert(numel(temp.ObjectiveContinuous) >= i, 'KroneckerBio:AbsTolLna:CellVectorTooShort', 'AbsTol is a struct and a cell vector is provided for "LnaObjectiveContinuous" which is required for experiment %i, but the cell vector is not long enough to provide for this experiment', i)
                        absTol(i) = temp.ObjectiveContinuous(i);
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
                abstol{i} = zeros(nx+nV+integrateObj(i),1) + absTol{i};
            elseif numel(absTol{i}) == nx+nV
                % AbsTol is not provided for continuous objective,
                % better hope it's not needed
                assert(~integrateObj(i), 'KroneckerBio:AbsTolLna:VectorWithContinuousObjective', 'Failed to specify AbsTol for continuous objective')
                abstol{i} = absTol{i};
            elseif numel(absTol{i}) == nx+nV+1
                % AbsTol is provided for 1 continuous objective,
                % use it as many times as needed
                abstol{i} = [absTol{i}(1:nx+nV);
                    repmat(absTol{i}(nx+nV+1), integrateObj(i),1)];
            elseif numel(absTol{i}) == nx+nV+integrateObj(i)
                % AbsTol is the correct length
                abstol = repmat({absTol(1:nx+nV+integrateObj)}, nCon,1);
            elseif numel(absTol{i}) >= nx+nV
                % AbsTol is not provided for continuous objective,
                % better hope it's not needed
                assert(~integrateObj(i), 'KroneckerBio:AbsTolLna:VectorWithContinuousObjective', 'Failed to specify AbsTol for continuous objective')
                abstol{i} = absTol{i}(1:nx+nV);
            else
                error('KroneckerBio:AbsTolLna:InvalidAbsTolLength', 'That is not a valid length for AbsTol')
            end
        end
    case 2
        % Gradient integration
        error('An unsupported order was passed to fixAbsTolLna.')
    case 3
        % Hessian integration
        error('An unsupported order was passed to fixAbsTolLna.')
    otherwise
        error('An unsupported order was passed to fixAbsTolLna.')
end
