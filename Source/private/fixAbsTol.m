function out = fixAbsTol(AbsTol, order, integrateObj, nx, nCon, UseAdjoint, UseParams, UseSeeds, UseInputControls, UseDoseControls)
%fixAbsTol Standardize the presentation of AbsTol
%
%   abstol = fixAbsTol(absTol, order, integrateObj, ns, nCon, UseAdjoint,
%   UseParams, UseSeeds, UseInputControls, UseDoseControls)
%
%   There are many ways to present the absolute integration tolerance to
%   KroneckerBio. This function is the processing center for these
%   different presentations. The many inputs define the type of problem
%   that the AbsTol is needed for; essentially, the length. The standard
%   presentation is a cell array of vectors of the correct length. Each
%   cell corresponds to an experimental condition.

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean up inputs
if nargin < 10
    UseDoseControls = {};
    if nargin < 9
        UseInputControls = {};
        if nargin < 8
            UseSeeds = [];
            if nargin < 7
                UseParams = [];
                if nargin < 6
                    UseAdjoint = [];
                end
            end
        end
    end
end

% The default AbsTol
absTolDefault = 1e-9;
if isempty(AbsTol)
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
nTq = sum(cat(1, UseInputControls{:}));
nTh = sum(cat(1, UseDoseControls{:}));
nT = nTk + nTs + nTq + nTh;

out = cell(nCon,1);

switch order
    case 1
        % System or Objective integration
        if isstruct(AbsTol)
            % It is a struct, extract the appropriate AbsTol
            temp = AbsTol;
            AbsTol = cell(nCon,1);
            for i = 1:nCon
                if ~integrateObj(i)
                    % System
                    assert(isfield(temp, 'System'), 'KroneckerBio:AbsTol:MissingStructField', 'AbsTol is a struct but experiment %i requires a "System" field, which does not exist on the struct', i)
                    if isscalar(temp)
                        AbsTol{i} = temp.System;
                    else
                        AbsTol{i} = temp(i).System;
                    end
                else
                    % Objective
                    assert(isfield(temp, 'Objective'), 'KroneckerBio:AbsTol:MissingStructField', 'AbsTol is a struct but experiment %i requires a "Objective" field, which does not exist on the struct', i)
                    if isscalar(temp)
                        AbsTol{i} = temp.Objective;
                    else
                        AbsTol{i} = temp(i).Objective;
                    end
                end
            end
        end
        
        if ~iscell(AbsTol)
            % It is numeric, copy it to every experiment
            AbsTol = repmat({AbsTol}, nCon,1);
        elseif numel(AbsTol) == 1
            % Only one cell, copy it to every experiment
            AbsTol = repmat(AbsTol, nCon,1);
        end
        
        for i = 1:nCon
            % AbsTol is the same for all conditions
            if isscalar(AbsTol{i})
                % Copy the value to all species and conditions
                out{i} = zeros(nx+integrateObj(i),1) + AbsTol{i};
            elseif numel(AbsTol{i}) == nx
                % AbsTol is not provided for continuous objective,
                % better hope it's not needed
                assert(~integrateObj(i), 'KroneckerBio:AbsTol:VectorWithContinuousObjective', 'Failed to specify AbsTol for continuous objective')
                out{i} = AbsTol{i};
            elseif numel(AbsTol{i}) == nx+1
                % AbsTol is provided for 1 continuous objective,
                % use it as many times as needed
                out{i} = [AbsTol{i}(1:nx);
                    repmat(AbsTol{i}(nx+1), integrateObj(i),1)];
            elseif numel(AbsTol{i}) == nx+integrateObj(i)
                % AbsTol is the correct length
                out = repmat({AbsTol(1:nx+integrateObj)}, nCon,1);
            elseif numel(AbsTol{i}) >= nx
                % AbsTol is not provided for continuous objective,
                % better hope it's not needed
                assert(~integrateObj(i), 'KroneckerBio:AbsTol:VectorWithContinuousObjective', 'Failed to specify AbsTol for continuous objective')
                out{i} = AbsTol{i}(1:nx);
            else
                error('KroneckerBio:AbsTol:InvalidAbsTolLength', 'That is not a valid length for AbsTol')
            end
        end
    case 2
        % Sensitivity or Gradient or Adjoint integration
        if isstruct(AbsTol)
            % It is a struct, extract the appropriate AbsTol
            temp = AbsTol;
            AbsTol = cell(nCon,1);
            for i = 1:nCon
                if ~UseAdjoint
                    if ~integrateObj(i)
                        % Sensitivity
                        assert(isfield(temp, 'Sensitivity'), 'KroneckerBio:AbsTol:MissingStructField', 'AbsTol is a struct but experiment %i requires a "Sensitivity" field, which does not exist on the struct', i)
                        if isscalar(temp)
                            AbsTol{i} = temp.Sensitivity;
                        else
                            AbsTol{i} = temp(i).Sensitivity;
                        end
                    else
                        % GradientContinuous
                        assert(isfield(temp, 'Gradient'), 'KroneckerBio:AbsTol:MissingStructField', 'AbsTol is a struct but experiment %i requires a "Gradient" field, which does not exist on the struct', i)
                        if isscalar(temp)
                            AbsTol{i} = temp.Gradient;
                        else
                            AbsTol{i} = temp(i).Gradient;
                        end
                    end
                else %useAdjoint
                    % Adjoint
                    assert(isfield(temp, 'Adjoint'), 'KroneckerBio:AbsTol:MissingStructField', 'AbsTol is a struct but experiment %i requires a "Adjoint" field, which does not exist on the struct', i)
                    if isscalar(temp)
                        AbsTol{i} = temp.Adjoint;
                    else
                        AbsTol{i} = temp(i).Adjoint;
                    end
                end
            end
        end
        
        if ~iscell(AbsTol)
            % It is numeric, copy it to every experiment
            AbsTol = repmat({AbsTol}, nCon,1);
        elseif numel(AbsTol) == 1
            % Only one cell, copy it to every experiment
            AbsTol = repmat(AbsTol, nCon,1);
        end
        
        % Process each vector in each cell
        assert(numel(AbsTol) >= nCon, 'KroneckerBio:AbsTol:CellVectorTooShort', 'AbsTol was provided as a cell vector of length %i, but the cell vector is too short for the number of experiments %i', numel(AbsTol), nCon)
        for i = 1:nCon
            inTs = nnz(UseSeeds(:,i));
            inTq = nnz(UseInputControls{i});
            inTh = nnz(UseDoseControls{i});
            
            inT = nTk + inTs + inTq + inTh;
            
            if UseAdjoint
                if isscalar(AbsTol{i})
                    % Copy the value to all species and conditions
                    out{i} = zeros(nx+integrateObj(i)+nx+nT,1) + AbsTol{i};
                else %isvector
                    % Verify vector length
                    assert(numel(AbsTol{i}) == nx+integrateObj(i)+nx+nT, 'KroneckerBio:AbsTol:InvalidAbsTolLength', 'That is not a valid length for AbsTol')
                    out{i} = AbsTol{i};
                end
            else %~useAdjoint
                if isscalar(AbsTol{i})
                    % AbsTol is starting scalar
                    out{i} = zeros(nx+integrateObj(i)+nx*inT+integrateObj(i)*inT,1) + AbsTol{i};
                elseif numel(AbsTol{i}) == nx+nx*inT
                    % AbsTol is not provided for continuous objective,
                    % better hope it's not needed
                    assert(~integrateObj(i), 'KroneckerBio:AbsTol:VectorWithContinuousObjective', 'Failed to specify AbsTol for continuous objective and gradient')
                    out{i} = AbsTol{i};
                elseif numel(AbsTol{i}) == nx+1+nx*inT+inT
                    % AbsTol is provided for 1 continuous objective,
                    % use it as many times as needed
                    out{i} = [AbsTol{i}(1:nx);
                        repmat(AbsTol{i}(nx+1), integrateObj(i),1);
                        AbsTol{i}(nx+1+1:nx+1+nx*inT);
                        repmat(AbsTol{i}(nx+1+nx*inT+1:nx+1+nx*inT+inT), integrateObj(i),1)];
                elseif numel(AbsTol) == nx+integrateObj(i)+nx*inT+inT*integrateObj(i)
                    % AbsTol is the correct length
                    out{i} = AbsTol{i};
                else
                    error('KroneckerBio:AbsTol:InvalidAbsTolLength', 'That is not a valid length for AbsTol.')
                end
            end
        end
    case 3
        % Curvature or Hessian integration
        if isstruct(AbsTol)
            % It is a struct, extract the appropriate AbsTol
            temp = AbsTol;
            AbsTol = cell(nCon,1);
            for i = 1:nCon
                if ~integrateObj(i)
                    % Sensitivity
                    assert(isfield(temp, 'Curvature'), 'KroneckerBio:AbsTol:MissingStructField', 'AbsTol is a struct but experiment %i requires a "Doublesensitivity" field, which does not exist on the struct', i)
                    if ~iscell(temp.Curvature)
                        % It is numeric, copy it to every experiment
                        AbsTol{i} = temp.Curvature;
                    elseif numel(temp.Curvature) == 1
                        % Only one cell, copy it to every experiment
                        AbsTol(i) = temp.Curvature;
                    else
                        % Multiple cells, extract the correct one
                        assert(numel(temp.Curvature) >= i, 'KroneckerBio:AbsTol:CellVectorTooShort', 'AbsTol is a struct and a cell vector is provided for "Doublesensitivity" which is required for experiment %i, but the cell vector is not long enough to provide for this experiment', i)
                        AbsTol(i) = temp.Curvature(i);
                    end
                else
                    % HessianContinuous
                    assert(isfield(temp, 'Hessian'), 'KroneckerBio:AbsTol:MissingStructField', 'AbsTol is a struct but experiment %i requires a "HessianContinuous" field, which does not exist on the struct', i)
                    if ~iscell(temp.Hessian)
                        % It is numeric, copy it to every experiment
                        AbsTol{i} = temp.Hessian;
                    elseif numel(temp.Hessian) == 1
                        % Only one cell, copy it to every experiment
                        AbsTol(i) = temp.Hessian;
                    else
                        % Multiple cells, extract the correct one
                        assert(numel(temp.Hessian) >= i, 'KroneckerBio:AbsTol:CellVectorTooShort', 'AbsTol is a struct and a cell vector is provided for "HessianContinuous" which is required for experiment %i, but the cell vector is not long enough to provide for this experiment', i)
                        AbsTol(i) = temp.Hessian(i);
                    end
                end
            end
        end
        
        if ~iscell(AbsTol)
            % It is numeric, copy it to every experiment
            AbsTol = repmat({AbsTol}, nCon,1);
        elseif numel(AbsTol) == 1
            % Only one cell, copy it to every experiment
            AbsTol = repmat(AbsTol, nCon,1);
        end
        
        % Process each vector in each cell
        assert(numel(AbsTol) >= nCon, 'KroneckerBio:AbsTol:CellVectorTooShort', 'AbsTol was provided as a cell vector of length %i, but the cell vector is too short for the number of experiments %i', numel(AbsTol), nCon)
        for i = 1:nCon
            inTs = nnz(UseSeeds(:,i));
            inTq = nnz(UseInputControls{i});
            inTh = nnz(UseDoseControls{i});
            
            inT = nTk + inTs + inTq + inTh;
            
            if isscalar(AbsTol{i})
                % AbsTol is starting scalar
                out{i} = zeros(nx+integrateObj(i)+nx*inT+integrateObj(i)*inT+nx*inT*inT+integrateObj(i)*inT*inT,1) + AbsTol{i};
            elseif numel(AbsTol{i}) == nx+nx*inT+nx*inT*inT
                % AbsTol is not provided for continuous objective, better hope it's not needed
                assert(~integrateObj(i), 'KroneckerBio:AbsTol:VectorWithContinuousObjective', 'Failed to specify AbsTol for continuous objective and gradient')
                out{i} = AbsTol{i};
            elseif numel(AbsTol{i}) == nx+1+nx*inT+inT+nx*inT*inT+inT*inT
                % AbsTol is provided for 1 continuous objective,
                % use it as many times as needed
                out{i} = [AbsTol{i}(1:nx);
                    repmat(AbsTol{i}(nx+1), integrateObj(i),1);
                    AbsTol{i}(nx+1+1:nx+1+nx*inT);
                    repmat(AbsTol{i}(nx+1+nx*inT+1:nx+1+nx*inT+inT), integrateObj(i),1);
                    AbsTol{i}(nx+1+nx*inT+inT+1:nx+1+nx*inT+inT+nx*inT*inT);
                    repmat(AbsTol(AbsTol{i}(nx+1+nx*inT+inT+nx*inT*inT+1:nx+1+nx*inT+inT+nx*inT*inT+inT*inT)), integrateObj(i),1)];
            elseif numel(AbsTol) == nx+integrateObj(i)+nx*nT+nT*integrateObj(i)+nx*inT*inT+inT*inT*integrateObj(i)
                % AbsTol is the correct length
                out{i} = AbsTol{i};
            else
                error('KroneckerBio:AbsTol:InvalidAbsTolLength', 'That is not a valid length for AbsTol.')
            end
        end
    otherwise
        error('An unsupported order was passed to fixAbsTol.')
end
