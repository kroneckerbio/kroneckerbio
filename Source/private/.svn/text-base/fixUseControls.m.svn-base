function [useControls, nTq] = fixUseControls(useControls, useModelInputs, nCon, nqModel, nqCon)
%[useControls, nTq] = fixUseControls(useControls, useModelInputs, nCon,
%nqModel, nqCon)
if useModelInputs
    if iscell(useControls)
        assert(numel(useControls)==1, 'KroneckerBio:UseControls:CellLength', 'If UseModelInputs is true, then UseControls cannot be a cell array of length > 1')
        useControls = useControls{1};
    end
    
    useControls = vec(useControls);
    
    if isnumeric(useControls)
        assert(all(floor(useControls) == useControls) && all(useControls >= 1), 'KroneckerBio:UseControls:InvalidValue', 'UseControls is an invalid linear index')
        assert(all(useControls <= nqModel), 'KroneckerBio:UseControls:LinearIndexOutOfRange', 'If UseModelInputs is true and UseControls is provided as a linear index, then no linear index can be larger than m.nq')
        temp = false(nqModel,1);
        temp(useControls) = true;
        useControls = temp;
    elseif islogical(useControls)
        assert(numel(useControls) <= nqModel, 'KroneckerBio:UseControls:InvalidLogicalLength', 'If UseModelInputs is true and UseControls is provided as a logical index, numel(UseControls) cannot be larger than m.nq')
        useControls = [useControls; false(nqModel - numel(useControls),1)];
    else
        error('KroneckerBio:UseControls:InvalidType', 'UseControls must be provided as logical or linear index vector into m.q')
    end
    
    nTq = nnz(useControls);
    useControls = {useControls};
else%~useModelInputs
    if ~iscell(useControls)
        % Repeat across experiments
        useControls = repmat({useControls}, nCon,1);
    end
    
    useControls = vec(useControls);
    
    % Repeat a singular array
    if numel(useControls) ~= nCon
        assert(numel(useControls) == 1 , 'KroneckerBio:UseControls:CellLength', 'If UseModelInputs is false and UseControls is provided as a cell array, its length must be equal to the length of con or 1')
        useControls = repmat(useControls, nCon,1);
    end
    
    nTq = 0;
    for iCon = 1:nCon
        if isnumeric(useControls{iCon})
            assert(all(floor(useControls{iCon}) == useControls{iCon}) && all(useControls{iCon} >= 1), 'KroneckerBio:UseControls:InvalidValue', 'UseControls is an invalid linear index')
            assert(all(useControls{iCon} <= nqCon(iCon)), 'KroneckerBio:UseControls:LinearIndexOutOfRange', 'If UseModelInputs is false and UseControls is provided as a linear index, then no linear index can be larger than the smallest con.nq')
            temp = false(nqCon(iCon),1);
            temp(useControls{iCon}) = true;
            useControls{iCon} = temp;
        elseif islogical(useControls{iCon})
            assert(numel(useControls{iCon}) <= nqCon(iCon), 'KroneckerBio:UseControls:InvalidLogicalLength', 'If UseModelInputs is false and UseControls is provided as a logical index, numel(UseControls) cannot be larger than the smallest con.nq')
            useControls{iCon} = [vec(useControls{iCon}); false(nqModel - numel(useControls{iCon}),1)];
        else
            error('KroneckerBio:UseControls:InvalidType', 'UseControls must be provided as logical or linear index vector into con.q')
        end
        nTq = nTq + nnz(useControls{iCon});
    end
end