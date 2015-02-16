function [UseControls, nTq] = fixUseControls(UseControls, nCon, nq)
%[useControls, nTq] = fixUseControls(useControls, nCon, nq)

if ~iscell(UseControls)
    % Repeat across experiments
    UseControls = repmat({UseControls}, nCon,1);
end

UseControls = vec(UseControls);

% Repeat a singular array
if numel(UseControls) ~= nCon
    assert(numel(UseControls) == 1 , 'KroneckerBio:UseControls:CellLength', 'UseControls, when a cell array, must have a length equal to nCon or 1')
    UseControls = repmat(UseControls, nCon,1);
end

nTq = 0;
for iCon = 1:nCon
    if isnumeric(UseControls{iCon}) && isscalar(UseControls{iCon}) && isnan(UseControls{iCon})
        % Default
        UseControls{iCon} = true(nq(iCon),1);
    elseif isnumeric(UseControls{iCon}) && all(floor(vec(UseControls{iCon})) == vec(UseControls{iCon})) && all(UseControls{iCon} >= 1)
        % Linear index
        UseControls{iCon} = vec(UseControls{iCon});
        assert(all(UseControls{iCon} <= nq(iCon)), 'KroneckerBio:UseControls:LinearIndexOutOfRange', 'UseControls, when a linear index, can have no index larger than con.nq')
        temp = false(nq(iCon),1);
        temp(UseControls{iCon}) = true;
        UseControls{iCon} = temp;
    elseif islogical(UseControls{iCon})
        assert(numel(UseControls{iCon}) <= nq(iCon), 'KroneckerBio:UseControls:InvalidLogicalLength', 'UseControls, when a logical index, cannot be longer than con.nq')
        UseControls{iCon} = [vec(UseControls{iCon}); false(nq - numel(UseControls{iCon}),1)];
    else
        error('KroneckerBio:UseControls:InvalidType', 'UseControls must be provided as logical or linear index vector into con.q')
    end
    nTq = nTq + nnz(UseControls{iCon});
end
