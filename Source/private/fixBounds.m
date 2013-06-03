function bounds = fixBounds(bounds, UseModelSeeds, UseModelInputs, UseParams, UseSeeds, UseControls)

% Constants
ns = size(UseSeeds, 1);
nk = numel(UseParams);
nq = numel(cat(1,UseControls{:}));
nTk = nnz(UseParams);
nTs = nnz(UseSeeds);
nTq = nnz(cat(1,UseControls{:}));
nT = nTk + nTs + nTq;
nCon = size(UseSeeds, 2);

if UseModelSeeds && UseModelInputs
    % Bounds can be nk+ns+nq, nT, nk, nk+ns, or 1
    l = numel(bounds);
    
    if l == 1
        bounds = zeros(nT,1) + bounds;
    elseif l == nk+ns+nq
        bounds = bounds([UseParams; UseSeeds; cat(1,UseControls{:})]);
    elseif l == nT
        %bounds = bounds;
    elseif l == nk && nTs == 0 && nTq == 0
        bounds = bounds(UseParams);
    elseif l == nk+ns && nTq == 0
        bounds = bounds([UseParams; UseSeeds]);
    else
        error('KroneckerBio:BoundSize', ...
            'LowerBound and UpperBound must be vectors the length of m.nk+m.ns+m.nq, number of variable parameters, m.nk if there are no variable ICs or controls, m.nk+m.ns if there are no variable controls, or scalar')
    end
elseif ~UseModelSeeds && UseModelInputs
    % Bounds can be nk+ns*nCon+nq, nT, nk, nk+ns*nCon, or 1
    l = numel(bounds);
    
    if l == 1
        bounds = zeros(nT,1) + bounds;
    elseif l == nk+ns*nCon+nq
        bounds = bounds([UseParams; vec(UseSeeds); cat(1,UseControls{:})]);
    elseif l == nT
        %bounds = bounds;
    elseif l == nk && nTs == 0 && nTq == 0
        bounds = bounds(UseParams);
    elseif l == nk+ns*nCon && nTq == 0
        bounds = bounds([UseParams; vec(UseSeeds)]);
    else
        error('KroneckerBio:BoundSize', ...
            'LowerBound and UpperBound must be vectors the length of m.nk+m.ns*numel(con)+m.nq, number of variable parameters, m.nk if there are no variable ICs or controls, m.nk+m.ns*numel(con) if there are no variable controls, or scalar')
    end
elseif UseModelSeeds && ~UseModelInputs
    % Bounds can be nk+ns+nq, nT, nk, nk+ns, or 1
    l = numel(bounds);
    
    if l == 1
        bounds = zeros(nT,1) + bounds;
    elseif l == nk+ns+nq
        bounds = bounds([UseParams; UseSeeds; cat(1,UseControls{:})]);
    elseif l == nT
        %bounds = bounds;
    elseif l == nk && nTs == 0 && nTq == 0
        bounds = bounds(UseParams);
    elseif l == nk+ns && nTq == 0
        bounds = bounds([UseParams; UseSeeds]);
    else
        error('KroneckerBio:BoundSize', ...
            'LowerBound and UpperBound must be vectors the length of m.nk+m.ns+m.nq, number of variable parameters, m.nk if there are no variable ICs or controls, m.nk+m.ns if there are no variable controls, or scalar')
    end
else%~useModelICs && ~useModelInputs
    % Bounds can be nk+ns*nCon+nq, nT, nk, nk+ns*nCon, or 1
    l = numel(bounds);
    
    if l == 1
        bounds = zeros(nT,1) + bounds;
    elseif l == nk+ns*nCon+nq
        bounds = bounds([UseParams; vec(UseSeeds); cat(1,UseControls{:})]);
    elseif l == nT
        %bounds = bounds;
    elseif l == nk && nTs == 0 && nTq == 0
        bounds = bounds(UseParams);
    elseif l == nk+ns*nCon && nTq == 0
        bounds = bounds([UseParams; vec(UseSeeds)]);
    else
        error('KroneckerBio:BoundSize', ...
            'LowerBound and UpperBound must be vectors the length of m.nk+m.ns*numel(con)+m.nq, number of variable parameters, m.nk if there are no variable ICs or controls, m.nk+m.ns*numel(con) if there are no variable controls, or scalar')
    end
end