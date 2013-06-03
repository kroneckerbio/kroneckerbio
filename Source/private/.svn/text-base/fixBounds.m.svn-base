function bounds = fixBounds(bounds, useModelICs, useModelInputs, useParams, useICs, useControls)

% Constants
nx = size(useICs, 1);
nk = numel(useParams);
nq = numel(cat(1,useControls{:}));
nTk = nnz(useParams);
nTx = nnz(useICs);
nTq = nnz(cat(1,useControls{:}));
nT = nTk + nTx + nTq;
nCon = size(useICs, 2);

if useModelICs && useModelInputs
    % Bounds can be nk+nx+nq, nT, nk, nk+nx, or 1
    l = numel(bounds);
    
    if l == 1
        bounds = zeros(nT,1) + bounds;
    elseif l == nk+nx+nq
        bounds = bounds([useParams; useICs; cat(1,useControls{:})]);
    elseif l == nT
        %bounds = bounds;
    elseif l == nk && nTx == 0 && nTq == 0
        bounds = bounds(useParams);
    elseif l == nk+nx && nTq == 0
        bounds = bounds([useParams; useICs]);
    else
        error('KroneckerBio:BoundSize', ...
            'LowerBound and UpperBound must be vectors the length of m.nk+m.nx+m.nq, number of variable parameters, m.nk if there are no variable ICs or controls, m.nk+m.nx if there are no variable controls, or scalar')
    end
elseif ~useModelICs && useModelInputs
    % Bounds can be nk+nx*nCon+nq, nT, nk, nk+nx*nCon, or 1
    l = numel(bounds);
    
    if l == 1
        bounds = zeros(nT,1) + bounds;
    elseif l == nk+nx*nCon+nq
        bounds = bounds([useParams; vec(useICs); cat(1,useControls{:})]);
    elseif l == nT
        %bounds = bounds;
    elseif l == nk && nTx == 0 && nTq == 0
        bounds = bounds(useParams);
    elseif l == nk+nx*nCon && nTq == 0
        bounds = bounds([useParams; vec(useICs)]);
    else
        error('KroneckerBio:BoundSize', ...
            'LowerBound and UpperBound must be vectors the length of m.nk+m.nx*numel(con)+m.nq, number of variable parameters, m.nk if there are no variable ICs or controls, m.nk+m.nx*numel(con) if there are no variable controls, or scalar')
    end
elseif useModelICs && ~useModelInputs
    % Bounds can be nk+nx+nq, nT, nk, nk+nx, or 1
    l = numel(bounds);
    
    if l == 1
        bounds = zeros(nT,1) + bounds;
    elseif l == nk+nx+nq
        bounds = bounds([useParams; useICs; cat(1,useControls{:})]);
    elseif l == nT
        %bounds = bounds;
    elseif l == nk && nTx == 0 && nTq == 0
        bounds = bounds(useParams);
    elseif l == nk+nx && nTq == 0
        bounds = bounds([useParams; useICs]);
    else
        error('KroneckerBio:BoundSize', ...
            'LowerBound and UpperBound must be vectors the length of m.nk+m.nx+m.nq, number of variable parameters, m.nk if there are no variable ICs or controls, m.nk+m.nx if there are no variable controls, or scalar')
    end
else%~useModelICs && ~useModelInputs
    % Bounds can be nk+nx*nCon+nq, nT, nk, nk+nx*nCon, or 1
    l = numel(bounds);
    
    if l == 1
        bounds = zeros(nT,1) + bounds;
    elseif l == nk+nx*nCon+nq
        bounds = bounds([useParams; vec(useICs); cat(1,useControls{:})]);
    elseif l == nT
        %bounds = bounds;
    elseif l == nk && nTx == 0 && nTq == 0
        bounds = bounds(useParams);
    elseif l == nk+nx*nCon && nTq == 0
        bounds = bounds([useParams; vec(useICs)]);
    else
        error('KroneckerBio:BoundSize', ...
            'LowerBound and UpperBound must be vectors the length of m.nk+m.nx*numel(con)+m.nq, number of variable parameters, m.nk if there are no variable ICs or controls, m.nk+m.nx*numel(con) if there are no variable controls, or scalar')
    end
end