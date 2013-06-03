function [m, conNew, obj] = updateAll(m, con, obj, T, useModelICs, useModelInputs, useParams, useICs, useControls)
%UPDATEALL Update Kronecker Bio structures when parameters change
%
%   [m, con, obj] = updateAll(m, con, obj, T, useParams, useICs, useModelICs)
%
%   This function performs the often needed task of updating m, con, and
%   obj whenever the active parameters T change.

% Constants
nx = m.nx;
nCon = size(con, 1);
nObj = size(obj, 1);
nTk = nnz(useParams);
nTx = nnz(useICs);
[useControls nTq] = fixUseControls(useControls, useModelInputs, nCon, m.nq, cat(1,con.nq));
nT = nTk + nTx + nTq;

% Update parameter sets
k = m.k;
k(useParams) = T(1:nTk);
if useModelICs
    x0 = m.x0;
    x0(useICs) = T(nTk+1:nTk+nTx);
else
    x0 = zeros(nx,nCon);
    for iCon = 1:nCon
        x0(:,iCon) = con(iCon).x0;
    end
    x0(useICs) = T(nTk+1:nTk+nTx);
end
if useModelInputs
    q = m.q;
    q(useControls) = T(nTk+nTx+1:nT);
else
    q = cell(nCon,1);
    endIndex = 0;
    for iCon = 1:nCon
        q{iCon} = con(iCon).q;
        startIndex = endIndex + 1;
        endIndex = endIndex + nnz(useControls{iCon});
        q{iCon}(useControls{iCon}) = T(nTk+nTx+startIndex:nTk+nTx+endIndex);
    end
end

% Update model
if useModelICs && useModelInputs
    m = m.Update(k, x0, q);
elseif useModelICs
    m = m.Update(k, x0, m.q);
elseif useModelInputs
    m = m.Update(k, m.x0, q);
else
    m = m.Update(k, m.x0, m.q);
end

% Update experimental conditions
if ~isnumeric(con)
    conNew = Uzero(nCon);
    if useModelICs && useModelInputs
        for iCon = 1:nCon
            conNew(iCon) = pastestruct(Uzero(m), con(iCon).Update(con(iCon).x0, con(iCon).q));
        end
    elseif useModelICs
        for iCon = 1:nCon
            conNew(iCon) = pastestruct(Uzero(m), con(iCon).Update(con(iCon).x0, q{iCon}));
        end
    elseif useModelInputs
        for iCon = 1:nCon
            conNew(iCon) = pastestruct(Uzero(m), con(iCon).Update(x0(:,iCon), con(iCon).q));
        end
    else
        for iCon = 1:nCon
            conNew(iCon) = pastestruct(Uzero(m), con(iCon).Update(x0(:,iCon), q{iCon}));
        end
    end
end

% Update objective functions
if ~isnumeric(obj)
    % refeshObj expects vector m
    obj = refreshObj(m, con, obj, useParams, useICs, useControls);
end