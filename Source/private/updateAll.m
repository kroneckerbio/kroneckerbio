function [m, conNew, obj] = updateAll(m, con, T, UseModelSeeds, useModelInputs, useParams, UseSeeds, useControls)
%updateAll Update Kronecker Bio structures when parameters change
%
%   [m, con, obj] = updateAll(m, con, obj, T, UseModelSeeds, useModelInputs, useParams, UseSeeds, useControls)
%
%   This function performs the often needed task of updating m, con, and
%   obj whenever the active parameters T change.

% Constants
ns = m.ns;
nCon = size(con, 1);
nTk = nnz(useParams);
nTs = nnz(UseSeeds);
[useControls nTq] = fixUseControls(useControls, useModelInputs, nCon, m.nq, cat(1,con.nq));
nT = nTk + nTs + nTq;

% Update parameter sets
k = m.k;
k(useParams) = T(1:nTk);

if UseModelSeeds
    s = m.s;
    s(UseSeeds) = T(nTk+1:nTk+nTs);
else
    s = zeros(ns,nCon);
    for iCon = 1:nCon
        s(:,iCon) = con(iCon).s;
    end
    s(UseSeeds) = T(nTk+1:nTk+nTs);
end

if useModelInputs
    q = m.q;
    q(useControls) = T(nTk+nTs+1:nT);
else
    q = cell(nCon,1);
    endIndex = 0;
    for iCon = 1:nCon
        q{iCon} = con(iCon).q;
        startIndex = endIndex + 1;
        endIndex = endIndex + nnz(useControls{iCon});
        q{iCon}(useControls{iCon}) = T(nTk+nTs+startIndex:nTk+nTs+endIndex);
    end
end

% Update model
if UseModelSeeds && useModelInputs
    m = m.Update(k, s, q);
elseif UseModelSeeds
    m = m.Update(k, s, m.q);
elseif useModelInputs
    m = m.Update(k, m.s, q);
else
    m = m.Update(k, m.s, m.q);
end

% Update experimental conditions
if ~isnumeric(con)
    conNew = Uzero(nCon);
    if UseModelSeeds && useModelInputs
        for iCon = 1:nCon
            conNew(iCon) = con(iCon).Update(con(iCon).s, con(iCon).q);
        end
    elseif UseModelSeeds
        for iCon = 1:nCon
            conNew(iCon) = con(iCon).Update(con(iCon).s, q{iCon});
        end
    elseif useModelInputs
        for iCon = 1:nCon
            conNew(iCon) = con(iCon).Update(s(:,iCon), con(iCon).q);
        end
    else
        for iCon = 1:nCon
            conNew(iCon) = con(iCon).Update(s(:,iCon), q{iCon});
        end
    end
end
