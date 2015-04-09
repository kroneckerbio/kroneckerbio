function [m, conNew] = updateAll(m, con, T, UseParams, UseSeeds, UseInputControls, UseDoseControls)
%updateAll Update Kronecker Bio structures when parameters change
%
%   [m, con] = updateAll(m, con, obj, T, UseParams, UseSeeds, UseInputControls, UseDoseControls)
%
%   This function performs the often needed task of updating m, con, and
%   obj whenever the active parameters T change.

% Constants
ns = m.ns;
nCon = size(con, 1);
nTk = nnz(UseParams);
nTs = nnz(UseSeeds);
[UseInputControls, nTq] = fixUseControls(UseInputControls, nCon, cat(1,con.nq));
[UseDoseControls, nTh] = fixUseControls(UseDoseControls, nCon, cat(1,con.nh));
nT = nTk + nTs + nTq + nTh;

% Update parameter sets
k = m.k;
k(UseParams) = T(1:nTk);

s = zeros(ns,nCon);
for iCon = 1:nCon
    s(:,iCon) = con(iCon).s;
end
s(UseSeeds) = T(nTk+1:nTk+nTs);

q = cell(nCon,1);
endIndex = 0;
for iCon = 1:nCon
    q{iCon} = con(iCon).q;
    startIndex = endIndex + 1;
    endIndex = endIndex + nnz(UseInputControls{iCon});
    q{iCon}(UseInputControls{iCon}) = T(nTk+nTs+startIndex:nTk+nTs+endIndex);
end

h = cell(nCon,1);
endIndex = 0;
for iCon = 1:nCon
    h{iCon} = con(iCon).h;
    startIndex = endIndex + 1;
    endIndex = endIndex + nnz(UseDoseControls{iCon});
    h{iCon}(UseDoseControls{iCon}) = T(nTk+nTs+nTq+startIndex:nTk+nTs+nTq+endIndex);
end

% Update model
m = m.Update(k);

% Update experimental conditions
if ~isnumeric(con)
    conNew = experimentZero(nCon);
    for iCon = 1:nCon
        conNew(iCon) = con(iCon).Update(s(:,iCon), q{iCon}, h{iCon});
    end
end
