function T = collectActiveParameters(m, con, UseModelSeeds, UseModelInputs, UseParams, UseSeeds, UseControls)

% Constants
nCon = size(con, 1);
ns   = m.ns;

% Store complete parameter sets
k = m.k;

if UseModelSeeds
    s = m.s;
else
    s = zeros(ns, nCon);
    for iCon = 1:nCon
        s(:,iCon) = con(iCon).s;
    end
end

if UseModelInputs
    q = m.q;
else
    nq = sum(cat(1,con.nq));
    q = zeros(nq,1);
    index = 0;
    for iCon = 1:nCon
        q(index+1:index+con(iCon).nq) = con(iCon).q;
    end
end

% Construct starting variable parameter set
T = [k(UseParams); vec(s(UseSeeds)); q(cat(1, UseControls{:}))];
