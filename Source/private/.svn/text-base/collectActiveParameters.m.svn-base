function T = collectActiveParameters(m, con, useModelICs, useModelInputs, useParams, useICs, useControls)

% Constants
nCon = size(con, 1);
nx   = m.nx;

% Store complete parameter sets
k = m.k;

if useModelICs
    x0 = m.x0;
else
    x0 = zeros(nx, nCon);
    for iCon = 1:nCon
        x0(:,iCon) = con(iCon).x0;
    end
end

if useModelInputs
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
T = [k(useParams); vec(x0(useICs)); q(cat(1, useControls{:}))];