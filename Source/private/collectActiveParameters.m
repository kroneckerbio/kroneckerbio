function T = collectActiveParameters(m, con, UseParams, UseSeeds, UseInputControls, UseDoseControls)

% Constants
nCon = size(con, 1);
ns   = m.ns;

% Store complete parameter sets
k = m.k;

s = zeros(ns, nCon);
for iCon = 1:nCon
    s(:,iCon) = con(iCon).s;
end

nq = sum(cat(1,con.nq));
q = zeros(nq,1);
index = 0;
for iCon = 1:nCon
    q(index+1:index+con(iCon).nq) = con(iCon).q;
    index = index+con(iCon).nq;
end

nh = sum(cat(1,con.nh));
h = zeros(nh,1);
index = 0;
for iCon = 1:nCon
    h(index+1:index+con(iCon).nh) = con(iCon).h;
    index = index+con(iCon).nh;
end

% Construct starting variable parameter set
T = [k(UseParams);
    vec(s(UseSeeds));
    q(cat(1, UseInputControls{:}));
    h(cat(1, UseDoseControls{:}))];
