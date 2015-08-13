function T = collectActiveParameters(m, con, UseParams, UseSeeds, UseInputControls, UseDoseControls)

% Constants
n_con = size(con, 1);
ns = m.ns;

% Store complete parameter sets
k = m.k;

s = zeros(ns, n_con);
for i_con = 1:n_con
    s(:,i_con) = con(i_con).s;
end

nq = sum(cat(1,con.nq));
q = zeros(nq,1);
index = 0;
for i_con = 1:n_con
    q(index+1:index+con(i_con).nq) = con(i_con).q;
    index = index + con(i_con).nq;
end

nh = sum(cat(1,con.nh));
h = zeros(nh,1);
index = 0;
for i_con = 1:n_con
    h(index+1:index+con(i_con).nh) = con(i_con).h;
    index = index + con(i_con).nh;
end

% Construct starting variable parameter set
T = [k(UseParams);
    vec(s(UseSeeds));
    q(cat(1, UseInputControls{:}));
    h(cat(1, UseDoseControls{:}))];
